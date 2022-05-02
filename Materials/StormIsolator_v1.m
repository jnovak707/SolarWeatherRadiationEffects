%% Set SHIELDOSE-2 Params and Cutoff
clear all
flux_start_cutoff=1; %if 8 MeV flux is above this number, start recording this as a storm. you define what is considered a storm.
flux_end_cutoff=1;
energy_cutoff=2; %1=6 MeV, 2=8MeV, 3=12MeV, 4=18MeV and so on, matching Egrid


onera_desp_lib_load('onera_desp_lib_Win64_x86.dll')%load this for SHIELDOSE-2

Egrid=[6 8 12 18 26 38 55 80 115 166 244]; %define that these are the energy ranges of GOES data. makes SHIELDOSE-2 inputting easier later

Target = struct('material','Si','unit','mils');
options = {'INUC',3,'NPTS',1001,'perYear',0};
Target.depth =  [2:2:160];

ProtSpect.E=[];%leave these all as 0 because we don't care about these. only solar protons.
ProtSpect.Flux=[];
ElecSpect.E=[];
ElecSpect.Flux=[];

%% Load and configure GOES data
load('D:\NOAADates');
load('D:\All_SEP_Events_noGeoShield');
load('D:\H_He_SPE_Fluxes_noGeoShield');
load('All_SEP_Events_noGeoShield.mat')
model(4:end)=[];
Egrid2=model(1).Egrid;

SolSpect.E=Egrid;%these are always the same so just do it here so it's not rewritten 4.5 million times during the loop
SolSpect.Erange=[.001,1e5];

SPE_Flux_H=renamevars(SPE_Flux_H,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12"],["Time","MeV6","MeV8","MeV12","MeV18","MeV26","MeV38","MeV55","MeV80","MeV115","MeV166","MeV244"]); %set headers to energy levels for easier referencing

storm.flux=[];

storm_counter=0;

storm_recording=0;%flag for if we're recording a storm or not

%% Capturing all the storms and their characteristics according to our definition

SPE_Flux_H_Array=table2array(SPE_Flux_H(:,2:end));%convert it to an array to make indexing and pulling values easier
tic
date_indices=[];
for i=1:length(NOAADates)
    date_indices=horzcat(date_indices,find(SPE_Flux_H.Time==NOAADates(i)));
end

storm_recording=0;
for i=1:length(NOAADates)
        
    finder_counter=0;
        while sum(SPE_Flux_H_Array(date_indices(i)+finder_counter,energy_cutoff:end))<flux_start_cutoff==1 %since the SPE dataset and GOES dataset aren't synced we have to find the storm
            finder_counter=finder_counter+1; %we do this by starting at the SPE dataset start date (typically a few hours behind the GOES data)
            if i==258 && finder_counter==18362
                break
            end
        end %and counting up until we get to the 10 PFU that indicates SPE start

    if i==259 && finder_counter==18362
        break
    end
    storm(i).Start=SPE_Flux_H.Time(date_indices(i)+finder_counter);
    if i>=2
    if storm(i).Start==storm(i-1).Start
        fprintf("%.f Not Found\n",i);
    end
    end
    storm(i).StartString=string(storm(i).Start);
    
    tracker_counter=0;
    while sum(SPE_Flux_H_Array(date_indices(i)+finder_counter+tracker_counter,energy_cutoff:end))>=flux_start_cutoff==1 %while the storm is above 10 PFU we keep recording it's data as an SPE
        storm(i).flux=vertcat(storm(i).flux,SPE_Flux_H_Array(date_indices(i)+finder_counter+tracker_counter,1:end));
        tracker_counter=tracker_counter+1;
    end
    storm(i).days=height(storm(i).flux)*5/60/24; %duration in days
    storm(i).hours=storm(i).days*24; %duration in hours
    storm(i).seconds=storm(i).hours*3600; %duration in seconds
end 

for i=1:length(storm) %storms that are not behind need the opposite process run, go backwards to find the start date and forward to find the end
    
    if isempty(storm(i).flux)==1;    
    finder_counter=0;
    while sum(SPE_Flux_H_Array(date_indices(i)-finder_counter,energy_cutoff:end))>flux_start_cutoff==1 %go backward until we're below 10 PFU
       finder_counter=finder_counter+1;
    end 
    storm(i).Start=SPE_Flux_H.Time(date_indices(i)-finder_counter);%when it exits the loop is when the storm actually starts
    storm(i).StartString=string(storm(i).Start);
    
    tracker_counter=0;
    while sum(SPE_Flux_H_Array(date_indices(i)-finder_counter+tracker_counter,energy_cutoff:end))>flux_start_cutoff==1 %then we record normally as with before: while the storm is above 10 PFU we keep recording it's data as an SPE
        storm(i).flux=vertcat(storm(i).flux,SPE_Flux_H_Array(date_indices(i)-finder_counter+tracker_counter,1:end));
        tracker_counter=tracker_counter+1;
    end
    storm(i).days=height(storm(i).flux)*5/60/24; %duration in days
    storm(i).hours=storm(i).days*24; %duration in hours
    storm(i).seconds=storm(i).hours*3600; %duration in seconds
    end
end

%% Extrapolate the Storms

for i=1:length(storm)%extrapolate energy spectra from GOES to 75 points
    placeholder=zeros(height(storm(i).flux),75);%can't combine arrays with different sizes (extrapolated w/non-extrapolated) so we use a placeholder
    for ii=1:height(storm(i).flux) %for each 5 min spectra
    corrected_flux=storm(i).flux(ii,:); %get the spectra
    corrected_flux(corrected_flux==0)=[];%so we can correct it by removing 0's that can't work in log-log space
    Const = polyfit(log10(Egrid(1:length(corrected_flux))),log10(corrected_flux),1);%fit the log-log space curve
        m = Const(1);
        k = Const(2);
        extrapolated_spectra = 10.^(m.*log10(Egrid2)+(k));%extrapolate the fit to the 75 points
    placeholder(ii,:)=extrapolated_spectra; %build out the placeholder
    end
    bad_vals=isnan(placeholder(:,1));%if the fit runs into any snags we just remove them
    placeholder(bad_vals,:)=[];
    storm(i).flux=placeholder;%replace the storm's time-flux profile with the extrapolated one in the placeholder now that everything is correct dimensions
end