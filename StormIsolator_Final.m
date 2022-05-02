%% Set SHIELDOSE-2 Params and Cutoff
clear storm
flux_start_cutoff=1; %if 8 MeV flux is above this number, start recording this as a storm. you define what is considered a storm.
flux_end_cutoff=1;
energy_cutoff=3; %1=6 MeV, 2=8MeV, 3=12MeV, 4=18MeV and so on, matching Egrid


onera_desp_lib_load('onera_desp_lib_Win64_x86.dll')%load this for SHIELDOSE-2

Egrid=[6.01 8.7 12.58 18.18 26.3 38.03 54.99 79.53 115 166.3 244.2]; %define that these are the energy ranges of GOES data. makes SHIELDOSE-2 inputting easier later

Target = struct('material','Si','unit','mils');
options = {'INUC',3,'NPTS',1001,'perYear',0};
Target.depth =  [2:2:160];

ProtSpect.E=[];%leave these all as 0 because we don't care about these. only solar protons.
ProtSpect.Flux=[];
ElecSpect.E=[];
ElecSpect.Flux=[];

% %Egrid2=logspace(log10(.1),log10(500),100);
% LowEnd=logspace(log10(.1),log10(10),50);
% HighEnd=logspace(log10(10),log10(1000),90);
% Egrid2=horzcat(LowEnd,HighEnd(2:end));%energy bins you want, usually want to keep log spaced

%load('D:\Materials\All_SEP_Events_noGeoShield.mat')
model(4:end)=[];
Egrid2=model(1).Egrid;
%Egrid2=Egrid;
%% Load and configure GOES data
load('D:\Materials\NOAADates.mat');
%load('D:\Materials\All_SEP_Events_noGeoShield');
load('D:\Materials\H_He_SPE_Fluxes_noGeoShield.mat');
TenMeVIndex=find(Egrid2==min(sort(abs(Egrid2-10)))+10)

SolSpect.E=Egrid;%these are always the same so just do it here so it's not rewritten 4.5 million times during the loop
SolSpect.Erange=[min(Egrid2),max(Egrid2)];

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
    consecutive_counter=0;%because we need the threshold to be exceeded for at least 3 periods
    finder_counter=0;
    while sum(SPE_Flux_H_Array(date_indices(i)+finder_counter,energy_cutoff:end))<flux_start_cutoff==1 && (date_indices(i)+finder_counter)<height(SPE_Flux_H) || consecutive_counter<3%since the SPE dataset and GOES dataset aren't synced we have to find the storm
        finder_counter=finder_counter+1; %we do this by starting at the SPE dataset start date (typically a few hours behind the GOES data)
        if date_indices(i)+finder_counter<height(SPE_Flux_H)%count up until the end of dataset
            if sum(SPE_Flux_H_Array(date_indices(i)+finder_counter,energy_cutoff:end))<flux_start_cutoff==0 %if we exceed the threshold
                consecutive_counter=consecutive_counter+1;%we start counting how many times in a row
            else
                consecutive_counter=0;%if we fail to exceed then we restart the counter 
            end
        end
        
        %         if i==258 && finder_counter==18362
%             break
%         end
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
    storm(i).GOES_Start=string(SPE_Flux_H.Time(date_indices(i)+finder_counter,1));
    storm(i).NOAA_Start=NOAADates(i);
    
    tracker_counter=0;
    while sum(SPE_Flux_H_Array(date_indices(i)+finder_counter+tracker_counter+25,energy_cutoff:end))>flux_start_cutoff==1 %while the storm is above 10 PFU we keep recording it's data as an SPE
        storm(i).flux=vertcat(storm(i).flux,SPE_Flux_H_Array(date_indices(i)+finder_counter+tracker_counter,1:end));
        tracker_counter=tracker_counter+1;
    end
    
    storm(i).days=height(storm(i).flux)*5/60/24; %duration in days
    storm(i).hours=storm(i).days*24; %duration in hours
    storm(i).seconds=storm(i).hours*3600; %duration in seconds
    storm(i).Datenum=datenum(storm(i).Start);
end
for i=1:length(storm) %storms that are not behind need the opposite process run, go backwards to find the start date and forward to find the end
    
    if isempty(storm(i).flux)==1;
        finder_counter=0;
        while sum(SPE_Flux_H_Array(date_indices(i)-finder_counter,energy_cutoff:end))>flux_start_cutoff==1 && (date_indices(i)+finder_counter)<height(SPE_Flux_H)%go backward until we're below 10 PFU
            finder_counter=finder_counter+1;
        end
        storm(i).Start=SPE_Flux_H.Time(date_indices(i)-finder_counter);%when it exits the loop is when the storm actually starts
        storm(i).StartString=string(storm(i).Start)
        
        tracker_counter=0;
        while sum(SPE_Flux_H_Array(date_indices(i)-finder_counter+tracker_counter,energy_cutoff:end))>flux_start_cutoff==1  %then we record normally as with before: while the storm is above 10 PFU we keep recording it's data as an SPE
            storm(i).flux=vertcat(storm(i).flux,SPE_Flux_H_Array(date_indices(i)-finder_counter+tracker_counter,1:end));
            tracker_counter=tracker_counter+1;
        end
        storm(i).days=height(storm(i).flux)*5/60/24; %duration in days
        storm(i).hours=storm(i).days*24; %duration in hours
        storm(i).seconds=storm(i).hours*3600; %duration in seconds
    end
    storm(i).GOES_Flux=storm(i).flux; %want to save the old fluxes in case we want to compare later
end

%% Extrapolate the Storms

for i=1:length(storm)%extrapolate energy spectra from GOES to 75 points
    placeholder=zeros(height(storm(i).flux),length(Egrid2));%can't combine arrays with different sizes (extrapolated w/non-extrapolated) so we use a placeholder
    for ii=1:height(storm(i).flux) %for each 5 min spectra
        %Linear fit for E<10
        corrected_flux=storm(i).flux(ii,:); %get the spectra
        corrected_flux(corrected_flux==0)=[];%so we can correct it by removing 0's that can't work in log-log space
        Const = polyfit(log10(Egrid(1:length(corrected_flux))),log10(corrected_flux),1);%fit the log-log space curve
        m = Const(1);
        k = Const(2);
        extrapolated_spectra_linear = 10.^(m.*log10(Egrid2(1:TenMeVIndex-1))+(k));%extrapolate the fit to the 75 points

        Const = fit(Egrid',storm(i).flux(ii,:)','power1');%fit the log-log space curve
        CL=confint(Const);%Get confidence level of fitted parameters
        Lower_Bound=CL(1,2);
        a=Const.a;
        b=(Const.b+Lower_Bound)/2;%Set b to average of fitted parameter and 95th CI lower bound, just fits well
        extrapolated_spectra_power = a.*(Egrid2(TenMeVIndex:end).^b);%extrapolate the fit to the 75 points
        
        extrapolated_spectra=vertcat(extrapolated_spectra_linear,extrapolated_spectra_power);
        
    placeholder(ii,:)=extrapolated_spectra; %build out the placeholder
    end
    bad_vals=isnan(placeholder(:,1));%if the fit runs into any snags we just remove them
    placeholder(bad_vals,:)=[];
    storm(i).flux=placeholder;%replace the storm's time-flux profile with the extrapolated one in the placeholder now that everything is correct dimensions
      
    [maximumFlux idx]=max(storm(i).GOES_Flux(:,2));
    Const=fit(Egrid',storm(i).GOES_Flux(idx,:)','power1');
    CL=confint(Const);%Get confidence level of fitted parameters
    Lower_Bound=CL(1,2);
    a=Const.a;
    b=(Const.b+Lower_Bound)/2;
    fun = @(E) (a.*(E.^b));%integral flux will give us PFU    
    storm(i).MaxPFU_Estimation=integral(fun,10,1000);
end

counter=0;
big_storm=0;
for i=1:length(storm)
    storm(i).diffPFU=abs(storm(i).unnamed-storm(i).MaxPFU_Estimation)
    if abs(storm(i).diffPFU)>storm(i).unnamed*.2 && storm(i).unnamed>1000
        counter=counter+1;
    end
    if storm(i).unnamed>1000
    big_storm_NOAA=big_storm+1;
    end
end

bad_array=[];
wondering=[];
big=[];
whatif=[];
howmany=0;
for i=1:length(storm)
    if storm(i).diffPFU>1000 && abs(storm(i).diffPFU)>storm(i).unnamed*.1
        bad_array=horzcat(bad_array,i);
        wondering=horzcat(wondering,height(storm(i).flux));
        big=horzcat(big,storm(i).unnamed);
        whatif=horzcat(whatif,storm(i).MaxPFU_Estimation);
    end
    if storm(i).diffPFU>1000
        howmany=howmany+1;
    end
end

v
% end

% figure()
% plot(log10(Egrid(1:length(corrected_flux))),log10(corrected_flux))
% hold on
% plot(log10(Egrid2),log10(extrapolated_spectra))
% 
% A=10^abs((-.36045/(m+1))*(10000*(((10000/10)^m)-1)))

% storm(1).intFlux
% storm(88).intFlux
% storm(5).intFlux
% save('Storms_Unpropagated.mat','storm')