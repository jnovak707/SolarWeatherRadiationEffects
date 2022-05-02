%% Set SHIELDOSE-2 Params and Cutoff
fprintf("Loading Parameters...\n")
clear storm
flux_start_cutoff=1; %if 8 MeV flux is above this number, start recording this as a storm. you define what is considered a storm.
flux_end_cutoff=1;
energy_cutoff=2; %1=6 MeV, 2=8MeV, 3=12MeV, 4=18MeV and so on, matching Egrid

%% Load and configure GOES data
load('D:\Materials\NOAADates.mat');
%load('D:\Materials\All_SEP_Events_noGeoShield');
load('D:\Materials\H_He_SPE_Fluxes_noGeoShield.mat');
load('D:\Materials\NOAA_PFU.mat');
TenMeVIndex=find(Egrid2==min(sort(abs(Egrid2-30)))+30)

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

fprintf("Isolating Storms...\n")
for i=1:length(NOAADates)-1
    idx1=date_indices(i);
    idx2=date_indices(i+1);
    
    peak_idx=idx1+find(SPE_Flux_H_Array(idx1:idx2,energy_cutoff)==max(SPE_Flux_H_Array(idx1:idx2,energy_cutoff)));
    peak_idx=peak_idx(1);
    
    storm(i).peakTime=string(SPE_Flux_H.Time(peak_idx,1));
    
    find_start=0;%count backward until no longer above threshold to find start
    while sum(SPE_Flux_H_Array(peak_idx-find_start,energy_cutoff:end))>flux_start_cutoff==1 && (peak_idx-find_start)>idx1 
        storm(i).flux=vertcat(SPE_Flux_H_Array(peak_idx-find_start,:),storm(i).flux)
        find_start=find_start+1;
    end
    storm(i).GOES_Start=SPE_Flux_H.Time(peak_idx-find_start,1);
    storm(i).GOES_StartString=string(SPE_Flux_H.Time(peak_idx-find_start,1));
    storm(i).NOAA_Start=NOAADates(i);
    storm(i).NOAA_StartString=string(NOAADates(i));
    
    find_end=0;%count forward until no longer above threshold to find end
    while sum(SPE_Flux_H_Array(peak_idx-find_end,energy_cutoff:end))>flux_start_cutoff==1 && (peak_idx+find_end)<idx2 
        storm(i).flux=vertcat(storm(i).flux,SPE_Flux_H_Array(peak_idx+find_end,:))
        find_end=find_end+1;
    end
    
    storm(i).GOES_End=SPE_Flux_H.Time(peak_idx+find_end,1);
    storm(i).GOES_EndString=string(SPE_Flux_H.Time(peak_idx+find_end,1));
    storm(i).DurationHours=(find_start+find_end)*5/60;
    storm(i).GOES_Flux=storm(i).flux;
end

fprintf("Extrapolating...\n")
for i=1:length(storm)%extrapolate energy spectra from GOES to 75 points
    placeholder=zeros(height(storm(i).flux),length(Egrid2));%can't combine arrays with different sizes (extrapolated w/non-extrapolated) so we use a placeholder
    for ii=1:height(storm(i).flux) %for each 5 min spectra
        %Linear fit for E<10
        corrected_flux=storm(i).flux(ii,:); %get the spectra
        corrected_flux(corrected_flux==0)=[];%so we can correct it by removing 0's that can't work in log-log space
        
        if length(corrected_flux)>2
        Const = polyfit(log10(Egrid(storm(i).flux(ii,:)~=0)),log10(corrected_flux),1);%fit the log-log space curve
        m = Const(1);
        k = Const(2);
        extrapolated_spectra_linear = 10.^(m.*log10(Egrid2(1:TenMeVIndex-1))+(k));%extrapolate the fit to the 75 points
        else
            extrapolated_spectra_linear = 10.^(m.*log10(Egrid2(1:TenMeVIndex-1))+(k));
        end
        
%         check_zero=length(storm(i).flux(ii,:));
%         if isempty(find(storm(i).flux(ii,:)==0))~=1 && length(find(storm(i).flux(ii,:)==0))<=length(storm(i).flux(ii,:))-2 %check if multiple zero values in the spectra, second part is to make sure we have at least 2 points to fit a curve to
%             check_zero=find(storm(i).flux(ii,:)==0);%sometimes the data has zeros for a few energy levels and then suddenly jump to a small number.
%             check_zero=check_zero(1);%This ruins the fit so we remove every 0 after the first
%         end

        corrected_flux=storm(i).flux(ii,:); %get the spectra
        corrected_flux(corrected_flux==0)=[];%so we can correct it by removing 0's that throw off a power curve fit
        
        if length(corrected_flux)>2 %if we have the 2 data points needed to make the fit, we do it     
        
        Const = fit(Egrid(storm(i).flux(ii,:)~=0)',corrected_flux','power1');%fit the power curve in linear space
        CL=confint(Const);%Get confidence level of fitted parameters
        Lower_Bound=CL(1,2);
        a=Const.a;
        b=(Const.b+Lower_Bound)/2;%Set b to average of fitted parameter and 95th CI lower bound, just fits well
        extrapolated_spectra_power = a.*(Egrid2(TenMeVIndex:end).^b);%extrapolate the fit to the 75 points
        
        else %otherwise we just use the previous fit for this one
            extrapolated_spectra_power = a.*(Egrid2(TenMeVIndex:end).^b);
        end
            
        
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
    storm(i).NOAA_PFU=NOAA_PFU_Array(i);
    fprintf("Finished storm number: %.f\n",i);
end
toc

save('Storms_Unpropagated.mat','storm')