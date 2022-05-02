%% Total Ionizing Dose Intro
%So far the flow is: 
% 1. Storm Isolator_Final - Retrieve NOAA storms from GOES data
% 2. Geomagnetic Shielding_Final - Apply geomagnetic shielding to storms

% And then there are a few .m and .mat files associated with the exploratory
% analysis such as:
% 1. Solar Cycle Comparison Plots - graphing solar cycles
% 2. SEPEM_events - dataset of SEPEM reference event list
% 3. NOAA_data - dataset of NOAA reference 
% 4. NOAADates - cleaned extraction of SPE event start dates from SPE site
% 5. H_He_SPE_Fluxes_noGeoShield - GOES 5 minute diffFlux dataset

%This code will use the shielded geomagnetic fluxes from Geomagnetic
%Shielding_Final to calculate the dose rates and total SPE dose from
%various SPEs at various shielding depths

%Note: The "storm" dataset is composed of SPEs which aren't unique to an
%individual satellite. The application of the geomagnetic shielding
%code is individualized to a satellite's location so that must be run first
%for the desired satellite's TLE and then this code follows.

%% Total Ionizing Dose Calculations
%% Reading in Files
load('D:\Propagated_Storms_SkySatShielded.mat');

onera_desp_lib_load('onera_desp_lib_Win64_x86.dll'); %load SHIELDOSE-2 files
Target = struct('material','Si','unit','mils');      %Silicon target, depth in mils
options = {'INUC',1,'NPTS',1001,'perYear',0};        %Pad the depths, result in rad/s
Target.depth =  horzcat(3.937, 5:5:100, 150, 200, 300, 400);   %Set the desired depths
depths=Target.depth; save('depths.mat','depths');
ProtSpect.E=[];     %For SPEs we only care about solar protons so these can all
ProtSpect.Flux=[];  %start and stay at 0
ElecSpect.E=[];
ElecSpect.Flux=[];

load('D:\Materials\All_SEP_Events_noGeoShield.mat')
model(4:end)=[];
Egrid2=model(1).Egrid;
SolSpect.E=Egrid2;  %This also never changes so we can leave it initialized only once

storm(1).TIDTimeProfile=[];
%% Total Ionizing Dose Calculations
tic
for i=1:length(storm)                       %for each storm
    for ii=1:width(storm(i).ShieldedFluxes) %for each 5 min diffFlux period
        SolSpect.Flux=storm(i).ShieldedFluxes(:,ii); %take that diffFlux and
        [ProtDose,ElecDose,BremDose,SolDose,TotDose] = onera_desp_lib_shieldose2(ProtSpect,ElecSpect,SolSpect,Target,options{:}); %Put it into SHIELDOSE-2
        storm(i).TIDTimeProfile=vertcat(storm(i).TIDTimeProfile,TotDose(:,3)'); %get the rad/s for each 5 minute interval for each shielding depth
    end
    if height(storm(i).TIDTimeProfile)==1
    storm(i).stormKrad=storm(i).TIDTimeProfile*.3;
    else
    storm(i).stormKrad=sum(storm(i).TIDTimeProfile)*.3; %sum rad/s for each 5 minute interval then multiple by .3 to get total rads (300 sec per 5min, 1000 rads per krad)
    end
    storm(i).maxKrad=storm(i).stormKrad(1,1);%krad per storm at minimum defined shielding
    i %Tells us which storm it just finished calculating
end %at the end we have TID/Time profiles and total TID for each storm for each shielding depth
toc

Peak_Flux=[];
for i=1:length(storm)
   Peak_Flux=horzcat(Peak_Flux,max(storm(i).flux(:,41)));
end

save('Storm_Doses_SkySatShielded.mat','storm')