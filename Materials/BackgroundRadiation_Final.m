%% Read in AE9AP9 Files and Interpolate Dose Rates at Each Shielding Level
addpath 'D:\Materials'

minFile='SkySatSolarMin.AE9.output_mean_TotaldoserateFullAvg.txt'; %Name of File at Solar Minimum (the higher dose rate)
maxFile='SkySatSolarMax.AE9.output_mean_TotaldoserateFullAvg.txt'; %Name of File at Solar Maximum (the lower dose rate)

MinFileID=fopen(minFile);
MaxFileID=fopen(maxFile);

g = textscan(MinFileID,'%s','delimiter','\n','HeaderLines',18);%skip everything until the shielding depths

Shielding=g{1,1}(1); %Extract the Shielding Depths
Shielding = cellfun(@(x) strsplit(x, ' '), Shielding, 'UniformOutput', false);
Shielding=str2double(Shielding{1,1}(5:end-1));%Convert it to an array

SolarMinDoseRates=g{1,1}(10);%Extract the minimum dose rates
SolarMinDoseRates = cellfun(@(x) strsplit(x, ','), SolarMinDoseRates, 'UniformOutput', false);
SolarMinDoseRates=str2double(SolarMinDoseRates{1,1}(5:end));
SolarMinExtrapolatedDoseRates=interp1(Shielding,SolarMinDoseRates,[4:1:400]);%interpolate these dose rates to all shielding levels


g = textscan(MaxFileID,'%s','delimiter','\n','HeaderLines',18);%Same as above
SolarMaxDoseRates=g{1,1}(10);
SolarMaxDoseRates = cellfun(@(x) strsplit(x, ','), SolarMaxDoseRates, 'UniformOutput', false);
SolarMaxDoseRates=str2double(SolarMaxDoseRates{1,1}(5:end));
SolarMaxExtrapolatedDoseRates=interp1(Shielding,SolarMaxDoseRates,[4:1:400]);

BackgroundSolarCycleDoses=((SolarMinExtrapolatedDoseRates(:)-SolarMaxExtrapolatedDoseRates(:)).*abs(cos([1:1:133]*pi/133))+SolarMaxExtrapolatedDoseRates(84))*60*60*24*365/30/1000; %converts to dose per monthly krad

save('BackgroundSolarCycleDoses_SkySat.mat','BackgroundSolarCycleDoses')