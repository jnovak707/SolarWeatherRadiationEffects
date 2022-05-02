%% Read in SEPEM Spectra Data
clc;
clear all;

addpath 'D:\'
addpath 'C:\Users\jonat\OneDrive\Documents\MATLAB\ThesisCode\ThesisDocs\space-environment\space-environment-master\Ae9Ap9\win64\bin\FlockThesis\FlockRun1'

fileID=fopen('FlockThesisRun.AE9.output_mean_TotaldoserateFullAvg.txt'); %replace with your AE9AP9 Dose Rate File, repeat for each confidence level if neccessary
g = textscan(fileID,'%s','delimiter','\n','HeaderLines',18);%skip everything until the shielding depths

Shielding=g{1,1}(1);                                                            %Extract the Shielding Depths
Shielding = cellfun(@(x) strsplit(x, ' '), Shielding, 'UniformOutput', false);
Shielding=str2double(Shielding{1,1}(5:end-1));

DoseRates=g{1,1}(10);
DoseRates = cellfun(@(x) strsplit(x, ','), DoseRates, 'UniformOutput', false);
DoseRates=str2double(DoseRates{1,1}(5:end));

ExtrapolatedDoseRates=interp1(Shielding,DoseRates,[4:1:400]);

% figure() %Validates that the interoplation looks good for the shielding levels chosen
% plot(Shielding,DoseRates)
% hold on
% plot([4:1:400],interp1(Shielding,DoseRates,[4:1:400]),'-r')
% legend('AE9AP9','Interpolation')

