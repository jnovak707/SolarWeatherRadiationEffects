onera_desp_lib_load('onera_desp_lib_Win64_x86.dll')
% shieldose2_test - test shieldose2 function

%     1     1    70     2
Target = struct('material','Si','unit','mils');
options = {'INUC',1,'NPTS',1001,'perYear',1};
Target.depth =  [4:2:250];

%  0.100  10000.000  0.100  10000.000  1001  0.050  10.000  1001
% GEOSTAT,35790 KM,INCL=0,PLONG=160W; SP:1AL(95%),TP:NONE,EL:AEI7-HI(79)
%     3     0    28 1.00000E+03 3.15360E+07
%SolSpect = struct('Erange',[0.1 1e4],'E0',2.65e+1,'N0',2.45e+10,'form','E');
% 0.0000E+00  0.0000E+00  0.0000E+00
%   2.6500E+01  2.4500E+10  0.0000E+00]';
SolSpect.N0 = SolSpect.N0/3.15360E+07; % convert to rate from fluence/duration
SolSpect.N0 = SolSpect.N0*1e3; % convert /keV to /MeV (don't do this for exponential spectrum)
ProtSpect = [];
SolSpect.E=[6.012736
8.695174
12.574353
18.184072
26.296427
38.027969
54.993218
79.527081
115.006277
166.313559
240.510291];

SolSpect.Flux = [
498740000
249700000
114910000
41950000
15943000
6258200
2180800
533620
155430
60188
23542
]';
SolSpect.Erange=[.05,500];
%ElecSpect(:,2) = ElecSpect(:,2)*1e3; % convert from /keV to /MeV
%ElecSpect = struct('E',ElecSpect(:,1),'Flux',ElecSpect(:,2),'Erange',[0.05 10]);

figure()
[ProtDose,ElecDose,BremDose,SolDose,TotDose] = onera_desp_lib_shieldose2(ProtSpect,ElecSpect,SolSpect,Target,options{:});
loglog(Target.depth,TotDose,Target.depth,ProtDose+ElecDose+BremDose+SolDose,'k.');
xlabel(sprintf('Depth %s (%s)',Target.material,Target.unit));
ylabel('Dose (rads/year)');
legend('DOSE IN SEMI-INFINITE ALUMINUM MEDIUM',...
    'DOSE AT TRANSMISSION SURFACE OF FINITE ALUMINUM SLAB SHIELDS',...
    '1/2 DOSE AT CENTER OF ALUMINUM SPHERES',...
    'location','northoutside');


