%% Initialize STK
if ~exist('scenObj')
    try
        uiapp = actxGetRunningServer('STK12.application'); % Grab an existing instance of STK
    catch
        uiapp = actxserver('STK12.application'); % STK is not running, launch new instance, Launch a new instance of STK and grab it
    end
    %get the root from the personality
    %it has two... get the second, its the newer STK Object Model Interface as
    %documented in the STK Help
    root = uiapp.Personality2;
    uiapp.visible = 1; % set visible to true (show STK GUI)
    try %close current scenario or open new one
        root.CloseScenario(); root.NewScenario('Flock_Scenario');
    catch
        root.NewScenario('Flock_Scenario');
    end
    scenObj = root.CurrentScenario; %get the scenario root, its of type IAgScenario
end

%% Initialize Parameters

%Setup time interval
params.epoch_start = datenum('01-Apr-2014'); %date number (for STK)
params.sim_duration = 7; %days
params.timestep = 60; %sec

%Setup sensor parameters
params.sens.vert_FOV_half_angle = 1.6; %deg
params.sens.horiz_FOV_half_angle = 1.6; %deg
params.sens.target_lighting = 'DirectSun'%DirectSun, PenumbraDirectSun
params.sens.LOS_sun_excl = 60; %deg
params.sens.A_FOV_pointing_azel = [0,90,params.sens.vert_FOV_half_angle,params.sens.horiz_FOV_half_angle];

%Setup constellation parameters
params.const.num_planes = 8;
params.const.num_sats_per_plane = 40;
params.const.F = 1;
params.const.inc = 97.5; %deg
params.const.alt = 450; %km

params.cov.granularity = 10; %deg

% constants
params.const.Re = 6378.1366; %km [SMAD-SME (IAU)]
params.const.mu = 398600.43560; %km3/sec2 [SMAD-SME (IAU)]
params.const.J2 = 0.0010826269;
params.const.Ts = 365.25636; %days [SMAD-SME (astronomical almanac)]
params.const.Tc = 365; %days [Calendar Year]
params.const.T_sideral_rotation = 86164.09053; %sec [SMAD-SME (astronomical almanac)]
params.const.w_earth = 2*pi/params.const.T_sideral_rotation; %rad/sec (inertial rotation rate)
params.const.g = 9.80665; %m/s/s [SMAD-SME]

coverage_run_type = 'static_point';
%% Set Simulation Time and Date Preferences

root.UnitPreferences.Item('DistanceUnit').SetCurrentUnit('km');
root.UnitPreferences.Item('TimeUnit').SetCurrentUnit('sec');
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('UTCG');
epoch = datestr((params.epoch_start), 'dd mmm yyyy HH:MM:SS.FFF');
scenObj.Epoch = epoch;
scenObj.StartTime = scenObj.Epoch;
scenObj.StopTime = datestr((params.epoch_start+params.sim_duration), 'dd mmm yyyy HH:MM:SS.FFF');

%% Setup Coverage Definition

% Create Coverage Definitions
cov_def = scenObj.Children.New('eCoverageDefinition', 'Space_Cov_Def');

% set up calculation grids
cov_def.Grid.BoundsType = 'eBoundsLatLonRegion';
cov_def.Grid.Bounds.MinLongitude = -180;
cov_def.Grid.Bounds.MaxLongitude = -90;
cov_def.Grid.Bounds.MinLatitude = -60;
cov_def.Grid.Bounds.MaxLatitude = 60;
cov_def.Grid.Resolution.LatLon = 3;

%setup FOM
% fom = cov_def.Children.New('eFigureOfMerit', 'average_revisits');
% fom.SetDefinitionType('eFmRevisitTime');
% fom.Definition.SetComputeType('eAverage');
%% Create Walker Constellation
walker_Np = params.const.num_planes;
walker_Nso = params.const.num_sats_per_plane;
F = params.const.F;

orbital_state = zeros(walker_Np*walker_Nso,6); %initialize satellite state tracker
count = 0;
init_RAAN = deg2rad(3.4467);
init_inc = params.const.inc;
RAANvect = zeros(length(orbital_state(:,1)),1);
incVect = zeros(length(orbital_state(:,1)),1);
Mvect = zeros(length(orbital_state(:,1)),1);
for k = 0:walker_Np-1
    for l = 0:walker_Nso-1
        count = count+1;
        RAANvect(count,1) = deg2rad(1*k)+init_RAAN;
        incVect(count,1) = (.5*k)+init_inc;
        Mvect(count,1) = 2*pi*l/walker_Nso+2*pi*k*F/(walker_Nso*walker_Np);
    end
end

orbital_state(:,1) = ones(length(orbital_state(:,1)),1)*(params.const.alt+params.const.Re); %semimajor axis (km)
orbital_state(:,2) = zeros(length(orbital_state(:,1)),1); %eccentricity
orbital_state(:,3) = incVect; %inclination (deg)
orbital_state(:,4) = RAANvect; %right ascention of the ascending node (rad)
orbital_state(:,5) = zeros(length(orbital_state(:,1)),1); %argument of perigee (rad)
orbital_state(:,6) = Mvect; %mean anomaly (rad)

%% Create & Assign Satellites
numPlane=params.const.num_planes;
satPlane=params.const.num_sats_per_plane;
root.BeginUpdate()
i=1;
for i=1:(numPlane*satPlane)
    satellite{i} = scenObj.Children.New('eSatellite', ['Satellite_',num2str(i)]);
    root.ExecuteCommand(['Graphics */Satellite/Satellite_',num2str(i),' Basic Color white']);
    root.ExecuteCommand(['Graphics */Satellite/Satellite_',num2str(i),' Pass2D GrndLead None GrndTrail None OrbitLead None OrbitTrail None']);
    root.ExecuteCommand(['Graphics */Satellite/Satellite_',num2str(i),' Label Off']);
    %create sensor A
    root.ExecuteCommand(['New / */Satellite/Satellite_',num2str(i),'/Sensor sens_A']);
    root.ExecuteCommand(['Define */Satellite/Satellite_',num2str(i),'/Sensor/sens_A Rectangular ',num2str(params.sens.A_FOV_pointing_azel(3)),' ',num2str(params.sens.A_FOV_pointing_azel(4))]);
    root.ExecuteCommand(['Point */Satellite/Satellite_',num2str(i),'/Sensor/sens_A Fixed AzEl ',num2str(params.sens.A_FOV_pointing_azel(1)),' ',num2str(params.sens.A_FOV_pointing_azel(2)),' Hold']);
    root.ExecuteCommand(['SetConstraint */Satellite/Satellite_',num2str(i),'/Sensor/sens_A LOSSunExclusion ',num2str(params.sens.LOS_sun_excl)]);
    root.ExecuteCommand(['Graphics */Satellite/Satellite_',num2str(i),'/Sensor/sens_A Show Off']);
    %set ephemeris
    root.ExecuteCommand(['SetState */Satellite/Satellite_',num2str(i),' Classical J2Perturbation UseScenarioInterval ',num2str(params.timestep),' J2000 "',epoch,'" ' ,num2str(orbital_state(i,1)*1000),' ', num2str(orbital_state(i,2)),' ', num2str(orbital_state(i,3)),' ', num2str(orbital_state(i,5)),' ', num2str(mod(orbital_state(i,4)*180/pi,360)),' ', num2str(mod(orbital_state(i,6)*180/pi,360))]);
    
    %assign assets
    root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Asset */Satellite/Satellite_',num2str(i),'/Sensor/sens_A Assign']);
    root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Asset */Satellite/Satellite_',num2str(i),'/Sensor/sens_A Activate']);
    root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Interval UseScenarioInterval']);
end

% Example of adding a specific object to a constellation...
myConstObj = scenObj.Children.New('eConstellation', 'Flock_Constellation');
availObjs = myConstObj.Objects.AvailableObjects;
for i=1:length(availObjs)
myConstObj.Objects.Add(strcat('*/',availObjs{i}));
end
root.EndUpdate()
%% Compute Coverage & Get Data
% IAgCoverageDefinition coverage: Coverage object
advanced = cov_def.Advanced;
advanced.AutoRecompute = true;
advanced.DataRetention = 'eAll';
advanced.SaveMode = 'eSaveAccesses';

root.ExecuteCommand('Cov */CoverageDefinition/Space_Cov_Def Asset */Constellation/Flock_Constellation Assign');
root.ExecuteCommand('Cov */CoverageDefinition/Space_Cov_Def Asset */Constellation/Flock_Constellation Separate');

% average revisit time
fom = cov_def.Children.New('eFigureOfMerit', 'RevisitTime');
fom.SetDefinitionType('eFmRevisitTime');
fom.Definition.SetComputeType('eMaximum');
root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Access ParallelCompute']);

% Find Average FOM value
overallValDP = fom.DataProviders.GetDataPrvFixedFromPath('Overall Value');
Result_1 = overallValDP.Exec();
avgRevisitTime = cell2mat(Result_1.DataSets.GetDataSetByName('Maximum').GetValues)/3600

%% 

%advanced.AutoRecompute = false;
% tic
% planarRevisits=[];
% for abbr=['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H']
%     root.BeginUpdate()
%     for satPlane=1:40
%         root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Asset */Satellite/Satellite_',num2str(i),'/Sensor/sens_A Deactivate']);
%         root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Asset */Satellite/Satellite_',num2str(i),' Deactivate']);
%     end
%     root.EndUpdate()
%     root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Access ParallelCompute']);
%     overallValDP = fom.DataProviders.GetDataPrvFixedFromPath('Overall Value');
%     Result_1 = overallValDP.Exec();
%     avgRevisitTime = cell2mat(Result_1.DataSets.GetDataSetByName('Average').GetValues)/3600
%     planarRevisits=horzcat(planarRevisits,avgRevisitTime);
% end
% toc
%% 

tic
PermRevisits=[];
revisits=[];
random_removals=randperm(numPlane*satPlane);
for remove=0:((numPlane*satPlane/10)-1)
    root.BeginUpdate()
    for i=random_removals(((10*remove)+1):(10*(remove+1)))
        root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Asset */Satellite/Satellite_',num2str(i),'/Sensor/sens_A Deactivate']);
        root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Asset */Satellite/Satellite_',num2str(i),' Deactivate']);
    end
    root.EndUpdate()    
    root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Access ParallelCompute']);
    overallValDP = fom.DataProviders.GetDataPrvFixedFromPath('Overall Value');
    Result_1 = overallValDP.Exec();
    avgRevisitTime = cell2mat(Result_1.DataSets.GetDataSetByName('Maximum').GetValues)/3600;
    revisits=horzcat(revisits,avgRevisitTime);
end
PermRevisits=vertcat(PermRevisits,revisits);
toc
%     root.BeginUpdate()
%     for i=random_removals(1:end)
%         try
%         root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Asset */Satellite/Satellite_',num2str(i),'/Sensor/sens_A Activate']);
%         root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Asset */Satellite/Satellite_',num2str(i),' Activate']);
%         catch
%         end
%     end
%     root.EndUpdate()
% PermRevisits=vertcat(PermRevisits,revisits);
% toc

revisits_Flock=mean(PermRevisits);
save('revisits_Flock.mat');
% advanced.AutoRecompute = true;
% root.EndUpdate()
% root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Access ParallelCompute']);
% overallValDP = fom.DataProviders.GetDataPrvFixedFromPath('Overall Value');
% Result_1 = overallValDP.Exec();
% avgRevisitTime = cell2mat(Result_1.DataSets.GetDataSetByName('Average').GetValues)/3600
% toc