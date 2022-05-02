% clc;
% clearvars;
format shortG
warning('off','all')
warning
figure()
tiledlayout(2,4)
load('cum_radiation_all_new.mat');

tStart=tic;
time_horizon=11;     %length of simulation in years
NumMonteCarlos=50; %Number of Simulations to Run

%% Pre-Allocate a Matrix for each KPP

NPV_All=zeros(1,NumMonteCarlos); %Preallocate an array for the NPV of each Monte Carlo simulation so we can build the histogram and target curve at the end
Performance_All=zeros(1,NumMonteCarlos);
Radiation_Failures=zeros(1,NumMonteCarlos);
demand_array_all=zeros(NumMonteCarlos,(time_horizon*12)+1);
resilience_array_all=zeros(NumMonteCarlos,(time_horizon*12)+1);
% S5=zeros(1,NumMonteCarlos);
% S4=zeros(1,NumMonteCarlos);
% S3=zeros(1,NumMonteCarlos);
% S2=zeros(1,NumMonteCarlos);
% S1=zeros(1,NumMonteCarlos);
% S0=zeros(1,NumMonteCarlos);
cum_radiation_total=zeros(1,NumMonteCarlos);
% Total_Storm_Events=zeros(1,NumMonteCarlos);
% Unshielded_Radiation_Total=zeros(1,NumMonteCarlos);
Run_Time=[0];

%% Initialize STK Integration - Step Through This Using Run and Advance to Give STK Time to Setup
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
        root.NewScenario('WorldView_Scenario');
    end
    scenObj = root.CurrentScenario; %get the scenario root, its of type IAgScenario
end

%% Initialize Parameters

%Setup time interval
params.epoch_start = datenum('01-Apr-2014'); %date number (for STK)
params.sim_duration = 7; %days
params.timestep = 60; %sec

%Setup sensor parameters
params.sens.vert_FOV_half_angle = 2; %deg
params.sens.horiz_FOV_half_angle = 1.6; %deg
params.sens.target_lighting = 'DirectSun'%DirectSun, PenumbraDirectSun
params.sens.LOS_sun_excl = 60; %deg
params.sens.A_FOV_pointing_azel = [0,90,params.sens.vert_FOV_half_angle,params.sens.horiz_FOV_half_angle];

%Setup constellation parameters
params.const.num_planes = 4;
params.const.num_sats_per_plane = 1;
params.const.F = 1;
params.const.inc = 97.5; %deg
params.const.alt = 617; %km

params.cov.granularity = 3; %deg

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

%% Create Satellites

    %WorldView-1
    satellite = scenObj.Children.New('eSatellite', 'Satellite_1');    
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeAltitude';  % Changes from Ecc/Inc to Perigee/Apogee Altitude
    keplerian.LocationType = 'eLocationTrueAnomaly'; % Makes sure True Anomaly is being used
    keplerian.SizeShape.PerigeeAltitude = 617;      % km
    keplerian.SizeShape.ApogeeAltitude = 617;       % km
    keplerian.Orientation.Inclination = 97.3853;    % deg
    keplerian.Orientation.ArgOfPerigee = 311.942;   % deg
    keplerian.Orientation.AscNode.Value = 233.378;  % deg
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    sensor = satellite.Children.New('eSensor','Sensor_A');
    sensor.CommonTasks.SetPatternSimpleConic(41.6,0.3);
    root.ExecuteCommand(['SetConstraint */Satellite/Satellite_1','/Sensor/Sensor_A LOSSunExclusion ',num2str(params.sens.LOS_sun_excl)]);
    
    %WorldView-2
    satellite = scenObj.Children.New('eSatellite', 'Satellite_2');    
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeAltitude';  % Changes from Ecc/Inc to Perigee/Apogee Altitude
    keplerian.LocationType = 'eLocationTrueAnomaly'; % Makes sure True Anomaly is being used
    keplerian.SizeShape.PerigeeAltitude = 617;      % km
    keplerian.SizeShape.ApogeeAltitude = 617;       % km
    keplerian.Orientation.Inclination = 98.4879;    % deg
    keplerian.Orientation.ArgOfPerigee = 128.971;   % deg
    keplerian.Orientation.AscNode.Value = 189.127;  % deg
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    sensor = satellite.Children.New('eSensor','Sensor_A');
    sensor.CommonTasks.SetPatternSimpleConic(41.6,0.3);
    root.ExecuteCommand(['SetConstraint */Satellite/Satellite_2','/Sensor/Sensor_A LOSSunExclusion ',num2str(params.sens.LOS_sun_excl)]);
    
    %WorldView-3
    satellite = scenObj.Children.New('eSatellite', 'Satellite_3');    
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeAltitude';  % Changes from Ecc/Inc to Perigee/Apogee Altitude
    keplerian.LocationType = 'eLocationTrueAnomaly'; % Makes sure True Anomaly is being used
    keplerian.SizeShape.PerigeeAltitude = 617;      % km
    keplerian.SizeShape.ApogeeAltitude = 617;       % km
    keplerian.Orientation.Inclination = 97.8561;    % deg
    keplerian.Orientation.ArgOfPerigee = 259.068;   % deg
    keplerian.Orientation.AscNode.Value = 189.155;  % deg
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    sensor = satellite.Children.New('eSensor','Sensor_A');
    sensor.CommonTasks.SetPatternSimpleConic(41.6,0.3);
    root.ExecuteCommand(['SetConstraint */Satellite/Satellite_3','/Sensor/Sensor_A LOSSunExclusion ',num2str(params.sens.LOS_sun_excl)]);
    
    %WorldView-4
    satellite = scenObj.Children.New('eSatellite', 'Satellite_4');    
    keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
    keplerian.SizeShapeType = 'eSizeShapeAltitude';  % Changes from Ecc/Inc to Perigee/Apogee Altitude
    keplerian.LocationType = 'eLocationTrueAnomaly'; % Makes sure True Anomaly is being used
    keplerian.SizeShape.PerigeeAltitude = 617;      % km
    keplerian.SizeShape.ApogeeAltitude = 617;       % km
    keplerian.Orientation.Inclination = 97.98;    % deg
    keplerian.Orientation.ArgOfPerigee = 0;   % deg
    keplerian.Orientation.AscNode.Value = 148;  % deg
    satellite.Propagator.InitialState.Representation.Assign(keplerian);
    satellite.Propagator.Propagate;
    sensor = satellite.Children.New('eSensor','Sensor_A');
    sensor.CommonTasks.SetPatternSimpleConic(41.6,0.3);
    root.ExecuteCommand(['SetConstraint */Satellite/Satellite_4','/Sensor/Sensor_A LOSSunExclusion ',num2str(params.sens.LOS_sun_excl)]); 
%% Setup Coverage Definition

% Create Coverage Definitions
cov_def = scenObj.Children.New('eCoverageDefinition', 'Space_Cov_Def');

% set up calculation grids
cov_def.Grid.BoundsType = 'eBoundsLatLonRegion';
cov_def.Grid.Bounds.MinLongitude = -180;
cov_def.Grid.Bounds.MaxLongitude = 0;
cov_def.Grid.Bounds.MinLatitude = -60;
cov_def.Grid.Bounds.MaxLatitude = 60;
cov_def.Grid.Resolution.LatLon = 3;

%% Make Constellation

% Example of adding a specific object to a constellation...
myConstObj = scenObj.Children.New('eConstellation', 'WorldView_Constellation');
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

root.ExecuteCommand('Cov */CoverageDefinition/Space_Cov_Def Asset */Constellation/WorldView_Constellation Assign');
root.ExecuteCommand('Cov */CoverageDefinition/Space_Cov_Def Asset */Constellation/WorldView_Constellation Separate');

% average revisit time
fom = cov_def.Children.New('eFigureOfMerit', 'RevisitTime');
fom.SetDefinitionType('eFmRevisitTime');
fom.Definition.SetComputeType('eMaximum');
root.ExecuteCommand(['Cov */CoverageDefinition/Space_Cov_Def Access ParallelCompute']);

% Find Average FOM value
overallValDP = fom.DataProviders.GetDataPrvFixedFromPath('Overall Value');
Result_1 = overallValDP.Exec();
Base_RR = cell2mat(Result_1.DataSets.GetDataSetByName('Maximum').GetValues)/3600

%% Input Variables
%Includes  spacecraft design variables as well as parameters for solar
%weather, costs, launches, etc
SatName='WorldView';
numSeeds=100;
ComponentType=1;

%Spacecraft and Orbit Parameters
sc_height=.1;        %spacecraft height in meters
sc_length=.1;        %spacecraft length in meters
sc_width=.3;         %spacecraft width in meters
sc_mass=5;           %spacecraft mass in kg
rad_tolerance=30000; %component radiation tolerance in rad
sc_shielding=80;     %spacecraft shielding in mils (must be between 3.937 and 200)
image_cap=1590000;   %imaging capacity per satellite per day (km^2)

manufac_time=1;      %manufacturing time per satellite in months
period=93.5;         %orbital period in minutes

%Imaging Parameters
marketable_ratio=0.1;%percentage of images that are marketable (aka not over ocean or uninhabited land)
download_ratio=0.3;  %percentage of images that can actually be downloaded without being written over or lost
cloud_coverage=0.6;  %percentage of images not ruined by cloud coverage

%Constellation Architecture Parameters
op_orbit_pop=200;    %number of satellites in constellation at t=0
park_orbit_pop=0;    %number of satellites in parking at t=0
ground_spare_pop=0;  %number of satellite ground spares at t=0;

%Launch Parameters
mtb_launches=67;     %mean time between launches in days
lop_time_mean=75;    %mean launch order processing time in days
lop_time_sdev=15;    %standard deviation of launch order processing time in days

%Cost Parameters
orbit_op_cost=10000; %operating cost per satellites in orbit in $
holding_cost=5000;   %holding cost per ground spare in $
manufac_cost=1790000; %manufacturing cost per satellite in $
fuel_cost=1000;      %cost for fuel for manuevering in $/kg
fixed_launch_cost=10000000; %10000000;  %fixed overhead cost for launching in $
variable_launch_cost=5000;  %variable cost for launching in $/kg
discount_rate=.25/12;    %discount rate for simulation (monthly)
learning_rate=.8;    %learning rate for simulation (monthly)
base_image_val=50000;  %value in dollars per Mkm^2 imaged

start_time=0; %months into solar cycle to start (0 if at solar minimum)


%% Base Calculations and Parameters
base_revisit_rate=period/op_orbit_pop; %revisit rate at t=0

storm_length_mean=2.4506;              %distributional mean fit for storm length in hours
storm_length_sdev=1.4658;              %distributional standard deviation for storm length in hours

full_period=[start_time:1:start_time+(time_horizon*12)];     %#ok<NBRAK> %Creates an array of months from 0 to the end of the time horizon

background_radiation=4269*(1/3)^(sc_shielding/10); %Calculates background radiation added per time step based on sc shielding

AoC_sat=image_cap*marketable_ratio*download_ratio*cloud_coverage*(365/12);

%Bathtub Curve based on Weibull Estimation of Lifecycle
PNZ=.9;
beta=.4032;
theta=365;
AverageReliability=PNZ*exp(-([1:548]/theta).^beta);
AverageReliabilityFlip=fliplr(-PNZ*exp(-([1:547]/theta).^beta))+(2*AverageReliability(end));
BathtubCurve=horzcat(AverageReliability,AverageReliabilityFlip);

shielding_mass=((sc_shielding*0.00254)*2.7)*(sc_width*sc_length)/1000;
%% Starting Populations

%We calculate the radiation twice. The first time is so we can estimate the
%leftover readiation on the first generation of satellites due to the previous solar cycle when this
%simulation begins at t=0. The second is to simulate the radiation the satellites will experience during the current solar cycle.
num_storms=zeros(1,length(full_period));
storm_intensity=zeros(10,length(full_period));
storm_length=zeros(10,length(full_period));
unshielded_radiation=zeros(1,length(full_period));
shielded_radiation=zeros(1,length(full_period));
cum_radiation=zeros(1,length(full_period));
Average_Performance=zeros(1,length(NumMonteCarlos));

for MonteCarlo=1:NumMonteCarlos
    
    seed=randi([1,100000],1,1); %instead of generating storms using the distributions for each monte carlo, we make 100,000 storm scenarios and choose a random one using a seed and the cum_radiation_all file
    
    cum_radiation=cum_radiation_all(seed,1:length(full_period));%pick a random storm scenario
    
    %Apply shielding to storm radiation
    if sc_shielding>0
        cum_radiation=cum_radiation.*(3.8254*(sc_shielding^-1.131)); %based on log-log space fit
    end
    
    %Apply shielding to background radiation
    for i=1:width(cum_radiation)    
    cum_radiation(i)=cum_radiation(i)+(background_radiation*i);
    end
    
    holding_pop=0;       %population of satellites in holding at t=0
    parking_pop=0;       %population of satellites in parking orbit at t=0
    Gen=[1:1:4]';
    NumSats=[42,42,43,43]';
    sats=table(Gen,NumSats);
    sats.Radiation=cum_radiation(length(full_period))-[cum_radiation(length(full_period)-27),cum_radiation(length(full_period)-18),cum_radiation(length(full_period)-9),cum_radiation(length(full_period))]';
    sats.Status=repmat(["Active"]',4,1);
    sats.Location=repmat(["Operational Orbit"]',4,1);
    sats.Age=[27,18,9,0]';
    operational_pop=sum(NumSats);
    %% Simulation Pt. 1 Generating Storms, Launch Availability, and Demand for the Solar Cycle
    
    %Preallocate Arrays for Speed
    sats_manufactured=zeros(1,length(full_period));
    AreaOfCoverage=zeros(1,length(full_period));
    Revisit=zeros(1,length(full_period));
    Performance=zeros(1,length(full_period));
    manufac_queue=cell2table(cell(0,2), 'VariableNames', {'NumBuilt','FinishMonth'});
    sats_launched=zeros(1,length(full_period));
    launch_queue=cell2table(cell(0,4), 'VariableNames', {'OpNumLaunched','OpFinishMonth','ParkNumLaunched','ParkFinishMonth'});
    op_orbit_pop_array=horzcat(op_orbit_pop,zeros(1,length(full_period)-1));
    park_orbit_pop_array=horzcat(park_orbit_pop,zeros(1,length(full_period)-1));
    ground_spare_pop_array=horzcat(ground_spare_pop,zeros(1,length(full_period)-1));
    OpOrbitPop=zeros(1,length(full_period));
    ParkOrbitPop=zeros(1,length(full_period));
    Holding_Costs_Array=zeros(1,length(full_period));
    Manufacturing_Costs_Array=zeros(1,length(full_period));%horzcat(max((ground_spare_pop*manufac_cost),1),zeros(1,length(full_period)-1));
    Launch_Costs_Array=zeros(1,length(full_period));
    Total_Costs_Array=zeros(1,length(full_period));
    Revenue_Array=zeros(1,length(full_period));
    Cashflow_Array=zeros(1,length(full_period));
    
    %Generate storms for solar cycle using seeds again
    seed=randi([1,100000],1,1); %instead of generating storms using the distributions for each monte carlo, we make 100,000 storm scenarios and choose a random one using a seed and the cum_radiation_all file
    
    cum_radiation=cum_radiation_all(seed,1:length(full_period));%pick a random storm scenario
    
    %Apply shielding to storm radiation
    if sc_shielding>0
        cum_radiation=cum_radiation.*(3.8254*(sc_shielding^-1.131)); %based on log-log space fit
    end
    
    %Apply shielding to background radiation
    for i=1:width(cum_radiation)    
    cum_radiation(i)=cum_radiation(i)+(background_radiation*i);
    end 

    %Launch Availability
    mtb_launch_array=-log(1-(.9*rand(length(full_period),1)))*mtb_launches; %Time needed for launch window
    lop_time_array=norminv(rand(length(full_period),1),lop_time_mean,lop_time_sdev); %Time needed to process and integrate
    tot_launch_time=ceil((mtb_launch_array+lop_time_array)/30); %Total number of months needed
    
    %Binomial Tree Brownian Motion Demand
    demand_array=zeros(1,length(full_period));
    demand_array(1)=op_orbit_pop*AoC_sat/1000000;
    volatility=.1/12;
    drift_param=.01/12;
    time_step=1;
    
    UpOrDown=rand();
    if UpOrDown>=.5
        UpOrDown=1;
    else
        UpOrDown=-1;
    end
    
    for i=2:length(full_period)
        if mod(i,12)==0 %this logic creates the binomial part of the binomial brownian motion
            UpOrDown=rand(); %every year, there is a 50/50 chance the demand for the next year will be
            if UpOrDown>=.5 %drifting upward or downward.
                UpOrDown=1;
            else
                UpOrDown=-1;
            end
        end
        drift=demand_array(i-1)*drift_param*time_step;
        uncertainty=norminv(rand(),0,1)*sqrt(time_step)*volatility*demand_array(i-1);
        change=(drift+uncertainty)*UpOrDown;
        demand_array(i)=demand_array(i-1)+change;
    end
    %% Simulation Pt. 2 Moving Through Time Steps
    
    for t=2:length(full_period)
        
        %% Failures
        Active=find(sats.Status=="Active");
        %Radiation
        %Update Radiation for Current Satellites
        sats.Radiation(Active)=sats.Radiation(Active)+cum_radiation(t)-cum_radiation(t-1);
        
        %Check Number of Radiation Failures
        sats.Status(find(sats.Radiation>=rad_tolerance))="Radiation Failure";
        sats.NumSats(find(sats.Status=="Radiation Failure"))=0;
        
        %Check SEU Failures
        SEU_Failures=poissinv(rand(),SEU_Rate); %SEUs
        Crit_SEU_Failures=sum([rand(1,SEU_Failures)]<=Crit_SEU_Rate); %Roll a random number for each SEU, for each one below the critical SEU rate, there is a permanent sat lost
        
        %Lifecycle Failures
        sats.NumSats(Active)=floor(sats.NumSats(Active)'.*BathtubCurve(sats.Age(Active)+1)); %Weibull Based Bathtub Curve
        sats.Status(find(((sats.Status~="Radiation Failure")==(sats.NumSats==0))))="Lifecycle Failure";
        
        %% Management Decisions
        %How many to manufacture and when, where to launch and when, when to
        %manuever from parking orbit to operational orbit
        %Both are based on decision rules (e.g. Launch when capacity<demand)
        
        %Manufacturing
        if ground_spare_pop<=80
            sats_manu_start(t)=5;
        else
            sats_manu_start(t)=2;
        end %Decision Rule for Start Manufacting
        
        manufac_update={sats_manu_start(t),t+manufac_time};
        
        manufac_queue=[manufac_queue;manufac_update]; %Update Queue with how many are being manufactured and when they finish based on manufacturing time
        
        if isempty(manufac_queue.NumBuilt(find(t==manufac_queue.FinishMonth)))==0
            ground_spare_pop_array(t)=ground_spare_pop_array(t-1)+sum(manufac_queue.NumBuilt(find(t==manufac_queue.FinishMonth))); %Update ground spare population when manufacturing is complete
        end
        
        %Launching
        %Decision Rule for When and Where to Launch
        if ground_spare_pop_array(t)>=40 %Need to calculate launch costs separately if there is a launch to operational orbit AND parking orbit
            op_launch_start(t)=40;
            parking_launch_start(t)=0;
        else
            op_launch_start(t)=0;
            parking_launch_start(t)=0;
        end
        
        %Updating Launch Queue
        launch_update={op_launch_start(t),t+tot_launch_time(t),parking_launch_start(t),tot_launch_time(t)};
        launch_queue=[launch_queue;launch_update];
        
        
        %Updating Launches that Have Arrived
        OpArrived=launch_queue.OpNumLaunched(find(t==launch_queue.OpFinishMonth));
        if isempty(OpArrived)==0
            op_orbit_pop_array(t)=op_orbit_pop_array(t-1)+sum(OpArrived);
            if sum(OpArrived)>0
                sats=[sats;{max(sats.Gen)+1,sum(OpArrived),0,"Active","Operational Orbit",0}];
            end
        end
        
        ParkArrived=launch_queue.ParkNumLaunched(find(t==launch_queue.ParkFinishMonth)); %Single Event Upsets don't matter because not operational but SAA radiation accumulates
        if isempty(ParkArrived)==0
            park_orbit_pop_array(t)=park_orbit_pop_array(t-1)+sum(ParkArrived);
            if sum(ParkArrived)>0
                sats=[sats;{max(sats.Gen)+1,sum(ParkArrived),0,"Active","Parking Orbit",0}];
            end
        end
        
        
        %Performance Calculations (Include Temporary SEUs Here)
        operationalPop=find((sats.Status=="Active")==(sats.Location=="Operational Orbit")); %Recheck Active Operational Orbit Sats
        OpOrbitPop(t)=sum(sats.NumSats(operationalPop));%Update Operational Orbit Population
        parkingPop=find((sats.Status=="Active")==(sats.Location=="Parking Orbit")); %Recheck Active Parking Orbit Sats
        ParkOrbitPop(t)=sum(sats.NumSats(parkingPop));%Update Parking Orbit Population
        AreaOfCoverage(t)=OpOrbitPop(t).*AoC_sat/1000000;
        Revisit(t)=OpOrbitPop(t).*Base_RR; %Need to update with variable
        Performance(t)=(.5*Revisit(t)/Revisit(2))+(.5*AreaOfCoverage(t)/AreaOfCoverage(2));
        %Financial Calculations
        Holding_Costs_Array(t)=((OpOrbitPop(t)+ParkOrbitPop(t))*orbit_op_cost)+(ground_spare_pop_array(t)*holding_cost);
        Manufacturing_Costs_Array(t)=sats_manu_start(t)*manufac_cost*((sum(sats_manu_start(1:t-1)+1))^(1-(log2(1/learning_rate))));
        Month_Launch_Cost=0;
        Launch_Costs_Array(t)=((launch_queue.OpNumLaunched(t-1)>0+launch_queue.ParkNumLaunched(t-1)>0)*fixed_launch_cost)+...
            (launch_queue.ParkNumLaunched(t-1)+launch_queue.OpNumLaunched(t-1))*variable_launch_cost*(sc_mass+shielding_mass);
        Total_Costs_Array(t)=(Holding_Costs_Array(t)+Manufacturing_Costs_Array(t)+Launch_Costs_Array(t))*-1;%add manuevering cost
       
        if Revisit(t)==0
            Revenue_Array(t)=0;
        else
            Revenue_Array(t)=AreaOfCoverage(t)*base_image_val*log(Revisit(t));
        end
        Cashflow_Array(t)=Revenue_Array(t)+Total_Costs_Array(t);
        %DCF_Array(t)=Cashflow_Array(t)/((1+discount_rate)^t)
        %Update Satellite Ages at End of Period
        sats.Age=sats.Age+1;
    end
    NPV_All(MonteCarlo)=pvvar(Cashflow_Array,discount_rate)+(ground_spare_pop*holding_cost); %Need to add back in costs at t=0. This assumes no launches, manuevering, or manufacturing happening at t=0 but that there is a starting population with associated holding costs.
    Performance_All(MonteCarlo)=mean(Performance);
    Radiation_Failures(MonteCarlo)=sum(sats.Status=="Radiation Failure");
    cum_radiation_total(MonteCarlo)=cum_radiation(end);
    demand_array_all(MonteCarlo,:)=demand_array;
%     S5(MonteCarlo)=length(find(storm_intensity>=100000));
%     S4(MonteCarlo)=length(find(storm_intensity>=10000))-S5(MonteCarlo);
%     S3(MonteCarlo)=length(find(storm_intensity>=1000))-S5(MonteCarlo)-S4(MonteCarlo);
%     S2(MonteCarlo)=length(find(storm_intensity>=100))-S5(MonteCarlo)-S4(MonteCarlo)-S3(MonteCarlo);
%     S1(MonteCarlo)=length(find(storm_intensity>=10))-S5(MonteCarlo)-S4(MonteCarlo)-S3(MonteCarlo)-S2(MonteCarlo);
%     S0(MonteCarlo)=length(find(storm_intensity>=1))-S5(MonteCarlo)-S4(MonteCarlo)-S3(MonteCarlo)-S2(MonteCarlo)-S1(MonteCarlo);
%     Total_Storm_Events(MonteCarlo)=length(find(storm_intensity>0));
%     Unshielded_Radiation_Total(MonteCarlo)=unshielded_radiation(end);
    tEnd=toc(tStart);
    if mod(MonteCarlo,100)==0
        Run_Time=horzcat(Run_Time,tEnd-Run_Time(end));
        fprintf("\n");
        toc(tStart)
        fprintf("Estimated Time Until Completion: %.2f Seconds (%.2f Minutes). \nTotal Monte Carlos Runs = %i / %i (%.2f Percent Complete).", ((NumMonteCarlos-MonteCarlo)/100*mean(Run_Time(2:end))),((NumMonteCarlos-MonteCarlo)/100*mean(Run_Time(2:end)))/60 ,MonteCarlo,NumMonteCarlos,MonteCarlo/NumMonteCarlos*100);
    end
end

    ENPV=mean(NPV_All)
    
    ax1=nexttile(1);
    histogram(NPV_All,100)
    title('Histogram of NPVs')
    
    ax2=nexttile(2);
    hold(ax2,'on')
    cdfplot(NPV_All)
    hold(ax2,'off')
    title('Target Curve for NPV')
    %legend(string([10:10:100]))
    
    ax3=nexttile(3);
    histogram(Performance_All)
    title('Histogram of Average Performance')
    
    ax4=nexttile(4);
    hold(ax4,'on')
    cdfplot(Performance_All)
    hold(ax4,'off')
    title('Target Curve for Performance')
    %legend(string([10:10:100]))
        
    ax5=nexttile(5);
    histogram(Radiation_Failures)
    title('Histogram of Radiation Failures')
    
    ax6=nexttile(6);
    hold(ax6,'on')
    cdfplot(Radiation_Failures)
    hold(ax2,'off')
    title('Target Curve of Radiation Failures')
    %legend(string([10:10:100]))
        
    ax7=nexttile(7);
    histogram(cum_radiation_total/1000)
    title('Histogram of Total Unshielded Radiation for Solar Cycle')
    
    ax8=nexttile(8);
    hold(ax8,'on')
    cdfplot(cum_radiation_total/1000)
    hold(ax8,'off')
    title('Target Curve of Total Unshielded Radiation for Solar Cycle')
    %legend(string([10:10:100]))
        
%     subplot(2,4,1)
%     histogram(NPV_All,100)
%     title('Histogram of NPVs')
%     
%     subplot(2,4,2)
%     cdfplot(NPV_All)
%     title('Target Curve for NPV')
%     
%     subplot(2,4,3)
%     histogram(Performance_All)
%     title('Histogram of Average Performance')
%     
%     subplot(2,4,4)
%     cdfplot(Performance_All)
%     title('Target Curve for Performance')
%     
%     subplot(2,4,5)
%     histogram(Radiation_Failures)
%     title('Histogram of Radiation Failures')
%     
%     subplot(2,4,6)
%     cdfplot(Radiation_Failures)
%     title('Target Curve of Radiation Failures')
%     
%     subplot(2,4,7)
%     histogram(cum_radiation_total)
%     title('Histogram of Total Unshielded Radiation for Solar Cycle')
%     
%     subplot(2,4,8)
%     cdfplot(cum_radiation_total)
%     title('Target Curve of Total Unshielded Radiation for Solar Cycle')
    
%     a=length(find(storm_intensity>=100000));
%     b=length(find(storm_intensity>=10000))-a;
%     c=length(find(storm_intensity>=1000))-a-b;
%     d=length(find(storm_intensity>=100))-a-b-c;
%     e=length(find(storm_intensity>=10))-a-b-c-d;
%     f=length(find(storm_intensity>=1))-a-b-c-d-e;
%     fprintf("Storms greater than 100,000: %i" ,mean(S5));
%     fprintf("\nStorms between 10,000 and 100,000: %i" ,mean(S4));
%     fprintf("\nStorms greater than 1,000: %i" ,mean(S3));
%     fprintf("\nStorms greater than 100: %i" ,mean(S2));
%     fprintf("\nStorms greater than 10: %i" ,mean(S1));
%     fprintf("\nStorms greater than 1: %i" ,mean(S0));
%     fprintf("\nNumber of Storms: %i" ,mean(Total_Storm_Events));
    fprintf("\nMean Number of Radiation Failures Over Solar Cycle (%i Simulations): %i\n" ,NumMonteCarlos, mean(Radiation_Failures));
    
    figure()
    hold on
    for i=1:height(demand_array_all)
        plot(demand_array_all(i,:))
    end
toc(tStart)