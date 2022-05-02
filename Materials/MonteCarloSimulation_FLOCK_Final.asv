% clc;
clearvars;
format shortG
warning('off','all')
warning
load('D:\Materials\revisits_Flock.mat');
revisits_array=fliplr(revisits_Flock);

tStart=tic;
time_horizon=11;     %length of simulation in years
NumMonteCarlos=1000; %Number of Simulations to Run

%% Pre-Allocate a Matrix for each KPP

NPV_All=zeros(1,NumMonteCarlos); %Preallocate an array for the NPV of each Monte Carlo simulation so we can build the histogram and target curve at the end
Performance_All=zeros(1,NumMonteCarlos);
Radiation_Failures=zeros(1,NumMonteCarlos);
%demand_array_all=zeros(NumMonteCarlos,(time_horizon*12)+1);
resilience_array_all=zeros(NumMonteCarlos,(time_horizon*12)+1);
cum_radiation_total=zeros(1,NumMonteCarlos);
total_cycle_SEUs_monte=zeros(1,NumMonteCarlos);
Lifecycle_Failures_monte=zeros(1,NumMonteCarlos);
Total_Costs_monte=zeros(1,NumMonteCarlos);
Sat_Age_monte=zeros(1,NumMonteCarlos);
Robustness_S4_monte=[];
Robustness_S5_monte=[];
Robustness_Max_monte=[];
Resilience_ShortTerm_S4=[];
Resilience_ShortTerm_S5=[];
No_Recoveries_S5=0;
No_Recoveries_S4=0;
Run_Time=[0];
Agility_S5=[];
Agility_S4=[];
No_Recoveries_S4_monte=[];
No_Recoveries_S5_monte=[];

%% Input Variables
%Includes  spacecraft design variables as well as parameters for solar
%weather, costs, launches, etc
SatName='Flock';
numSeeds=10000;
Component_Type=3;%1 = Space-Grade, 2 = Space-Tolerant, 3 = COTS - only used for SEUs

%Spacecraft and Orbit Parameters
sc_length=.1;        %spacecraft length in meters
sc_width=.1;         %spacecraft width in meters
sc_mass=5;         %spacecraft mass in kg
rad_tolerance=25;     %component radiation tolerance in krad
sc_shielding=80;       %spacecraft shielding in mils (must be between 5 and 400)
image_cap=1590000;    %imaging capacity per satellite per day (km^2)
manufac_time=1;      %manufacturing time per satellite in months

%Imaging Parameters
marketable_ratio=0.1; %percentage of images that are marketable (aka not over ocean or uninhabited land)
download_ratio=0.3;   %percentage of images that can actually be downloaded without being written over or lost
cloud_coverage=0.6;   %percentage of images not ruined by cloud coverage

%Constellation Architecture Parameters
op_orbit_pop=200;       %number of satellites in constellation at t=0
park_orbit_pop=0;     %number of satellites in parking at t=0
ground_spare_pop=0;   %number of satellite ground spares at t=0;

%Launch Parameters
mtb_launches=67;      %mean time between launches in days
lop_time_mean=75;     %mean launch order processing time in days
lop_time_sdev=15;     %standard deviation of launch order processing time in days

%Cost Parameters
orbit_op_cost=1000;  %operating cost per satellites in orbit in $
holding_cost=10000;  %holding cost per ground spare in $
manufac_cost=1700000;%manufacturing cost per satellite in $
fuel_cost=1000;        %cost for fuel for manuevering in $/kg
fixed_launch_cost=10000000;%fixed overhead cost for launching in $
variable_launch_cost=5000;%variable cost for launching in $/kg
discount_rate=.20/12;  %discount rate for simulation (monthly)
learning_rate=.80;     %learning rate for simulation (monthly)
base_image_val=50000;  %value in dollars per Mkm^2 imaged

%Weibull Parameters
PNZ=.9;
K=0.9057;
beta1=0.4154;
theta1=558612.12;
beta2=5.06;
theta2=122.52;

start_time=0; %months into solar cycle to start (0 if at solar minimum)
%% Generate Storm Scenario Seeds
full_period=[start_time:1:start_time+(time_horizon*12)];     %#ok<NBRAK> %Creates an array of months from 0 to the end of the time horizon

[SEU_Probs,Rad_Accum,S4_Flags,S5_Flags]=StormSeedGenerator_Final(sc_shielding,numSeeds,Component_Type,full_period,SatName); %Comment if none of the inputs change to prevent re-calculating all seeds

%% Base Calculations and Parameters

AoC_sat_initial=image_cap*marketable_ratio*download_ratio*cloud_coverage*(365/12);
AoC_sat=image_cap*marketable_ratio*download_ratio*cloud_coverage*(365/12);

%Bathtub Curve based on Weibull Estimation of Lifecycle
Reliability=zeros(1,length(full_period));
Reliability(1)=((K*exp(-(1/theta1)^beta1))+((1-K)*exp(-(1/theta2)^beta2)));
for i=2:length(full_period)
    Reliability(i)=(((K*exp(-((i)/theta1)^beta1))+((1-K)*exp(-((i)/theta2)^beta2))))/Reliability(i-1);
end
shielding_mass=((sc_shielding*0.00254)*2.7)*(sc_width*sc_length)/1000;
%% Starting Populations

%We calculate the radiation twice. The first time is so we can estimate the
%leftover readiation on the first generation of satellites due to the previous solar cycle when this
%simulation begins at t=0. The second is to simulate the radiation the satellites will experience during the current solar cycle.
Performance=zeros(1,length(NumMonteCarlos));
total_cycle_radiation_monte=zeros(1,NumMonteCarlos);
for MonteCarlo=1:NumMonteCarlos
    
    seed=randi([1,numSeeds],1,1); %instead of generating storms using the distributions for each monte carlo, we make 100,000 storm scenarios and choose a random one using a seed and the cum_radiation_all file
    
    cum_radiation=Rad_Accum(seed,1:length(full_period));%pick a random storm scenario
    
    holding_pop=0;       %population of satellites in holding at t=0
    parking_pop=0;       %population of satellites in parking orbit at t=0
    Gen=[1:1:5]';
    NumSats=[40,40,40,40,40]';
    sats=table(Gen,NumSats);
    sats.Radiation=[sum(cum_radiation(length(full_period)-36:end)),sum(cum_radiation(length(full_period)-27:end)),sum(cum_radiation(length(full_period)-18:end)),sum(cum_radiation(length(full_period)-9:end)),cum_radiation(length(full_period))]';
    sats.Status=repmat(["Active"]',length(NumSats),1);
    sats.Location=repmat(["Operational Orbit"]',length(NumSats),1);
    sats.Age=[36,27,18,9,0]';%in months
    operational_pop=sum(NumSats);
    %% Simulation Pt. 1 Generating Storms, Launch Availability, and Demand for the Solar Cycle
    
    %Preallocate Arrays for Speed
    sats_manufactured=zeros(1,length(full_period));
    AreaOfCoverage=zeros(1,length(full_period));
    Revisit=zeros(1,length(full_period));
    Revisit(1)=revisits_array(21);
    Performance=zeros(1,length(full_period));
    manufac_queue=cell2table(cell(0,2), 'VariableNames', {'NumBuilt','FinishMonth'});
    sats_launched=zeros(1,length(full_period));
    launch_queue=cell2table(cell(0,4), 'VariableNames', {'OpNumLaunched','OpFinishMonth','ParkNumLaunched','ParkFinishMonth'});
    OpOrbitPop=horzcat(op_orbit_pop,zeros(1,length(full_period)-1));
    park_orbit_pop_array=horzcat(park_orbit_pop,zeros(1,length(full_period)-1));
    ground_spare_pop_array=horzcat(ground_spare_pop,zeros(1,length(full_period)-1));
    ParkOrbitPop=zeros(1,length(full_period));
    Holding_Costs_Array=zeros(1,length(full_period));
    Manufacturing_Costs_Array=zeros(1,length(full_period));%horzcat(max((ground_spare_pop*manufac_cost),1),zeros(1,length(full_period)-1));
    Launch_Costs_Array=zeros(1,length(full_period));
    Total_Costs_Array=zeros(1,length(full_period));
    Revenue_Array=zeros(1,length(full_period));
    Cashflow_Array=zeros(1,length(full_period));
    Lifecycle_Failures=zeros(1,length(full_period));
    Tot_SEU_Failures=zeros(1,length(full_period));
    sats_manu_start=zeros(1,length(full_period));
    Tot_Failures=zeros(1,length(full_period));
    total_cycle_radiation=zeros(1,length(full_period));
    Robustness=zeros(1,length(full_period));
    Num_Radiation_Failures=zeros(1,length(full_period));
    %Generate storms for solar cycle using seeds again
    seed=randi([1,numSeeds],1,1); %instead of generating storms using the distributions for each monte carlo, we make 100,000 storm scenarios and choose a random one using a seed and the cum_radiation_all file
    
    cum_radiation=Rad_Accum(seed,1:length(full_period));%pick a random storm scenario
    
    %Launch Availability
    mtb_launch_array=-log(1-(.9*rand(length(full_period),1)))*mtb_launches; %Time needed for launch window
    lop_time_array=norminv(rand(length(full_period),1),lop_time_mean,lop_time_sdev); %Time needed to process and integrate
    tot_launch_time=ceil((mtb_launch_array+lop_time_array)/30); %Total number of months needed
    
    %% Simulating Demand - Leaving in Because May be Useful in Some Applications
    %Is currently set up as random brownian motion for amount of area
    %demanded in Mkm^2 over time
    
    %     %Binomial Tree Brownian Motion Demand
    %     demand_array=zeros(1,length(full_period));
    %     demand_array(1)=op_orbit_pop*AoC_sat/1000000;
    %     volatility=.1/12;
    %     drift_param=.01/12;
    %     time_step=1;
    %
    %     UpOrDown=rand();
    %     if UpOrDown>=.5
    %         UpOrDown=1;
    %     else
    %         UpOrDown=-1;
    %     end
    %
    %     for i=2:length(full_period)
    %         if mod(i,12)==0 %this logic creates the binomial part of the binomial brownian motion
    %             UpOrDown=rand(); %every year, there is a 50/50 chance the demand for the next year will be
    %             if UpOrDown>=.5 %drifting upward or downward.
    %                 UpOrDown=1;
    %             else
    %                 UpOrDown=-1;
    %             end
    %         end
    %         drift=demand_array(i-1)*drift_param*time_step;
    %         uncertainty=norminv(rand(),0,1)*sqrt(time_step)*volatility*demand_array(i-1);
    %         change=(drift+uncertainty)*UpOrDown;
    %         demand_array(i)=demand_array(i-1)+change;
    %     end
    %% Simulation Pt. 2 Moving Through Time Steps
    SEU_Fails=zeros(sum(sats.NumSats),1);
    num_SEU_Fails=0;
    for t=2:length(full_period)
        
        %% Reactivate Temporary SEU Failures
        sats.Status(sats.Status=="SEU Failure")="Active";
        Active=find(sats.Status=="Active");
        
        %% Failures
        Initial_Active=sum(sats.NumSats(Active));
        %Radiation
        %Update Radiation for Current Satellites
        sats.Radiation(Active)=sats.Radiation(Active)+cum_radiation(t);
        total_cycle_radiation(t)=total_cycle_radiation(t-1)+cum_radiation(t);
        
        %Check Number of Radiation Failures
        sats.Status(find(sats.Radiation>=rad_tolerance))="Radiation Failure";
        sats.NumSats(find(sats.Status=="Radiation Failure"))=0;
        Active=find(sats.Status=="Active");
        Num_Radiation_Failures(t)=Initial_Active-sum(sats.NumSats(Active));
        
        %Lifecycle Failures
        for p=1:length(sats.NumSats(Active))
            Lifecycle_Fails=rand(sum(sats.NumSats(Active(p))),1)>Reliability(t);
            sats.NumSats(Active(p))=sats.NumSats(Active(p))-sum(Lifecycle_Fails);
            Lifecycle_Failures(t)=Lifecycle_Failures(t)+sum(Lifecycle_Fails);
        end
        %sats.NumSats(Active)=~Lifecycle_Fails;
        %sats.Status(Lifecycle_Fails)="Lifecycle Failure";
        %Lifecycle_Failures(t)=sum(Lifecycle_Fails);
        %         Lifecycle_Fails=0; %Uncomment if you want to disable lifecycle failures
        %         Lifecycle_Failures(t)=0;
        sats.Status(find(((sats.Status~="Radiation Failure")==(sats.NumSats==0))))="Lifecycle Failure";
        Active=find(sats.Status=="Active");
        
        num_SEU_Fails=0;
        %Check SEU Failures
        if isempty(SEU_Probs{seed,t})==0
            SEU_Fails=(sum((rand(sum(sats.NumSats(Active)),length(SEU_Probs(seed,t)))>SEU_Probs{seed,t}),2))>1;%Counts the number of sats that experience a failure per SPE SEU Probability
            %sats.Status(logical(SEU_Fails))="SEU Failure";
            num_SEU_Fails=sum(SEU_Fails);
        end
        %sats.Status(logical(SEU_Fails))="SEU Failure";
        Tot_SEU_Failures(t)=num_SEU_Fails;
        Active=find(sats.Status=="Active");
        End_Active=sum(sats.NumSats(Active));
        
        Tot_Failures(t)=Initial_Active-End_Active;
        %% Management Decisions
        %How many to manufacture and when, where to launch and when, when to
        %manuever from parking orbit to operational orbit
        %Both are based on decision rules (e.g. Launch when capacity<demand)
        OpOrbitPop(t)=max(sum(sats.NumSats(sats.Status=="Active")),0);
        %Manufacturing
        if (ground_spare_pop_array(t)<80 && (OpOrbitPop(t)+Tot_SEU_Failures(t))<(1.1*op_orbit_pop) && ((sum(sats.NumSats))+sum(sats_manu_start))<(1.1*op_orbit_pop)+sum(Lifecycle_Failures)+sum(Num_Radiation_Failures))
            sats_manu_start(t)=40;
        elseif ground_spare_pop_array(t)<100
            sats_manu_start(t)=5;
        else
            sats_manu_start(t)=0;
        end %Decision Rule for Start Manufacting
        
        manufac_update={sats_manu_start(t),t+manufac_time};
        
        manufac_queue=[manufac_queue;manufac_update]; %Update Queue with how many are being manufactured and when they finish based on manufacturing time
        
        if isempty(manufac_queue.NumBuilt(find(t==manufac_queue.FinishMonth)))==0
            ground_spare_pop_array(t)=ground_spare_pop_array(t-1)+sum(manufac_queue.NumBuilt(find(t==manufac_queue.FinishMonth))); %Update ground spare population when manufacturing is complete
        end
        
        %Launching
        %Decision Rule for When and Where to Launch
        lop_time_mean=75;
        if ground_spare_pop_array(t)>=40  && (OpOrbitPop(t)+Tot_SEU_Failures(t))<(1.2*op_orbit_pop) && (sum(launch_queue.OpNumLaunched)+sum(sats.NumSats))<(1.2*op_orbit_pop)+sum(Lifecycle_Failures)+sum(Num_Radiation_Failures)%Need to calculate launch costs separately if there is a launch to operational orbit AND parking orbit
            if (OpOrbitPop(t)+Tot_SEU_Failures(t))<(.8*op_orbit_pop)
                lop_time_mean=45;
            end
            op_launch_start(t)=40;
            parking_launch_start(t)=0;
        else
            op_launch_start(t)=0;
            parking_launch_start(t)=0;
        end
        
        %Updating Launch Queue
        launch_update={op_launch_start(t),t+tot_launch_time(t),parking_launch_start(t),tot_launch_time(t)};
        ground_spare_pop_array(t)=ground_spare_pop_array(t)-op_launch_start(t);
        launch_queue=[launch_queue;launch_update];
        
        
        %Updating Launches that Have Arrived
        OpArrived=launch_queue.OpNumLaunched(find(t==launch_queue.OpFinishMonth));
        OpOrbitPop(t)=max(sum(sats.NumSats(sats.Status=="Active"))+PNZ*sum(OpArrived)-Tot_SEU_Failures(t),0);
        if sum(OpArrived)>0
            sats=[sats;{max(sats.Gen)+1,PNZ*sum(OpArrived),0,"Active","Operational Orbit",0}];
        end
        
        ParkArrived=launch_queue.ParkNumLaunched(find(t==launch_queue.ParkFinishMonth)); %Single Event Upsets don't matter because not operational but SAA radiation accumulates
        if isempty(ParkArrived)==0
            park_orbit_pop_array(t)=park_orbit_pop_array(t-1)+sum(ParkArrived);
            if sum(ParkArrived)>0
                sats=[sats;{max(sats.Gen)+1,sum(ParkArrived),0,"Active","Parking Orbit",0}];
            end
        end
        
        
        %Performance Calculations (Include Temporary SEUs Here)
        operationalPop=find((sats.Status=="Active")==(sats.Location=="Operational Orbit"))-Tot_SEU_Failures(t); %Recheck Active Operational Orbit Sats
        parkingPop=find((sats.Status=="Active")==(sats.Location=="Parking Orbit")); %Recheck Active Parking Orbit Sats
        ParkOrbitPop(t)=sum(sats.NumSats(parkingPop));%Update Parking Orbit Population
        AreaOfCoverage(t)=OpOrbitPop(t).*AoC_sat/1000000;
        Revisit(t)=interp1([0:10:320],revisits_array,(OpOrbitPop(t)));
        Performance(t)=(.5*Revisit(1)/Revisit(t))+(.5*AreaOfCoverage(t)/(AoC_sat_initial*op_orbit_pop/1000000));
        Robustness(t)=Performance(t-1)-Performance(t);
        %Financial Calculations
        Holding_Costs_Array(t)=((OpOrbitPop(t)+ParkOrbitPop(t))*orbit_op_cost)+(ground_spare_pop_array(t)*holding_cost);
        Manufacturing_Costs_Array(t)=sats_manu_start(t)*manufac_cost*((sum(sats_manu_start(1:t-1))+1)^(-(log2(1/learning_rate))));
        Month_Launch_Cost=0;
        Launch_Costs_Array(t)=((launch_queue.OpNumLaunched(t-1)>0+launch_queue.ParkNumLaunched(t-1)>0)*fixed_launch_cost)+...
            (launch_queue.ParkNumLaunched(t-1)+launch_queue.OpNumLaunched(t-1))*variable_launch_cost*(sc_mass+shielding_mass);
        Total_Costs_Array(t)=(Holding_Costs_Array(t)+Manufacturing_Costs_Array(t)+Launch_Costs_Array(t))*-1;%add manuevering cost
        
        if Revisit(t)==168
            Revenue_Array(t)=0;
        else
            Revenue_Array(t)=AreaOfCoverage(t)*base_image_val*log(Revisit(t));
        end
        Cashflow_Array(t)=Revenue_Array(t)+Total_Costs_Array(t);
        %DCF_Array(t)=Cashflow_Array(t)/((1+discount_rate)^t)
        %Update Satellite Ages at End of Period
        sats.Age(Active)=sats.Age(Active)+1;
    end
    NPV_All(MonteCarlo)=pvvar(Cashflow_Array,discount_rate)+(ground_spare_pop*holding_cost); %Need to add back in costs at t=0. This assumes no launches, manuevering, or manufacturing happening at t=0 but that there is a starting population with associated holding costs.
    Performance_All(MonteCarlo)=mean(Performance);
    Radiation_Failures(MonteCarlo)=sum(sats.Status=="Radiation Failure");
    total_cycle_radiation_monte(MonteCarlo)=total_cycle_radiation(end);
    total_cycle_SEUs_monte(MonteCarlo)=sum(Tot_SEU_Failures);
    Lifecycle_Failures_monte(MonteCarlo)=sum(Lifecycle_Failures);
    Total_Costs_monte(MonteCarlo)=sum(Total_Costs_Array);
    Sat_Age_monte(MonteCarlo)=mean(sats.Age(:));
    Max_Extensibility_monte(MonteCarlo)=max(Performance)-1;
    
    
    %Ilities
    No_Recoveries_S4=0;
    No_Recoveries_S4=0;
    Robustness_monte=zeros(1,NumMonteCarlos);
    if ~isempty(S4_Flags{seed})
        for i=1:length(S4_Flags{seed})
            try
                Avg_Start_Performance=mean(Performance(S4_Flags{seed}(i)-3:1:S4_Flags{seed}(i)));
            catch
                Avg_Start_Performance=mean(Performance(1:1:S4_Flags{seed}(i)));
            end
            
            try
                Post_Performance=Performance(S4_Flags{seed}(i):1:S4_Flags{seed}(i)+6);
            catch
                Post_Performance=Performance(S4_Flags{seed}(i):1:end);
            end
            Robustness_S4_monte=horzcat(Robustness_S4_monte,(Avg_Start_Performance-min(Post_Performance)));
            
            Resilience_ShortTerm_S4=horzcat(Resilience_ShortTerm_S4,mean(Post_Performance)/Avg_Start_Performance);
            Improvement_S4=find(Performance(S4_Flags{seed}(i):end)>=(Avg_Start_Performance));
            if ~isempty(Improvement_S4)
                Agility_S4=horzcat(Agility_S4,Improvement_S4(1));
            else
                No_Recoveries_S4=No_Recoveries_S4+1;
            end
            
        end
        
        No_Recoveries_S4_monte=horzcat(No_Recoveries_S4_monte,(No_Recoveries_S4/length(S4_Flags{seed})));
    end
    if ~isempty(S5_Flags{seed})
        for i=1:length(S5_Flags{seed})
            try
                Avg_Start_Performance=mean(Performance(S5_Flags{seed}(i)-3:1:S5_Flags{seed}(i)));
            catch
                Avg_Start_Performance=mean(Performance(1:1:S5_Flags{seed}(i)));
            end
            
            try
                Post_Performance=Performance(S5_Flags{seed}(i):1:S5_Flags{seed}(i)+6);
            catch
                Post_Performance=Performance(S5_Flags{seed}(i):1:end);
            end
            Robustness_S5_monte=horzcat(Robustness_S5_monte,(Avg_Start_Performance-min(Post_Performance)));
        end
        
        Robustness_Max_monte=horzcat(Robustness_Max_monte,min(Robustness(10:end)));
        Resilience_ShortTerm_S5=horzcat(Resilience_ShortTerm_S5,mean(Post_Performance)/Avg_Start_Performance);
        Improvement_S5=find(Performance(S5_Flags{seed}(i):end)>=(Avg_Start_Performance));
        if ~isempty(Improvement_S5)
            Agility_S5=horzcat(Agility_S5,Improvement_S5(1));
        else
            No_Recoveries_S5=No_Recoveries_S5+1;
        end
        
        
        No_Recoveries_S5_monte=horzcat(No_Recoveries_S5_monte,(No_Recoveries_S5/length(S5_Flags{seed})));
    end
    
    %demand_array_all(MonteCarlo,:)=demand_array;
    tEnd=toc(tStart);
    if mod(MonteCarlo,100)==0
        Run_Time=horzcat(Run_Time,tEnd-Run_Time(end));
        fprintf("\n");
        toc(tStart)
        fprintf("Estimated Time Until Completion: %.2f Seconds (%.2f Minutes). \nTotal Monte Carlos Runs = %i / %i (%.2f Percent Complete).", ((NumMonteCarlos-MonteCarlo)/100*mean(Run_Time(2:end))),((NumMonteCarlos-MonteCarlo)/100*mean(Run_Time(2:end)))/60 ,MonteCarlo,NumMonteCarlos,MonteCarlo/NumMonteCarlos*100);
    end
end
%% Analyze Monte Carlo Results

ENPV=mean(NPV_All)

%% "Ilities" Chart
figure()
tiledlayout(4,4)

ax1=nexttile(1);
hold(ax1,'on')
cdfplot(Robustness_S4_monte)
hold(ax1,'off')
title('Target Curve of Robustness to S4 Events Over Solar Cycle (Performance Lost)')

ax2=nexttile(2);
hold(ax2,'on')
cdfplot(Robustness_S5_monte)
hold(ax2,'off')
title('Target Curve of Robustness to S5 Events Over Solar Cycle (Performance Lost)')

ax3=nexttile(3);
hold(ax3,'on')
cdfplot(Robustness_Max_monte)
hold(ax3,'off')
title('Target Curve of Maximum Robustness in Each Simulation (Most Performance Lost)')

ax4=nexttile(4);
hold(ax4,'on')
Resilience_Truncation_S4=sort(Resilience_ShortTerm_S4);
cdfplot(Resilience_Truncation_S4(1:floor(.95*length(Resilience_ShortTerm_S4))))
hold(ax4,'off')
title('Target Curve of Short Term S4 Event Resilience')

ax5=nexttile(5);
hold(ax5,'on')
Resilience_Truncation_S5=sort(Resilience_ShortTerm_S5);
cdfplot(Resilience_Truncation_S5(1:floor(.95*length(Resilience_ShortTerm_S5))))
hold(ax5,'off')
title('Target Curve of Short Term S5 Event Resilience')

ax6=nexttile(6);
hold(ax6,'on')
cdfplot(Agility_S4)
hold(ax6,'off')
title('CDF of Agility for S4 Events (Months to Return to Pre-Event Performance)')

ax7=nexttile(7);
hold(ax7,'on')
cdfplot(Agility_S5)
hold(ax7,'off')
title('CDF of Agility for S5 Events (Months to Return to Pre-Event Performance)')

ax8=nexttile(8);
hold(ax8,'on')
cdfplot(Max_Extensibility_monte)
hold(ax8,'off')
title('Target Curve of Maximum Extensibility')
%% Auxiliary Statistics
figure()
tiledlayout(4,4)

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
histogram(total_cycle_radiation_monte)
title('Histogram of Total Radiation for Solar Cycle')

ax8=nexttile(8);
hold(ax8,'on')
cdfplot(total_cycle_radiation_monte)
hold(ax8,'off')
title('Target Curve of Total Radiation for Solar Cycle')

ax9=nexttile(9);
histogram(Lifecycle_Failures_monte)
title('Histogram of Lifecycle Failures Over Solar Cycle')

ax10=nexttile(10);
hold(ax10,'on')
cdfplot(Lifecycle_Failures_monte)
hold(ax10,'off')
title('Target Curve of Lifecycle Failures Over Solar Cycle')

ax11=nexttile(11);
histogram(Total_Costs_monte)
title('Histogram of Total Costs Over Solar Cycle')

ax12=nexttile(12);
hold(ax12,'on')
cdfplot(Total_Costs_monte)
hold(ax12,'off')
title('Target Curve of Total Cost Over Solar Cycle')

ax13=nexttile(13);
histogram(Sat_Age_monte)
title('Histogram of Average Satellite Age at Failure Over Solar Cycle')

ax14=nexttile(14);
hold(ax14,'on')
cdfplot(Sat_Age_monte)
hold(ax14,'off')
title('Target Curve of Average Satellite Age at Failure Over Solar Cycle')

ax15=nexttile(15);
hold(ax15,'on')
cdfplot(No_Recoveries_S4_monte)
hold(ax15,'off')
title('CDF of Percent of S4 Events That Never Return to Pre-Event Performance')

ax16=nexttile(16);
hold(ax6,'on')
cdfplot(No_Recoveries_S5_monte)
hold(ax6,'off')
title('CDF of Percent of S5 Events That Never Return to Pre-Event Performance')
%% Sample Run Print-Out

figure()
tiledlayout(2,2)
ax1=nexttile(1);
plot(OpOrbitPop)
title('Operational Orbit Population Over Solar Cycle')

ax2=nexttile(2);
hold(ax2,'on')
plot(Performance)
hold(ax2,'off')
title('Overall Performance Over Solar Cycle')

ax3=nexttile(3);
plot(AreaOfCoverage)
title('Area of Coverage Over Solar Cycle')

ax4=nexttile(4);
hold(ax4,'on')
plot(Revisit)
hold(ax4,'off')
title('Revisit Rates (Hours) Over Solar Cycle')

mean(Performance_All)
mean(Resilience_Truncation_S4(1:floor(.95*length(Resilience_ShortTerm_S4))))
mean(Resilience_Truncation_S5(1:floor(.95*length(Resilience_ShortTerm_S5))))