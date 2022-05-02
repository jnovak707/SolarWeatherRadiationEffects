addpath('D:\Materials\')
load('D:\Materials\H_He_SPE_Fluxes_noGeoShield');
load('D:\Materials\All_SEP_Events_noGeoShield.mat')
model(4:end)=[];
Egrid2=model(1).Egrid;
%% Propagate Satellite to get all Lat, Lon, and Alt
 files='SkySatShieldedTLE.txt';    
    t_start=datetime('2021-07-12 12:00:00','TimeZone','UTC');
    t_interval=minutes(100);
    t_step=seconds(10);
    t_now=[t_start:t_step:(t_start+(t_interval))];
     
    for iii=1:length(t_now)
        
        sats=sgp4_file(t_now(iii),files);

        allECEF=[];
        for ii=1:length(sats)
            ecf=eci2ecf(t_now(iii),sats(ii).pos');
            sats(ii).ecf=ecf;
            allECEF(ii,:) = ecf;
            [lla_params(ii,1,iii),lla_params(ii,2,iii),lla_params(ii,3,iii),lla_params(ii,4,iii)] = ijk2ll(ecf);
            lla_params(ii,1:3,iii)=rad2deg([lla_params(ii,1:3,iii)]);
        end
        location.alt=squeeze(lla_params(:,4,:));location.lat=squeeze(lla_params(:,2,:));location.lon=squeeze(lla_params(:,3,:));
        
        %Plot the satellite track
        %plot3(allECEF(:,1),allECEF(:,2),allECEF(:,3),'.','linewidth',1,"Color",'#FF69B4')
    end
    
    for qq=1:length(sats)
        location_final{qq}.alt=location.alt(qq,:);
        location_final{qq}.lat=location.lat(qq,:);
        location_final{qq}.lon=location.lon(qq,:);
    end
    
%% Load the Unpropagated NOAA_Storms from Storm_Isolator
load('D:\Materials\Storms_Unpropagated')

%% Apply Geomagnetic Shielding 
%  This gets applied to every 5 minute differential flux in all of
%  storm.flux which means about ~96k runs
alt=location.alt;
lat=location.lat;
lon=location.lon;
%Load all the relevant shielding datasets
%Load Kp6 Model
L_shells=load('Lshell_T89_kp6.mat');
L_shells=L_shells.LmL_m3;
grid_pnts = load(strcat('lat_long_alt.mat'));
ALT_L = grid_pnts.ALT_L;
LAT_L = grid_pnts.LAT_L;
LON_L = grid_pnts.LON_L;
Ai = load(strcat('elements_A.mat'));   %atomic masses (amu)
lon = wrapTo180(lon);

%L Shells
Ls = interp3(LAT_L,LON_L,ALT_L,L_shells,lat,lon,alt);
reg1 = Ls<=2.9;
reg2 = Ls>2.9 & Ls<=5;

R_cv = zeros(size(Ls));

R_cv(reg1) =  14.5./Ls(reg1).^2;

R_cv(reg2) = 15.062./Ls(reg2).^2-0.363; % Sampex Correction

% calculate cuttoff energy for each element
m_c2 = .931494;    %rest energy of 1 amu (GeV)
Zi = 1;%:92;         %atomic number, uncomment to extend to heavy ions
Ai = Ai.Ai';        
Ai=Ai(1);           %comment this to extend to all heavy ions

for i=1:length(storm)
    Flux=[];
    for ii=1:height(storm(i).flux)
E_cv = arrayfun(@(r) 1000*m_c2*(sqrt((r*Zi./(Ai*m_c2)).^2+1)-1),R_cv,'UniformOutput',false);

% adjust flux given cutoff energy for each element
d_f = cellfun(@(ecv,h) (Egrid2>ecv).*(storm(i).flux(ii,:)'*0.5*(1+sqrt((6378.1+h)^2-(6378.1)^2)/(6378.1+h))),E_cv,num2cell(alt),'UniformOutput',false);
d_f=reshape(cell2mat(d_f),length(Egrid2),length(alt));
d_f(d_f==0)=NaN;
OrbitAveragedFlux=mean(d_f,2,'omitnan');
Flux=horzcat(Flux,OrbitAveragedFlux);
    end
    storm(i).ShieldedFluxes=Flux;
    i
end

for i=1:length(storm)
    storm(i).ShieldedFluxes(:,isnan(sum(storm(i).ShieldedFluxes)))=[];%remove columns with NaN values
end
save('Propagated_Storms_SkySatShielded.mat','storm')