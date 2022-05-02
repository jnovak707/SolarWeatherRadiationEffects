%% Plot SpaceX constellation
clc; clear all;
addpath(genpath('/Users/jo21938/Documents/GitHub/space-environment/'))

% include azel plot for a site
azel_plot=false;

figure(1); clf; 
if azel_plot; subplot(1,2,1); hold on; end
h=plotearth('MapType','bluemarble','SampleStep',5,'AxeHandle',gca);

t_now=datetime('23-May-2019 11:30:00','TimeZone','UTC')-tzoffset(datetime('today','TimeZone','America/New_York'));%datetime('2-June-2019 00:00:00','TimeZone','UTC');%
sats=sgp4_file(t_now,'spaceX.tle');

satListstr = [];
for ii=1:length(sats) 
    ecf=eci2ecf(t_now,sats(ii).pos');
    sats(ii).ecf=ecf;
    allECEF(ii,:) = ecf;    
end
% Site Location
site.lat = 42.3581; site.lon = -71.0636;
site.h=0; 

% Site Vector
wgs84 = wgs84Ellipsoid('kilometers');
[site.x,site.y,site.z] = geodetic2ecef(wgs84,site.lat,site.lon,site.h);

plot3(site.x,site.y,site.z,'ro')
plot3(allECEF(:,1),allECEF(:,2),allECEF(:,3),'.g')

view([site.x,site.y,site.z])
camzoom(1.5)

% Plot visibility
% Kinda sorta topocentric Az./El. angular measurments
% [az,el,rho]=Observe(repmat([site.x,site.y,site.z],length(sats),1)',allECEF');

% better matlab
[az,el,rho]=ecef2aer(allECEF(:,1),allECEF(:,2),allECEF(:,3),site.lat,site.lon,site.h,wgs84);
above_horizon=find(el>0);

if azel_plot
figure(1); 
subplot(1,2,2)
polarplot(deg2rad(az(above_horizon)),90-el(above_horizon),'o')
rticks([0:10:90])
rticklabels(split(num2str(90:-10:0)))
rlim([0 90])
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise')
% 39 25
end
sgtitle(['Line of Site to Starlink at ',datestr(t_now+tzoffset(datetime('today','TimeZone','America/New_York'))),' ET'])