% Convert ECI (J2000) to ECEF
% SatECI: with velocity: n x 6, X, X-dot, Y, Y-dot, Z, Z-dot
%   without velocity: n x 3, X, Y, Z
% mjd: n x 1
%
% Output:
%   SatECEF: n x 6, X, X-dot, Y, Y-dot, Z, Z-dot
function SatEcef=eci2ecf(t,SatECI)

% inherently use matlab datetime, convert to mjd here
mjd=juliandate(t)-2400000.5;

% this is using WGS-84, be congnizant of that
earth=wgs84Ellipsoid('kilometers');

if(size(SatECI,2)==3)
    x=1;
    y=2;
    z=3;
    calc_vel=false;
else
    x=1;
    y=3;
    z=5;
    calc_vel=true;
end

% convert to geodetic to do rotation shift and convert back
SatGeodetic=zeros(length(mjd),3);
[SatGeodetic(:,1),SatGeodetic(:,2),SatGeodetic(:,3)]=...
            ecef2geodetic(SatECI(:,x),SatECI(:,y),SatECI(:,z),earth);

SatGeodetic(:,[1,2])=SatGeodetic(:,[1,2])*180./pi;
change=mod(SatGeodetic(:,2)-(280.4606+360.9856473*(mjd-51544.5)),360)-SatGeodetic(:,2);
SatGeodetic(:,2)=mod(SatGeodetic(:,2)-(280.4606+360.9856473*(mjd-51544.5)),360);
RotateAngle=change;

if(calc_vel)
    SatEcef=zeros(size(SatGeodetic,1),6);
    [SatEcef(:,x),SatEcef(:,y),SatEcef(:,z)]=geodetic2ecef(SatGeodetic(:,1)*pi/180,SatGeodetic(:,2)*pi/180,SatGeodetic(:,3),earth);
    for i=1:size(SatEcef,1)
        Rot=[cosd(RotateAngle(i)) -sind(RotateAngle(i)) 0;
            sind(RotateAngle(i)) cosd(RotateAngle(i)) 0;
            0 0 1];
        SatEcef(i,[2,4,6])=SatECI(i,[2,4,6])*Rot';
    end
else
    SatEcef=zeros(size(SatGeodetic,1),3);
    [SatEcef(:,x),SatEcef(:,y),SatEcef(:,z)]=geodetic2ecef(SatGeodetic(:,1)*pi/180,SatGeodetic(:,2)*pi/180,SatGeodetic(:,3),earth);
end