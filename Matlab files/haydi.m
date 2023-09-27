function [AZ]=haydi(lat,lon,A1)

%figure
surfc(lon*180/pi,lat*180/pi,A1');
SHADING INTERP 
cb=colorbar('v');
%set(get(cb,'title'),'string','Titan stress');
ylabel('latitude')
xlabel('longitude')
ylim([0 180])
xlim([0 360])
zlim([-40 40])