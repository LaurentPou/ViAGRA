function [Y20,dY20,d2Y20]=Y20(theta)

Y20=5/16/pi*(3*cos(lat).^2-1);
dY20=-3/2*sqrt(5/pi)*cos(theta).*sin(theta);
d2Y20=-3/2*sqrt(5/pi)*(2.*cos(theta).^2-1);