% MESH
Dlat=180/N_lat;
Dlon=360/N_lon;
lat=(0:Dlat:180)*pi/180;
lon=(0:Dlon:360)*pi/180;
Nlat=length(lat);
Nmiddle=(length(lat)+1)/2;
Nlon=length(lon);
lat(1)=lat(2)/Frac_lat;
lat(end)=pi-lat(1);
lat(Nmiddle)=lat(Nmiddle)+lat(1);
%Nt = Nt+1;

% Harmonics Initialisation - not defined because sparse thus faster

% Y20=zeros(Nlon,Nlat);
% dY20dt=zeros(Nlon,Nlat);
% d2Y20dt=zeros(Nlon,Nlat);
% dY20dp=zeros(Nlon,Nlat);
% d2Y20dp=zeros(Nlon,Nlat);
% d2Y20dpt=zeros(Nlon,Nlat);
% 
% Y21=zeros(Nlon,Nlat);
% 
% Y22=zeros(Nlon,Nlat);
% dY22dt=zeros(Nlon,Nlat);
% d2Y22dt=zeros(Nlon,Nlat);
% dY22dp=zeros(Nlon,Nlat);
% d2Y22dp=zeros(Nlon,Nlat);
% d2Y22dpt=zeros(Nlon,Nlat);
% 
% Y22n=zeros(Nlon,Nlat);
% dY22ndt=zeros(Nlon,Nlat);
% d2Y22ndt=zeros(Nlon,Nlat);
% dY22ndp=zeros(Nlon,Nlat);
% d2Y22ndp=zeros(Nlon,Nlat);
% d2Y22ndpt=zeros(Nlon,Nlat);

clear i
for ix=1:Nlon
    
    sinlon=sin(lon(ix));
    coslon=cos(lon(ix));
    sin2lon=sin(2*lon(ix));
    cos2lon=cos(2*lon(ix));
    Ei=exp(lon(ix)*1i);
    E2i=exp(2*lon(ix)*1i);
    Ein=conj(Ei);
    E2in=conj(E2i);
    %Ein=exp(-2*lon(ix)*i);
    
    for j=1:Nlat
        sinlat=sin(lat(j));
        coslat=cos(lat(j));
        P20=1/2*(3*(cos(lat(j)))^2-1);
        P21=-3*sin(lat(j))*cos(lat(j));
        P22=3*(sin(lat(j)))^2;
        
        
        
        %Y20
        Y20(ix,j)=sqrt(5/16/pi)*(3*coslat^2-1);
        dY20dt(ix,j)=-3/2*sqrt(5/pi)*coslat*sinlat;
        d2Y20dt(ix,j)=-3/2*sqrt(5/pi)*(2*coslat^2-1);
        dY20dp(ix,j)=0;
        d2Y20dp(ix,j)=0;
        d2Y20dpt(ix,j)=0;
        %Y21
        Y21(ix,j)=0.5*sqrt(5/6/pi)*(-3)*sinlat*coslat*coslon;
        %Y21(ix,j)=-sqrt(15/8/pi)*sinlat*coslat*Ei;
        
        % Y22
        Y22(ix,j)=(sqrt(15/32/pi)*sinlat^2*E2i);
        dY22dt(ix,j)=(1/4*sqrt(30/pi)*sinlat*E2i*coslat);
        d2Y22dt(ix,j)=(1/4*sqrt(30/pi)*E2i*(2*coslat^2-1));
        dY22dp(ix,j)=1/4*sqrt(30/pi)*sinlat^2*1i*E2i;
        d2Y22dp(ix,j)=-1/2*sqrt(30/pi)*sinlat^2*E2i;
        d2Y22dpt(ix,j)=1/2*sqrt(30/pi)*1i*E2i*coslat*sinlat;
        %Y22n
        Y22n(ix,j)=(sqrt(15/32/pi)*sinlat^2*E2in);
        dY22ndt(ix,j)=1/4*sqrt(30/pi)*sinlat*E2in*coslat;
        d2Y22ndt(ix,j)=1/4*sqrt(30/pi)*E2in*(2*coslat^2-1);
        dY22ndp(ix,j)=-1/4*1i*sqrt(30/pi)*sinlat^2*E2in;
        d2Y22ndp(ix,j)=-1/2*sqrt(30/pi)*sinlat^2*E2in;
        d2Y22ndpt(ix,j)=-1/2*1i*sqrt(30/pi)*E2in*coslat*sinlat;

%         % Y22
%         Y22(ix,j)=(sqrt(15/32/pi)*sinlat^2*Ei);
%         dY22dt(ix,j)=(1/4*sqrt(30/pi)*sinlat*Ei*coslat);
%         d2Y22dt(ix,j)=(1/4*sqrt(30/pi)*Ei*(2*coslat^2-1));
%         dY22dp(ix,j)=1/4*sqrt(30/pi)*sinlat^2*1i*Ei;
%         d2Y22dp(ix,j)=-1/2*sqrt(30/pi)*sinlat^2*Ei;
%         d2Y22dpt(ix,j)=1/2*sqrt(30/pi)*1i*Ei*coslat*sinlat;
%         %Y22n
%         Y22n(ix,j)=(sqrt(15/32/pi)*sinlat^2*Ein);
%         dY22ndt(ix,j)=1/4*sqrt(30/pi)*sinlat*Ein*coslat;
%         d2Y22ndt(ix,j)=1/4*sqrt(30/pi)*Ein*(2*coslat^2-1);
%         dY22ndp(ix,j)=-1/4*1i*sqrt(30/pi)*sinlat^2*Ein;
%         d2Y22ndp(ix,j)=-1/2*sqrt(30/pi)*sinlat^2*Ein;
%         d2Y22ndpt(ix,j)=-1/2*1i*sqrt(30/pi)*Ein*coslat*sinlat;
        
        %Y22(ix,j)=sqrt(15/32/pi)*sinlat^2*cos2lon;
        %dY22dt(ix,j)=1/4*sqrt(30/pi)*sinlat*cos2lon*coslat;
        %d2Y22dt(ix,j)=1/4*sqrt(30/pi)*cos2lon*(2*coslat^2-1);
        %dY22dp(ix,j)=1/4*sqrt(30/pi)*sinlat^2*sin2lon;
        %d2Y22dp(ix,j)=1/2*sqrt(30/pi)*cos2lon*sinlat^2;
        %d2Y22dpt(ix,j)=-1/2*sqrt(30/pi)*sin2lon*coslat*sinlat;
        %Y22n
        %Y22n(ix,j)=sqrt(15/32/pi)*sinlat^2*sin2lon;
        %dY22ndt(ix,j)=1/4*sqrt(30/pi)*sinlat*sin2lon*coslat;
        %d2Y22ndt(ix,j)=1/4*sqrt(30/pi)*sin2lon*(2*coslat^2-1);
        %dY22ndp(ix,j)=-1/4*sqrt(30/pi)*sinlat^2*cos2lon;
        %d2Y22ndp(ix,j)=1/2*sqrt(30/pi)*cos2lon*sinlat^2;
        %d2Y22ndpt(ix,j)=1/2*sqrt(30/pi)*cos2lon*coslat*sinlat;
    end
end
    
% TEST IF IT IS OKEY
%---------------------------
%Y22=real(Y22);
%dY22dt=real(dY22dt);
%d2Y22dt=real(d2Y22dt);
%dY22dp=real(dY22dp);
%d2Y22dp=real(d2Y22dp);
%d2Y22dpt=real(d2Y22dpt);

%Y22n=real(Y22n);
%dY22ndt=real(dY22ndt);
%d2Y22ndt=real(d2Y22ndt);
%dY22ndp=real(dY22ndp);
%d2Y22ndp=real(d2Y22ndp);
%d2Y22ndpt=real(d2Y22ndpt);
%---------------------------


