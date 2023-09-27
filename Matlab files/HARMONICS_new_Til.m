%function [disp,lat,lon,time]=PLOTS(A1,per,kpot,Vpot)


% POTENTIAL
nn=sqrt((G*(M1+m2)/(a^3)));
epot=(nn*ref(1))^2*e; %potentiel subit par le coprs ?tudi? (eccentricity)
opot=(nn*ref(1))^2*obl; % (obliquity)
    
% Values Initialisation
Time=zeros(Nt);

D2YdT=zeros(Nt,Nlon,Nlat);

ur=zeros(Nt,Nlon,Nlat);
ut=zeros(Nt,Nlon,Nlat);
up=zeros(Nt,Nlon,Nlat);

Pot=zeros(Nt,Nlon,Nlat);
disp2=zeros(Nt,Nlon,Nlat);
deltag=zeros(Nt,Nlon,Nlat);
dPot=zeros(Nt,Nlon,Nlat);
deltilt=zeros(Nt,Nlon,Nlat);
    
% Time loop
Nt=72;
Tf=per;
t=0;
for tt=1:Nt
    Time(tt)=t;
    cos_nt=cos(nn*t);%pulsation
    sin_nt=sin(nn*t);
    A20=-3*sqrt(pi/5)*cos_nt*epot;
    A21=2*sqrt(6*pi/5)*sin_nt*opot;
    %A22=3*sqrt(6*pi/5)*cos_nt*epot;
    %A22n=4*sqrt(6*pi/5)*sin_nt*epot;
    
    A22=epot*1/2*sqrt(6*pi/5)*(3*cos_nt+4/i*sin_nt);
    A22n=epot*1/2*sqrt(6*pi/5)*(3*cos_nt-4/i*sin_nt);
    
    for ix=1:Nlon
        
        
        for j=1:Nlat
            sinlat=sin(lat(j));
            coslat=cos(lat(j));
            
            P20=1/2*(3*(cos(lat(j)))^2-1);
            P21=-3*sin(lat(j))*cos(lat(j));
            P22=3*(sin(lat(j)))^2;
            
            %Anm
            
            
            D2YdT(tt,ix,j)=(A20*d2Y20dt(ix,j)+A22*d2Y22dt(ix,j)+A22n*d2Y22ndt(ix,j));
            
            %displacement
            ur(tt,ix,j)=(yi(1,end)*(A20*Y20(ix,j)+A22*Y22(ix,j)+A22n*Y22n(ix,j)));
            %ur_21(tt,i,j)=real(yi(1,end))*(A20*Y20(i,j)+A21*Y21(i,j)+A22*Y22(i,j)+A22n*Y22n(i,j));
            ut(tt,ix,j)=(yi(3,end)*(A20*dY20dt(ix,j)+A22*dY22dt(ix,j)+A22n*dY22ndt(ix,j)));
            up(tt,ix,j)=(yi(3,end)/sinlat*(A20*dY20dp(ix,j)+A22*dY22dp(ix,j)+A22n*dY22ndp(ix,j)));
            
            
            % POTENTIAL
            
            
            Pot(tt,ix,j)= epot*(-3/2*P20*cos_nt+1/4*P22*(3*cos(2*lon(ix))*cos_nt+4*sin(2*lon(ix))*sin_nt));
            %Pot_St(tt,i,j)= kpot/e*(1+3*e*cos_nt)*(-1/2*P20+1/4*P22*(cos(2*lon(i))+4*e*sin(2*lon(i))*sin_nt));
            %Pot1_1(tt,i,j)= kpot*(-3/2*P20*cos_nt);
            %Pot1_2(tt,i,j)= kpot*(1/4*P22*3*cos(2*lon(i))*cos_nt);
            %Pot1_3(tt,i,j)= kpot*(1/4*P22*4*sin(2*lon(i))*sin_nt);
            %Pot1_4(tt,i,j)=opot*P21*obl*cos(lon(i))*sin_nt;
            %Pot_rad(tt,i,j)=kpot*3/2*(3*(sin(lat(j)))^2*(cos(lon(i)))^2-1)*cos_nt;
            %Pot_lib(tt,i,j)=kpot*3*(sin(lat(j)))^2*sin(2*lon(i))*sin_nt;
            %Pot_sears(tt,i,j)=kpot*(3*P20*cos_nt+2*P21*cos(lon(i))*sin_nt);
            %PotS1_1(tt,i,j)=kpot*(9/8+3/2*cos(lon(i))*cos(2*lat(j)-n*t));
            %PotS1_2(tt,i,j)=kpot*(9/8-3/2*cos(lon(i))*cos(2*lat(j)+n*t));
            %PotS1_3(tt,i,j)=kpot*3/4*cos(n*t);
            
            disp2(tt,ix,j)=(LN(1))*Pot(tt,ix,j)/ref(4);
            factg= 1 + LN(1) - 3/2* LN(3);
            deltag(tt,ix,j) = real(factg) * Pot(tt,ix,j) * ref(4) /ref(1);
            
            %deltilt(tt,ix,j)=- real( 1 + LN(3) -  LN(1) ) * (2*coslat *P20- 2*P21) * Pot(tt,ix,j) / ref(4) /ref(1)/sinlat;
            %derivative wrt lat
            dP20= (-3)*cos(lat(j))*sin(lat(j));
            dP21=3*sin(lat(j))^2 - 3*cos(lat(j))^2;
            dP22= 6*cos(lat(j))*sin(lat(j));
            dPot(tt,ix,j)= epot*(-3/2*dP20*cos_nt+1/4*P22*(3*cos(2*lon(ix))*cos_nt+4*sin(2*lon(ix))*sin_nt));
            deltilt(tt,ix,j)=- real( 1 + LN(3) -  LN(1) ) /ref(4) /ref(1) * dPot(tt,ix,j);
        end
    end

    t=t+Tf/(Nt-1);
    
end

%disp=disp*kpot/Vpot;
d1=max(max(max(ur)));
d1c=max(max(max(disp2)));
disp(strcat('u_r =',num2str(d1)))
disp(strcat('u_r2 =',num2str(d1c)))
d2=max(max(max(ut)));
d3=max(max(max(up)));
%d12=max(max(max(ur_21)))
%disp(strcat('u_r =',num2str(d12)))
%d2=max(max(max(disp2)))


%Radial Surface Displacements (m)
if bool_surface_disp_maps
    Cizdir9(lon,lat,real(ur),Nt,Time,per);
end
%h=colorbar('v');

%figure
%Cizdir(lon,lat,Pot_St);

 
figure
title('Radial & Lateral Surface Displacements (in meters)','FontSize',18)
subplot(311)
Cizdir1(lon,lat,real(ur),Nt,Time,per,1);
% title(' u_r','FontSize',18)
% h=colorbar('v');
subplot(312)
Cizdir1(lon,lat,real(ut),Nt,Time,per,1);
% title('u_\theta','FontSize',18)
% h=colorbar('v');
subplot(313)
Cizdir1(lon,lat,real(up),Nt,Time,per,1);
% title('u_\phi','FontSize',18)
% h=colorbar('v');
%subplot(224)
%Cizdir1(lon,lat,Pot,Nt,Time,per,1);
%title('\Phi','FontSize',18)
%h=colorbar('h');

% Equator
figure
plot(Time/(24*60*60),real(ur(:,1,19)), '-k')
xlim([min(Time/(24*60*60)) max(Time/(24*60*60))])
xlabel('Time (Days)')
ylabel('Surface Displacement (m)')

% Equator
figure
plot(Time/(24*60*60),real(deltilt(:,1,19)), '-k')
xlim([min(Time/(24*60*60)) max(Time/(24*60*60))])
xlabel('Time (Days)')
ylabel('Tilt (rad)')

% Equator
figure
plot(Time/(24*60*60),real(deltag(:,1,19)*1e5), '-k')
xlim([min(Time/(24*60*60)) max(Time/(24*60*60))])
xlabel('Time (Days)')
ylabel('g (mgal)')

figure
plot(Time/(24*60*60),real(deltag(:,1,10)), '-k')
xlim([min(Time/(24*60*60)) max(Time/(24*60*60))])
xlabel('Time (Days)')
ylabel('g (m/s2)')


surface_displacemement.radial=real(ur);
surface_displacemement.theta=real(ut);
surface_displacemement.phi=real(ut);
surface_displacemement.lon=lon;
surface_displacemement.lat=lat;
surface_displacemement.Nt=Nt;
surface_displacemement.Time=Time;
surface_displacemement.per=per;
delta_g=real(deltag);
delta_tilt=real(deltilt);

nur=ur;
ndeltilt=deltilt;
ndeltag=deltag;


save('dispnoocean.mat', 'nur', 'lon', 'lat' , 'Time')
save('Tiltnoocean.mat', 'ndeltilt', 'lat' , 'Time')
save('gravnoocean.mat', 'ndeltag', 'lon', 'lat' , 'Time')


%save('Europadispoft.mat', 'ur', 'lon', 'lat' , 'Time')
%save('EuropaTiltsoft.mat', 'deltilt', 'lat' , 'Time')
%save('Europagravsoft.mat', 'deltag', 'lon', 'lat' , 'Time')

%save('dispocean.mat', 'ur', 'lon', 'lat' , 'Time')
%save('Tiltocean.mat', 'deltilt', 'lat' , 'Time')
%save('gravocean.mat', 'deltag', 'lon', 'lat' , 'Time')