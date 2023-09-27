
N_lat=36;
N_lon=48;
Frac_lat=5; %1/fraclat*lat(2) is the first point of the mesh to avoid 0 in denominator
Nt=72;

% POTENTIAL
nn=sqrt((G*(M1+m2)/(a^3)));
epot=(nn*ref(1))^2*e; %potentiel subit par le coprs étudié (eccentricity)
opot=(nn*ref(1))^2*obl; % (obliquity)

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

clear i
for ix=1:Nlon
    
    sinlon=sin(lon(ix)); 
    coslon=cos(lon(ix));
    sin2lon=sin(2*lon(ix)); 
    cos2lon=cos(2*lon(ix)); 
    Ei=exp(2*lon(ix)*i);
    Ein=conj(Ei);
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
            % Y22
            
            
            Y22(ix,j)=(sqrt(15/32/pi)*sinlat^2*Ei);
            dY22dt(ix,j)=(1/4*sqrt(30/pi)*sinlat*Ei*coslat);
            d2Y22dt(ix,j)=(1/4*sqrt(30/pi)*Ei*(2*coslat^2-1));
            dY22dp(ix,j)=1/4*sqrt(30/pi)*sinlat^2*i*Ei;
            d2Y22dp(ix,j)=-1/2*sqrt(30/pi)*sinlat^2*Ei;
            d2Y22dpt(ix,j)=1/2*sqrt(30/pi)*i*Ei*coslat*sinlat;
            %Y22n
            Y22n(ix,j)=(sqrt(15/32/pi)*sinlat^2*Ein);
            dY22ndt(ix,j)=1/4*sqrt(30/pi)*sinlat*Ein*coslat;
            d2Y22ndt(ix,j)=1/4*sqrt(30/pi)*Ein*(2*coslat^2-1);
            dY22ndp(ix,j)=-1/4*i*sqrt(30/pi)*sinlat^2*Ein;
            d2Y22ndp(ix,j)=-1/2*sqrt(30/pi)*sinlat^2*Ein;
            d2Y22ndpt(ix,j)=-1/2*i*sqrt(30/pi)*Ein*coslat*sinlat;
            
      end
    end
    
    
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
        end
    end
    t=t+Tf/(Nt-1);
end  

surface_displacemement.radial=real(ur);
surface_displacemement.theta=real(ut);
surface_displacemement.phi=real(ut);
surface_displacemement.lon=lon;
surface_displacemement.lat=lat;
surface_displacemement.Nt=Nt;
surface_displacemement.time=Time;
surface_displacemement.per=per;
surface_displacemement.pot=Pot;
