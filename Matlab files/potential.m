function [Pot,lat,lon,Thn]=decomp_pot(yi,rp,LN)
Prot=15.945; %jr
n=4.56E-6; %m.s moyen
R=2575000;
e=0.0292;


%yi le long de r ;cf alterman; matrice nb_couches * 6
%rp rayons correspondant aux yi 
%LN nombre de love : h,l,k de degre 2

Tf=Prot*24*60*60;
Nlat=18*2;
Nlon=36*2;

Dlat=180/Nlat;
Nt=100;

lat(1:Nlat+1)=(0:Dlat:180)*pi/180;
lon(1:Nlon+1)=(0:Dlat:360)*pi/180;

t=0;
for tt=1:Nt
    Time(tt)=t;
    for i=1:Nlon+1
        for j=1:Nlat+1
           Pot(tt,i,j)= (9/8+3/2*cos(lat(j)))*cos(2*lon(i)-n*t)+ (9/8-3/2*cos(lat(j)))*cos(2*lon(i)+n*t)+3/4*cos(n*t);
           dPotdlat(tt,i,j)=-3/2*sin(lat(j))*cos(2*lon(i)-n*t)+3/2*sin(lat(j))*cos(2*lon(i)+n*t); %ajoute derivé seconde lat
           dPotdlon(tt,i,j)=-(9/8+3/2*cos(lat(j)))*2*sin(2*lon(i)-n*t)-(9/8-3/2*cos(lat(j)))*2*sin(2*lon(i)+n*t);%ajoute derivé seconde lon
       end
    end
    t=t+Tf/(Nt-1);
    
end
       Th=Time/24/60/60;
       Thn=Th/Prot;
       
       lon=lon*180/pi;
       lat=lat*180/pi;
      
       figure
      %colorbar('vert')
      grid on
      title('Potential')
      zlabel ('Longitude','FontSize',16);
      ylabel ('\sigma','FontSize',16);
      %set(gca,'yscale','log')
      xlabel ('time','FontSize',16)
      contourf(lon,Thn,Pot(:,:,1))
      colorbar('vert')
      grid on
      
      figure
      surf(lon,Thn,Pot(:,:,1)) 
      xlim([0 360])
      ylim([0:1])
      SHADING INTERP 
    
      hold on
      contourf(lon,Thn,Pot(:,:,1))
      xlabel ('Longitude (Degree)','FontSize',16);
      ylabel ('Time/Period','FontSize',16);
      zlabel ('\psi/(nR)^2e+3)')
      colorbar('vert')
      grid on
      
      figure
      surfc(lon,Thn,Pot(:,:,1)) 
      xlim([0 360])
      ylim([0:1])
      SHADING INTERP 
      xlabel ('Longitude (Degree)','FontSize',16);
      ylabel ('Time/Period','FontSize',16);
      zlabel ('\psi/(nR)^2e+3)')
      colorbar('vert')
      grid on
      
      
      figure
      meshc(lon,Thn,Pot(:,:,1))
      SHADING INTERP
      xlabel ('Longitude (Degree)','FontSize',16);
      ylabel ('Time/Period','FontSize',16);
      zlabel ('\psi/(nR)^2e+3)')
       xlim([0 360])
      ylim([0:1])
      SHADING INTERP 
    
      colorbar('vert')
      
      
for i=1:100
    i=i-mod(i,1)
contourf(lat,lon,squeeze(Pot(i,:,:)));
xlabel('lat')
ylabel('lon')
ylim([0 360])
view(2);
frames(:,i)=getframe;
end
map=colormap;
mpgwrite(frames,map,'Titan_tide2')

