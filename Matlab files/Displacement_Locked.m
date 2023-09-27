%function [disp,lat,lon,time]=PLOTS(A1,per,kpot,Vpot)


%% POTENTIAL
nn=sqrt((G*(M1+m2)/(a^3)));%2*pi/per;

%% Constant terms (order 0 in e and o), linked to V2200 and V2010

%% Variable terms (order 1 in e and o, linked to V2201, V220-1, V201-1, V2011)
%epot=(nn*ref(1))^2*e; %potentiel subit par le corps ?tudi? (eccentricity)
epot = e*((G*M1*(ref(1))^2)/(a^3)); %
%opot=(nn*ref(1))^2*obl; % (obliquity)
opot = obl*((G*M1*(ref(1))^2)/(a^3));

if bool_retrograde
    nn = -nn;
end


    
% Values Initialisation - not defined because sparse thus faster
% Time=zeros(1,Nt);
% 
% D2YdT=zeros(Nt,Nlon,Nlat);
% 
% ur=zeros(Nt,Nlon,Nlat);
% ut=zeros(Nt,Nlon,Nlat);
% up=zeros(Nt,Nlon,Nlat);
% 
% Pot=zeros(Nt,Nlon,Nlat);
% disp2=zeros(Nt,Nlon,Nlat);
% deltag=zeros(Nt,Nlon,Nlat);
% dPot=zeros(Nt,Nlon,Nlat);
% deltilt=zeros(Nt,Nlon,Nlat);


% Time loop
Tf=per;
t=0;
for tt=1:Nt+1
    Time(tt)=t;
    
    cos_nt=cos(nn*t);%pulsation
    sin_nt=sin(nn*t);
    A20=-3*sqrt(pi/5)*cos_nt*epot;
    A21=2*sqrt(6*pi/5)*sin_nt*opot;
    %A22=3*sqrt(6*pi/5)*cos_nt*epot;
    %A22n=4*sqrt(6*pi/5)*sin_nt*epot;
    
    A22=epot*1/2*sqrt(6*pi/5)*(3*cos_nt+4/1i*sin_nt);
    A22n=epot*1/2*sqrt(6*pi/5)*(3*cos_nt-4/1i*sin_nt);
   
    
    dA20dt = 3*nn*sqrt(pi/5)*sin_nt*epot;
    dA21dt = 2*nn*sqrt(6*pi/5)*cos_nt*opot;
    dA22dt = epot*nn/2*sqrt(6*pi/5)*(-3*sin_nt+4/1i*cos_nt);
    dA22ndt = epot*nn/2*sqrt(6*pi/5)*(-3*sin_nt-4/1i*cos_nt);
    
    d2A20dt = 3*nn^2*sqrt(pi/5)*cos_nt*epot;
    d2A21dt = -2*nn^2*sqrt(6*pi/5)*sin_nt*opot;
    d2A22dt = -epot*nn^2/2*sqrt(6*pi/5)*(3*cos_nt+4/1i*sin_nt);
    d2A22ndt = epot*nn^2/2*sqrt(6*pi/5)*(-3*cos_nt+4/1i*sin_nt);
    
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
            
            
            vr(tt,ix,j)=(yi(1,end)*(dA20dt*Y20(ix,j)+dA22dt*Y22(ix,j)+dA22ndt*Y22n(ix,j)));
            vt(tt,ix,j)=(yi(3,end)*(dA20dt*dY20dt(ix,j)+dA22dt*dY22dt(ix,j)+dA22ndt*dY22ndt(ix,j)));
            
            
            ar(tt,ix,j)=(yi(1,end)*(d2A20dt*Y20(ix,j)+d2A22dt*Y22(ix,j)+d2A22ndt*Y22n(ix,j)));
            at(tt,ix,j)=(yi(3,end)*(d2A20dt*dY20dt(ix,j)+d2A22dt*dY22dt(ix,j)+d2A22ndt*dY22ndt(ix,j)));
            
            
            if bool_greff_lefftz
                up(tt,ix,j)=(yi(7,end)/sinlat*(A20*dY20dp(ix,j)+A22*dY22dp(ix,j)+A22n*dY22ndp(ix,j)));
                vp(tt,ix,j)=(yi(7,end)/sinlat*(dA20dt*dY20dp(ix,j)+dA22dt*dY22dp(ix,j)+dA22ndt*dY22ndp(ix,j)));
                ap(tt,ix,j)=(yi(7,end)/sinlat*(d2A20dt*dY20dp(ix,j)+d2A22dt*dY22dp(ix,j)+d2A22ndt*dY22ndp(ix,j)));
            else
                up(tt,ix,j)=(yi(3,end)/sinlat*(A20*dY20dp(ix,j)+A22*dY22dp(ix,j)+A22n*dY22ndp(ix,j)));
                vp(tt,ix,j)=(yi(3,end)/sinlat*(dA20dt*dY20dp(ix,j)+dA22dt*dY22dp(ix,j)+dA22ndt*dY22ndp(ix,j)));
                ap(tt,ix,j)=(yi(3,end)/sinlat*(d2A20dt*dY20dp(ix,j)+d2A22dt*dY22dp(ix,j)+d2A22ndt*dY22ndp(ix,j)));
            end
            
            
            
            if bool_rotation_speed
                w_t(tt,ix,j) = vt(tt,ix,j)/ref(1);
                w_p(tt,ix,j) = vp(tt,ix,j)/ref(1);
                
                rotation_speed_unit = 'rad s-1';
            end
            
            
            
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
            dP20dj= (-3)*cos(lat(j))*sin(lat(j));
            dP21dj=3*sin(lat(j))^2 - 3*cos(lat(j))^2;
            dP22dj= 6*cos(lat(j))*sin(lat(j));
            dPotdj(tt,ix,j)= epot*(-3/2*dP20dj*cos_nt+1/4*dP22dj*(3*cos(2*lon(ix))*cos_nt+4*sin(2*lon(ix))*sin_nt));
            deltilt_south(tt,ix,j)= -real( 1 + LN(3) -  LN(1) ) /ref(4) /ref(1) * dPotdj(tt,ix,j);
            
            %derivative wrt lon
            dPotdx(tt,ix,j)= epot*(-3/2*P20*cos_nt+1/4*P22*(-6*sin(2*lon(ix))*cos_nt+8*cos(2*lon(ix))*sin_nt));
            deltilt_east(tt,ix,j)= -real( 1 + LN(3) -  LN(1) ) /ref(4) /ref(1) /sinlat * dPotdx(tt,ix,j);
            
        end
    end

    t=t+Tf/Nt;
    
end

if bool_print_tidal_outputs

    %disp=disp*kpot/Vpot;
    d1=max(max(max(ur)));
    d1c=max(max(max(disp2)));
    disp(['u_r =',num2str(d1),' m']);
    disp(['u_r2 =',num2str(d1c),' m']);
    d2=max(max(max(ut)));
    d3=max(max(max(up)));
    %d12=max(max(max(ur_21)))
    %disp(strcat('u_r =',num2str(d12)))
    %d2=max(max(max(disp2)))
    
    end

%Radial Surface Displacements (m)
if bool_surface_disp_maps
    if bool_sweep_lon
        Cizdir9_New_LonSweep(lon,lat,real(ur)*disp_factor,Nt,time_plot,per,sprintf(['Radial surface displacement (' disp_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);
    elseif bool_sweep_time
        Cizdir9_New_TimeSweep(Time,lat,real(ur)*disp_factor,Nt,lon_plot,per,sprintf(['Radial surface displacement ' disp_unit]),bool_color_num,bool_save,model_name,bool_maps_type);
    end
end
%h=colorbar('v');

%h=figure;
%Cizdir(lon,lat,Pot_St);

if bool_plot_displacement
    h=figure;
    set(gcf, 'renderer', 'zbuffer'); % Prevents colorbar from repeating
    %suptitle(sprintf(['Radial & Lateral Surface Displacements (in ' disp_unit ')']),'FontSize',18)
    suptitle(sprintf(['Radial & Lateral Surface Displacements (in ' disp_unit ')']))
    subplot(2,2,1)
    Cizdir1(lon,lat,real(ur(time_plot,:,:))*disp_factor,Nt,time_plot,per,1);
    %title(sprintf(['u_r (' disp_unit ')']),'FontSize',28)
    title('u_r','FontSize',28)
    colorbar;
    subplot(2,2,2)
    Cizdir1(lon,lat,real(ut(time_plot,:,:))*disp_factor,Nt,time_plot,per,1);
    title('u_\theta','FontSize',28)
    colorbar;
    subplot(2,2,3)
    Cizdir1(lon,lat,real(up(time_plot,:,:))*disp_factor,Nt,time_plot,per,1);
    title('u_\phi','FontSize',28)
    colorbar;
    subplot(2,2,4)
    Cizdir1(lon,lat,Pot(time_plot,:,:),Nt,time_plot,per,1);
    title('\Phi','FontSize',28)
    h2=colorbar('h');
    
    A_SavePlot(bool_save,h,sprintf(['Radial & Lateral Surface Displacements (in ' disp_unit ') ' model_name]));
    
    if bool_plot_displacement_time

        % Equator
        h=figure;
        plot(Time/(60*60),real(ur(:,1,floor(Nlat/2)+1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Radial surface Displacement at Equator (' disp_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Radial surface Displacement at Equator (' disp_unit ') ' model_name]));
        
        h=figure;
        plot(Time/(60*60),real(yi(1,end)*Pot(:,1,floor(Nlat/2)+1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Radial surface Displacement at Equator from Potential (' disp_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Radial surface Displacement at Equator from Potential (' disp_unit ') ' model_name]));
        
    end

    if bool_plot_all_displacement
        
        h=figure;
        plot(Time/(60*60),real(ut(:,1,floor(Nlat/2)+1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Colatitudinal surface Displacement at Equator (' disp_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Colatitudinal surface Displacement at Equator (' disp_unit ') ' model_name]));
        
        h=figure;
        plot(Time/(60*60),real(up(:,1,floor(Nlat/2)+1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Longitudinal surface Displacement at Equator (' disp_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Longitudinal surface Displacement at Equator (' disp_unit ') ' model_name]));
        
    end
    
    if bool_plot_acc
        
        % h=figure;
        % plot(Time/(60*60),real(vr(:,1,floor(Nlat/2)+1)), '-k')
        % xlim([min(Time/(60*60)) max(Time/(60*60))])
        % xlabel('Time (Hours)')
        % ylabel('Radial surface Speed at Equator (m/s)')

        h=figure;
        plot(Time/(60*60),real(ar(:,1,floor(Nlat/2)+1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Radial surface Acceleration at Equator (' disp_unit '.s^-^2)']))
        
        A_SavePlot(bool_save,h,sprintf(['Radial surface Acceleration at Equator (' disp_unit '.s^-^2) ' model_name]));

        h=figure;
        plot(Time/(60*60),real(at(:,1,floor(Nlat/2)+1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Southward surface Acceleration at Equator (' disp_unit '.s^-^2)']));
        
        A_SavePlot(bool_save,h,sprintf(['Southward surface Acceleration at Equator (' disp_unit '.s^-^2) ' model_name]));

        h=figure;
        plot(Time/(60*60),real(ap(:,1,floor(Nlat/2)+1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Eastward surface Acceleration at Equator (' disp_unit '.s^-^2)']));
        
        A_SavePlot(bool_save,h,sprintf(['Eastward surface Acceleration at Equator (' disp_unit '.s^-^2) ' model_name]));
    
    end
    
    if bool_plot_tilt
    
        h=figure;
        plot(Time/(60*60),real(ar(:,1,floor(Nlat/2)+1))/ref(4), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel('Tilt seen by seismometer at Equator (m.s^-^2)')
        
        A_SavePlot(bool_save,h,sprintf(['Tilt seen by seismometer at Equator (m.s^-^2) ' model_name]));

        % % Equator
        % h=figure;
        % plot(Time/(60*60),real(disp2(:,1,floor(Nlat/2)+1)), '-k')
        % xlim([min(Time/(60*60)) max(Time/(60*60))])
        % xlabel('Time (Hours)')
        % ylabel('Surface Displacement (m)')

        h=figure;
        plot(Time/(60*60),real(deltilt_east(:,1,floor(Nlat/2)+1)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel('Eastward tilt at Equator (rad)')
        
        A_SavePlot(bool_save,h,sprintf(['Eastward tilt at Equator (rad) ' model_name]));

        h=figure;
        plot(Time/(60*60),real(deltilt_south(:,1,floor(Nlat/2)+1)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel('Southward tilt at Equator (rad)')
        
        A_SavePlot(bool_save,h,sprintf(['Southward tilt at Equator (rad) ' model_name]));
        
    end
    
    
    if bool_plot_gravity
    
%         h=figure;
%         plot(Time/(60*60),real(deltag(:,1,floor(Nlat/2)+1)*1e5), '-k')
%         xlim([min(Time/(60*60)) max(Time/(60*60))])
%         xlabel('Time (Hours)')
%         ylabel('g at Equator (mgal)')

        h=figure;
        plot(Time/(60*60),real(deltag(:,1,floor(Nlat/2)+1))*gravity_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['\\Delta g at Equator (' gravity_unit ' s^-^2) (local g = ' ...
            num2str(ref(4)*gravity_ref_factor) ' ' gravity_ref_unit '.s^-^2)']));
        
        A_SavePlot(bool_save,h,sprintf(['Delta g at Equator (' gravity_unit ' s-2) ' model_name]));
        
    end
    
    if bool_plot_displacement_time

        % Pole
        h=figure;
        plot(Time/(60*60),real(ur(:,1,1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Radial surface Displacement at North Pole (' disp_unit ') ' model_name]));
        
        A_SavePlot(bool_save,h,sprintf(['Radial surface Displacement at North Pole (' disp_unit ') ' model_name]));
        
        h=figure;
        plot(Time/(60*60),real(yi(1,end)*Pot(:,1,1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Radial surface Displacement at North Pole from Potential (' disp_unit ') ' model_name]));
        
        A_SavePlot(bool_save,h,sprintf(['Radial surface Displacement at North Pole from Potential (' disp_unit ') ' model_name]));

    end
    
    % h=figure;
    % plot(Time/(60*60),real(vr(:,1,1)), '-k')
    % xlim([min(Time/(60*60)) max(Time/(60*60))])
    % xlabel('Time (Hours)')
    % ylabel('Radial surface Speed at North Pole (m/s)')
    
    if bool_plot_all_displacement
        
        h=figure;
        plot(Time/(60*60),real(ut(:,1,1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Longitudinal surface Displacement at North Pole (' disp_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Longitudinal surface Displacement at North Pole (' disp_unit ') ' model_name]));
        
        h=figure;
        plot(Time/(60*60),real(up(:,1,1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Longitudinal surface Displacement at North Pole (' disp_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Longitudinal surface Displacement at North Pole (' disp_unit ') ' model_name]));
        
    end
    
    if bool_plot_acc
        
        h=figure;
        plot(Time/(60*60),real(ar(:,1,1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Radial surface Acceleration at North Pole (' disp_unit '.s^-^2)']));
        
        A_SavePlot(bool_save,h,sprintf(['Radial surface Acceleration at North Pole (' disp_unit '.s^-^2) ' model_name]));

        h=figure;
        plot(Time/(60*60),real(at(:,1,1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Southward surface Acceleration at North Pole (' disp_unit '.s^-^2) ' model_name]));
        
        A_SavePlot(bool_save,h,sprintf(['Southward surface Acceleration at North Pole (' disp_unit '.s^-^2) ' model_name]));

        h=figure;
        plot(Time/(60*60),real(ap(:,1,1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Eastward surface Acceleration at North Pole (' disp_unit '.s^-^2)']))
        
        A_SavePlot(bool_save,h,sprintf(['Eastward surface Acceleration at North Pole (' disp_unit '.s^-^2) ' model_name]));
        
    end
    
    if bool_plot_tilt
    
        h=figure;
        plot(Time/(60*60),real(ar(:,1,1))/ref(4), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel('Tilt seen by seismometer at North Pole (m.s^-^2)')
        
        A_SavePlot(bool_save,h,sprintf(['Tilt seen by seismometer at North Pole (m.s^-^2) ' model_name]));

        h=figure;
        plot(Time/(60*60),real(deltilt_east(:,1,1)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel('Eastward tilt at North Pole (rad)')
        
        A_SavePlot(bool_save,h,sprintf(['Eastward tilt at North Pole (rad) ' model_name]));

        h=figure;
        plot(Time/(60*60),real(deltilt_south(:,1,1)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel('Southward tilt at North Pole (rad)')
        
        A_SavePlot(bool_save,h,sprintf(['Southward tilt at North Pole (rad) ' model_name]));
        
    end
    
    if bool_plot_gravity
    
%         h=figure;
%         plot(Time/(60*60),real(deltag(:,1,1)*1e5), '-k')
%         xlim([min(Time/(60*60)) max(Time/(60*60))])
%         xlabel('Time (Hours)')
%         ylabel('g at North Pole (mgal)')

        h=figure;
        plot(Time/(60*60),real(deltag(:,1,1))*disp_factor, '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['\\Delta g at Pole (' gravity_unit ' s^-^2) (local g = ' ...
            num2str(ref(4)*gravity_ref_factor) ' ' gravity_ref_unit '.s^-^2)']));
        
        A_SavePlot(bool_save,h,sprintf(['Delta g at Pole (' gravity_unit ' s-2) ' model_name]));
        
    end
    
    if bool_plot_rotation_speed
        
        h=figure;
        plot(Time/(60*60),real(w_t(:,1,floor(Nlat/2)+1)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Colatitudinal rotational speed at Equator (' rotation_speed_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Colatitudinal rotational speed at Equator (' rotation_speed_unit ') ' model_name]));
        
        h=figure;
        plot(Time/(60*60),real(w_t(:,1,1)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Colatitudinal rotational speed at North Pole (' rotation_speed_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Colatitudinal rotational speed at North Pole (' rotation_speed_unit ') ' model_name]));
        
        
        h=figure;
        plot(Time/(60*60),real(w_p(:,1,floor(Nlat/2)+1)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Longitudinal rotational speed at Equator (' rotation_speed_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Longitudinal rotational speed at Equator (' rotation_speed_unit ') ' model_name]));
        
        h=figure;
        plot(Time/(60*60),real(w_p(:,1,1)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Longitudinal rotational speed at North Pole (' rotation_speed_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Longitudinal rotational speed at North Pole (' rotation_speed_unit ') ' model_name]));
        
        
        h=figure;
        plot(Time/(60*60),real(sqrt(w_t(:,1,floor(Nlat/2)+1).^2 + w_p(:,1,floor(Nlat/2)+1).^2)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Horizontal rotational speed at Equator (' rotation_speed_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Horizontal rotational speed at Equator (' rotation_speed_unit ') ' model_name]));
        
        h=figure;
        plot(Time/(60*60),real(sqrt(w_t(:,1,1).^2 + w_p(:,1,1).^2)), '-k')
        xlim([min(Time/(60*60)) max(Time/(60*60))])
        xlabel('Time (Hours)')
        ylabel(sprintf(['Horizontal rotational speed at North Pole (' rotation_speed_unit ')']))
        
        A_SavePlot(bool_save,h,sprintf(['Horizontal rotational speed at North Pole (' rotation_speed_unit ') ' model_name]));
        
    end
    
end


% surface_displacemement.radial=real(ur);
% surface_displacemement.theta=real(ut);
% surface_displacemement.phi=real(ut);
% surface_displacemement.lon=lon;
% surface_displacemement.lat=lat;
% surface_displacemement.Nt=Nt;
% surface_displacemement.Time=Time;
% surface_displacemement.per=per;
% delta_g=real(deltag);
% delta_tilt=real(deltilt);
% 
% nur=ur;
% ndeltilt=deltilt;
% ndeltag=deltag;


% save('dispnoocean.mat', 'nur', 'lon', 'lat' , 'Time')
% save('Tiltnoocean.mat', 'ndeltilt', 'lat' , 'Time')
% save('gravnoocean.mat', 'ndeltag', 'lon', 'lat' , 'Time')


%save('Europadispoft.mat', 'ur', 'lon', 'lat' , 'Time')
%save('EuropaTiltsoft.mat', 'deltilt', 'lat' , 'Time')
%save('Europagravsoft.mat', 'deltag', 'lon', 'lat' , 'Time')

%save('dispocean.mat', 'ur', 'lon', 'lat' , 'Time')
%save('Tiltocean.mat', 'deltilt', 'lat' , 'Time')
%save('gravocean.mat', 'deltag', 'lon', 'lat' , 'Time')