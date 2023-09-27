%close all;
%press_s=[press(1) reshape([press;press],1,nbcouches*2)];
%press_s=[press(1) reshape([press;press],1,nbcouches*2)];
g_s =gravity(rho_s,r_s,2*nbcouches+1,G);% calcul de g pour chaque couche
press_s = Pressure(r_s,rho_s,g_s);

failure_reached = zeros(Ntimeloop+1,Nlon,Nlat,Nradius);

failure_radius = zeros(1,Nradius);
failure_radius_idx = zeros(1,Nradius);

failure_time = zeros(1,Ntimeloop+1);
failure_time_idx = zeros(1,Ntimeloop+1);

failure_lon = zeros(1,Nlon);
failure_lon_idx = zeros(1,Nlon);

failure_lat = zeros(1,Nlat);
failure_lat_idx = zeros(1,Nlat);

for tt=1:1:Ntimeloop+1
    for ix=1:1:Nlon
        for j=1:Nlat
            for radius=1:1:Nradius

                tau_m(tt,ix,j,radius) = abs(s1(tt,ix,j,radius)-s3(tt,ix,j,radius))/2; % Mohr radius
                
                if bool_pressure_failure
                    c_m(tt,ix,j,radius) = abs(s1(tt,ix,j,radius)+s3(tt,ix,j,radius))/2+press_s(radius); % Abscissa of center of Mohr circle
                else
                    c_m(tt,ix,j,radius) = abs(s1(tt,ix,j,radius)+s3(tt,ix,j,radius))/2;
                end
                % Mohr criterion : y = c + mu*x, with x = normal stress sigma_m and y = shear stress tau_m
                d_y = 1;
                d_x = sin(friction);%tan(friction);
                d_c = cohe*cos(friction);%cohe;
                
                C_m_c(tt,ix,j,radius) = abs(d_x*c_m(tt,ix,j,radius)+cohe)/sqrt(d_y^2+d_x^2); % distance between criterion line and center
                criterion(tt,ix,j,radius) = C_m_c(tt,ix,j,radius);

                % Angle between criterion line and x axis
                slope_angle = friction;
                criterion_angle(tt,ix,j,radius) = (pi/2 - friction); % 2theta in classical notation
                 
                % Comparison of distance between criterion line and center
                if criterion(tt,ix,j,radius) - tau_m(tt,ix,j,radius) <= 0

                    failure_reached(tt,ix,j,radius) = 1;
                    
                    % Individual dim failure plots
                    failure_radius(radius) = failure_radius(radius) + 1;
                    failure_radius_idx(radius) = radius;

                    failure_time(tt) = failure_time(tt) + 1;
                    failure_time_idx(tt) = tt;

                    failure_lon(ix) = failure_lon(ix) + 1;
                    failure_lon_idx(ix) = ix;

                    failure_lat(j) = failure_lat(j) + 1;
                    failure_lat_idx(j) = j;

%                     % If failure is reached, there is possibly several angles possible for fault orientation
%                     % Equation: x^2 (1+tan(friction)^2) + 2x (c*tan(friction) - c_m) + c^2 + c_m^2 = tau_m^2
%                     tau_m0 = (tau_m(tt,ix,j,radius));
%                     c_m0 = c_m(tt,ix,j,radius);
%                     
%                     a_2 = 1 + (tan(friction))^2;
%                     b_2 = 2*(cohe*tan(friction) - c_m0);
%                     c_2 = cohe^2 + (c_m0)^2 - (tau_m0)^2;
% 
%                     x_sol = roots([a_2 b_2 c_2]);
%                     x_sol1 = min(x_sol);
% 
% %                     if x_sol1 <= 0
% %                         MohrCoulombVisu;
% %                         pause(2);
% %                     end
% 
%                     y_sol1 = cohe + x_sol1*tan(friction);
%                     %y_sol2 = cohe + x_sol2*tan(friction);
%                     
%                     alpha_angle = asin(y_sol1/tau_m0);
%                     theta_sol1(tt,ix,j,radius) = alpha_angle;
% 
%                     %if x_sol2 < c_m(tt,ix,j,radius)
%                         %theta_sol2bis(tt,ix,j,radius) = asin(y_sol2/tau_m(tt,ix,j,radius));
%                     %else
%                         %theta_sol2bis(tt,ix,j,radius) = pi - asin(y_sol2/tau_m(tt,ix,j,radius));
%                     %end
%                     
%                     delta_angle = (pi/2 - friction) - alpha_angle; % Value from sol1
%                     theta_sol2(tt,ix,j,radius) = (pi/2 - friction) + delta_angle;

                else
% 
%                     theta_sol1(tt,ix,j,radius) = nan;
%                     theta_sol2(tt,ix,j,radius) = nan;
%                     %theta_sol2bis(tt,ix,j,radius) = 0;

                end
            end
        end
    end
end

% Relative value for individual failure plots
failure_radius_rel = 100*failure_radius/((Ntimeloop+1)*Nlon*Nlat); % /((Ntimeloop+1)*Nradius*Nlon*Nlat);
failure_time_rel = 100*failure_time/(Nradius*Nlon*Nlat); %/((Ntimeloop+1)*Nradius*Nlon*Nlat);
failure_lon_rel = 100*failure_lon/((Ntimeloop+1)*Nradius*Nlat); %/((Ntimeloop+1)*Nradius*Nlon*Nlat);
failure_lat_rel = 100*failure_lat/((Ntimeloop+1)*Nradius*Nlon); %/((Ntimeloop+1)*Nradius*Nlon*Nlat);

Failure_Depth;

% % Relative value for individual failure plots at sc (radius)
% failure_radius_rel_top = failure_radius_rel(top_layers_idx);
% failure_radius_rel_bot = failure_radius_rel(bot_layers_idx);

%% Coordinates plots

if bool_plot_failure
    h=figure;
    plot(squeeze(criterion(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor,'r');
    %plot(InterpCriterion(squeeze(criterion(time_plot,lon_plot,lat_plot,:)),r_s/1000,0),r_s/1000,'r');
    hold on;
    plot(squeeze(tau_m(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor,'b');
    hold off;
    xlabel(sprintf(['C_m_c and \\tau_m (' stress_unit ')']));
    ylabel(sprintf(['Radius (' depth_unit ')']));
    legend(sprintf(['C_m_c for Failure Criterion (' stress_unit ')']),sprintf(['Maximal shear stress \\tau_m (' stress_unit ')']));
    
    A_SavePlot(bool_save,h,sprintf(['Radius tau_m and C_m_c ' model_name]));

    h=figure;
    plot(squeeze(criterion(time_plot,lon_plot,:,radius_plot))*stress_factor,lat*180/pi,'r');
    hold on;
    plot(squeeze(tau_m(time_plot,lon_plot,:,radius_plot))*stress_factor,lat*180/pi,'b');
    hold off;
    xlabel(sprintf(['C_m_c and \\tau_m (' stress_unit ')']));
    ylabel('Colatitude (°)');
    legend(sprintf(['C_m_c for Failure Criterion (' stress_unit ')']),sprintf(['Maximal shear stress \\tau_m (' stress_unit ')']));

    A_SavePlot(bool_save,h,sprintf(['Colat tau_m and C_m_c ' model_name]));
    
    h=figure;
    plot(squeeze(criterion(time_plot,:,lat_plot,radius_plot))*stress_factor,lon*180/pi,'r');
    hold on;
    plot(squeeze(tau_m(time_plot,:,lat_plot,radius_plot))*stress_factor,lon*180/pi,'b');
    hold off;
    xlabel(sprintf(['C_m_c and \\tau_m (' stress_unit ')']));
    ylabel('Longitude (°)');
    legend(sprintf(['C_m_c for Failure Criterion (' stress_unit ')']),sprintf(['Maximal shear stress \\tau_m (' stress_unit ')']));
    
    A_SavePlot(bool_save,h,sprintf(['Long tau_m and C_m_c ' model_name]));

    h=figure;
    plot(squeeze(criterion(:,lon_plot,lat_plot,radius_plot))*stress_factor,Time/3600,'r');
    hold on;
    plot(squeeze(tau_m(:,lon_plot,lat_plot,radius_plot))*stress_factor,Time/3600,'b');
    hold off;
    xlabel(sprintf(['C_m_c and \\tau_m (' stress_unit ')']));
    ylabel('Time (hours)');
    legend(sprintf(['C_m_c for Failure Criterion (' stress_unit ')']),sprintf(['Maximal shear stress \\tau_m (' stress_unit ')']));
    
    A_SavePlot(bool_save,h,sprintf(['Time tau_m and C_m_c ' model_name]));
    
end

%% Sweep plots

if bool_surface_Mohr_maps
    if bool_sweep_lon
        
        Cizdir9_New_LonSweep(lon,lat,tau_m(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['Max shear stress \\\\tau_m (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);
        Cizdir9_New_LonSweep(lon,lat,criterion(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['Distance C_m_c between Criterion and Mohr circle center \\n(' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        if bool_color_num
            Cizdir9_New_LonSweep(lon,lat,(- criterion(:,:,:,radius_plot) + tau_m(:,:,:,radius_plot))*stress_factor,Nt,time_plot,per,sprintf(['Difference between Max shear stress \\\\tau_m and Failure \\ncriterion C_m_c (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);
        else
            Cizdir9_New_LonSweep(lon,lat,(- criterion(:,:,:,radius_plot) + tau_m(:,:,:,radius_plot))*stress_factor,Nt,time_plot,per,sprintf(['Areas where Failure criterion is reached (' stress_unit ')\\n']),bool_color_num,bool_save,model_name,bool_maps_type);
        end
    
    elseif bool_sweep_time
        
        Cizdir9_New_TimeSweep(Time,lat,tau_m(:,:,:,radius_plot)*stress_factor,Nt,lon_plot,per,sprintf(['Max shear stress \\\\tau_m (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);
        Cizdir9_New_TimeSweep(Time,lat,criterion(:,:,:,radius_plot)*stress_factor,Nt,lon_plot,per,sprintf(['Distance C_m_c between Criterion and Mohr circle center \\n(' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);
        
        if bool_color_num
            Cizdir9_New_TimeSweep(Time,lat,(- criterion(:,:,:,radius_plot) + tau_m(:,:,:,radius_plot))*stress_factor,Nt,lon_plot,per,sprintf(['Difference between Max shear stress \\\\tau_m and Failure \\ncriterion C_m_c (' stress_unit ')\\n']),bool_color_num,bool_save,model_name,bool_maps_type);
        else
            Cizdir9_New_TimeSweep(Time,lat,(- criterion(:,:,:,radius_plot) + tau_m(:,:,:,radius_plot))*stress_factor,Nt,lon_plot,per,sprintf(['Areas where Failure criterion is reached (' stress_unit ')\\n'])',bool_color_num,bool_save,model_name,bool_maps_type);
        end
        
    end
end