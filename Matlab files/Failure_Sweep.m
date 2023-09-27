%close all;
%press_s=[press(1) reshape([press;press],1,nbcouches*2)];
%press_s=[press(1) reshape([press;press],1,nbcouches*2)];
g_s =gravity(rho_s,r_s,2*nbcouches+1,G);% calcul de g pour chaque couche
press_s = Pressure(r_s,rho_s,g_s);


%% Sweep on cohesion and friction
tic

cohe = sort(cohe,'descend'); % For plotting purposes

n_cohe = numel(cohe);
n_friction = numel(friction);

failure_radius_rel = cell(n_cohe,n_friction);
failure_time_rel = cell(n_cohe,n_friction);
failure_lon_rel = cell(n_cohe,n_friction);
failure_lat_rel = cell(n_cohe,n_friction);

failure_radius_tot = zeros(n_cohe,n_friction);
failure_radius_tot_bis = zeros(n_cohe,n_friction);

% Weights to convert index to actual depths
depth_interval = r_s;
depth_interval(2:end) = r_s(2:end) - r_s(1:end-1);
depth_interval = depth_interval/sum(depth_interval);


hbar = waitbar(0, 'Cohesion and friction sweep...');

for i_c = 1:numel(cohe)

    waitbar(i_c/n_cohe, hbar, sprintf('Cohesion loop, at step %n', i_c));

    for i_f = 1:numel(friction)

        cohe_current = cohe(i_c);
        friction_current = friction(i_f);

        % Reinitialization to allow stacking
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
                            %c_m(tt,ix,j,radius) = abs(s1(tt,ix,j,radius)+s3(tt,ix,j,radius))/2+press_s(radius); % Abscissa of center of Mohr circle
                            c_m(tt,ix,j,radius) = (s1(tt,ix,j,radius)+s3(tt,ix,j,radius))/2+press_s(radius); % Abscissa of center of Mohr circle
                        else
                            %c_m(tt,ix,j,radius) = abs(s1(tt,ix,j,radius)+s3(tt,ix,j,radius))/2;
                            c_m(tt,ix,j,radius) = (s1(tt,ix,j,radius)+s3(tt,ix,j,radius))/2;
                        end
                        % Mohr criterion : y = c + mu*x, with x = normal stress sigma_m and y = shear stress tau_m
                        d_y = 1;
                        d_x = sin(friction_current);%tan(friction_current);
                        d_c = cohe_current*cos(friction_current);%cohe_current;

                        C_m_c(tt,ix,j,radius) = abs(d_x*c_m(tt,ix,j,radius)+cohe_current)/sqrt(d_y^2+d_x^2); % distance between criterion line and center
                        criterion(tt,ix,j,radius) = C_m_c(tt,ix,j,radius);

                        % Angle between criterion line and x axis
                        slope_angle = friction_current;
                        criterion_angle(tt,ix,j,radius) = (pi/2 - friction_current); % 2theta in classical notation

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
        failure_radius_rel{i_c,i_f} = 100*failure_radius/((Ntimeloop+1)*Nlon*Nlat); % /((Ntimeloop+1)*Nradius*Nlon*Nlat);
        failure_time_rel{i_c,i_f} = 100*failure_time/(Nradius*Nlon*Nlat); %/((Ntimeloop+1)*Nradius*Nlon*Nlat);
        failure_lon_rel{i_c,i_f} = 100*failure_lon/((Ntimeloop+1)*Nradius*Nlat); %/((Ntimeloop+1)*Nradius*Nlon*Nlat);
        failure_lat_rel{i_c,i_f} = 100*failure_lat/((Ntimeloop+1)*Nradius*Nlon); %/((Ntimeloop+1)*Nradius*Nlon*Nlat);

        % Summation
        %failure_radius_tot(i_c,i_f) = sum(failure_radius_rel{i_c,i_f})/numel(r_s); % Check that both are the same size
        failure_radius_tot_bis(i_c,i_f) = sum(failure_radius_rel{i_c,i_f}.*depth_interval);

        % Failure_Depth;
        %
        % % % Relative value for individual failure plots at sc (radius)
        % % failure_radius_rel_top = failure_radius_rel(top_layers_idx);
        % % failure_radius_rel_bot = failure_radius_rel(bot_layers_idx);
        %
%% Individual plots

        % if bool_plot_failure
        %     h=figure;
        %     plot(squeeze(criterion(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor,'r');
        %     %plot(InterpCriterion(squeeze(criterion(time_plot,lon_plot,lat_plot,:)),r_s/1000,0),r_s/1000,'r');
        %     hold on;
        %     plot(squeeze(tau_m(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor,'b');
        %     hold off;
        %     xlabel(sprintf(['C_m_c and \\\\tau_m (' stress_unit ')']));
        %     ylabel(sprintf(['Radius (' depth_unit ')']));
        %     legend(sprintf(['C_m_c for Failure Criterion (' stress_unit ')']),sprintf(['Maximal shear stress \\tau_m (' stress_unit ')']));
        %
        %     A_SavePlot(bool_save,h,sprintf(['Radius tau_m and C_m_c ' model_name]));
        %
        %     h=figure;
        %     plot(squeeze(criterion(time_plot,lon_plot,:,radius_plot))*stress_factor,lat*180/pi,'r');
        %     hold on;
        %     plot(squeeze(tau_m(time_plot,lon_plot,:,radius_plot))*stress_factor,lat*180/pi,'b');
        %     hold off;
        %     xlabel(sprintf(['C_m_c and \\tau_m (' stress_unit ')']));
        %     ylabel('Colatitude (°)');
        %     legend(sprintf(['C_m_c for Failure Criterion (' stress_unit ')']),sprintf(['Maximal shear stress \\tau_m (' stress_unit ')']));
        %
        %     A_SavePlot(bool_save,h,sprintf(['Colat tau_m and C_m_c ' model_name]));
        %
        %     h=figure;
        %     plot(squeeze(criterion(time_plot,:,lat_plot,radius_plot))*stress_factor,lon*180/pi,'r');
        %     hold on;
        %     plot(squeeze(tau_m(time_plot,:,lat_plot,radius_plot))*stress_factor,lon*180/pi,'b');
        %     hold off;
        %     xlabel(sprintf(['C_m_c and \\tau_m (' stress_unit ')']));
        %     ylabel('Longitude (°)');
        %     legend(sprintf(['C_m_c for Failure Criterion (' stress_unit ')']),sprintf(['Maximal shear stress \\tau_m (' stress_unit ')']));
        %
        %     A_SavePlot(bool_save,h,sprintf(['Long tau_m and C_m_c ' model_name]));
        %
        %     h=figure;
        %     plot(squeeze(criterion(:,lon_plot,lat_plot,radius_plot))*stress_factor,Time/3600,'r');
        %     hold on;
        %     plot(squeeze(tau_m(:,lon_plot,lat_plot,radius_plot))*stress_factor,Time/3600,'b');
        %     hold off;
        %     xlabel(sprintf(['C_m_c and \\tau_m (' stress_unit ')']));
        %     ylabel('Time (hours)');
        %     legend(sprintf(['C_m_c for Failure Criterion (' stress_unit ')']),sprintf(['Maximal shear stress \\tau_m (' stress_unit ')']));
        %
        %     A_SavePlot(bool_save,h,sprintf(['Time tau_m and C_m_c ' model_name]));
        %
        % end

    end
end

close(hbar);

Failure_Sweep_Plots;

%% Save files
if bool_data_save

    save(sprintf([nom_file '.mat']),'failure_radius_rel','failure_radius_tot','r_s','cohe','friction');
    A_SavePlot(1,h,sprintf([nom_file '_failure_map']));

    save(sprintf([nom_file '_CriterionStress.mat']),'Ntimeloop','Nlon','Nlat','Nradius',...
        'Time','lon','lat','criterion','tau_m');
end

toc
