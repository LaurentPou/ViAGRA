%% Plot time plot tests


% % Missing plots parameters
% nom_file='Moon_Weber_2011';
% stress_factor = 1;
% depth_factor = 1000;
% stress_unit = 'Pa';
% depth_unit = 'km';
% 
% % For running Failure_Sweep
% %bool_pressure_failure = 1;
% %legend_pressure = 'with pressure in failure criterion';
% %bool_data_save = 1;

%
time_span = 1:1:Ntimeloop;
lonlat_plot = {[1 1];[1 2];[1 3];[1 4];[1 5];[2 2];[2 3];[2 4];[3 2];[3 3];[3 4];...
    [7 2];[7 3];[7 4];[8 2];[8 3];[8 4]};
radius_plot = 1;

xmin = min(min(tau_m,[],'all'),min(sigma_m,[],'all'));
xmax = max(max(tau_m,[],'all'),max(sigma_m,[],'all'));


for ii = 1:numel(lonlat_plot)

    lon_i = lonlat_plot{ii}(1);
    lat_i = lonlat_plot{ii}(2);

    % Save as video
    aviobj = VideoWriter(sprintf([nom_file '_Stress_lon_' num2str(lon(lon_i)*180/pi) '_lat_' num2str(lat(lat_i)*180/pi) '.avi']));
    open(aviobj);

    f1 = figure;
    set(f1,'Units','Normalized','OuterPosition',[0 0.5 0.5 0.5]); % pos x (px) pos y (px) width x (%) width y (%)
    set(f1,'Units','Inches');

    for tt = time_span

        plot(squeeze(tau_m(tt,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
        hold on;
        plot(squeeze(sigma_m(tt,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
        plot(squeeze(criterion(tt,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
        plot(squeeze(-criterion(tt,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
        hold off;
        xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
        ylabel(sprintf(['Radius (' depth_unit ')']));
        %legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
        %    sprintf(['Failure Criterion C_{mc} (' stress_unit ')']),sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
        title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
        xlim([xmin xmax])

        % Freeze frame for longer
        for j = 1:1
            F = getframe(f1);
            writeVideo(aviobj,F)
        end
    end

    close(f1)
    close(aviobj);

end

% %% Plot orbital movement
% close all;
% 
% bool_exageration = 10;
% Ntimeloop = 400;
% 
% time_axis = 0:1:Ntimeloop;
% M = (2*pi/(Ntimeloop))*time_axis; % in rad
% 
% time_axis_plot = 0:0.1:Ntimeloop;
% M_plot = (2*pi/(Ntimeloop))*time_axis_plot; % in rad
% 
% e=0.0549;%Moon orbit eccentricity
% a=3.84399e8;
% 
% %e=0.0093;%Europe orbit eccentricity
% %a=670.9e6;
% 
% if bool_exageration
%     e = 0.3;%e*bool_exageration;
% end
% 
% % Basic eccentric anomaly solver: M = E - e sin E
% Etest = M + e*sin(M);
% Etest_c = M + e*sin(Etest);
% tolerance = 1e-5;
% 
% tic
% while(sum(Etest) - sum(Etest_c) >= tolerance)
%     Etest = Etest_c;
%     Etest_c = M + e*sin(Etest);
% end
% E = Etest_c;
% toc
% 
% % Plot orbits
% M1_e = [a*e 0];
% 
% x_c = a*cos(M_plot) + M1_e(1);
% y_c = a*sin(M_plot) + M1_e(2);
% 
% b = a*sqrt(1-e^2);
% x_e = a*cos(M_plot);
% y_e = b*sin(M_plot);
% 
% tic
% f_orbit = figure;
% hold on;
% for i_plot = 1:numel(M_plot)
%     plot(x_c(i_plot),y_c(i_plot),'ko','MarkerSize',0.5);
%     plot(x_e(i_plot),y_e(i_plot),'bo','MarkerSize',0.5);
% end
% plot(M1_e(1),M1_e(2),'go','MarkerSize',10)
% xline(M1_e(1),'k:');
% xline(0,'b');
% yline(0,'b');
% 
% toc
% ax_orbit = gca;
% 
% tic
% 
% %% Plot movement
% % Save as video
% aviobj_motion = VideoWriter('Orbital motion');
% open(aviobj_motion);
% 
% for j_plot = 1:numel(M)
% %     f_copy = figure;
% %     ax_copy = copyobj(ax_orbit,f_copy);
% %     hold on;
%     plot(a*cos(M(j_plot)) + M1_e(1),a*sin(M(j_plot)) + M1_e(2),'ko','MarkerSize',6)
%     plot(a*cos(E(j_plot)),b*sin(E(j_plot)),'bo','MarkerSize',6)
% 
%     % Freeze frame for longer
%     for j = 1:1
%         %F_motion = getframe(f_copy);
%         F_motion = getframe(f_orbit);
%         writeVideo(aviobj_motion,F_motion)
%     end
% 
%     %close(f_copy);
% 
% end
% 
% % Change framerate
% %aviobj_motion.FrameRate = 1;
% 
% close(aviobj_motion);
% close(f_orbit);
% 
% 
% toc
% 
% 
% % From the point of view of the satellite

