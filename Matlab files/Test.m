close all;

%% Colormap tests
% h=figure;
% imagesc(x_axis,y_axis,failure_reached_array);
% xlabel('Cohesion (Pa)');
% ylabel('Angle of friction (°)');
% xlim([min(x_axis) max(x_axis)]);
% ylim([min(y_axis) max(y_axis)]);
%
% figure;
% imagesc(failure_reached_array);
% axes = get(gcf,'CurrentAxes');
% xlim([1 numel(x_axis)]);
% ylim([1 numel(y_axis)]);

% axes.XTick = x_axis;
% axes.XTickLabel = {num2str(x_axis(1)),num2str(x_axis(2))};
%
% yy_axis = sort(y_axis);
% axes.YTick = yy_axis;
% axes.YTickLabel = {num2str(yy_axis(1)),num2str(yy_axis(2))};

% set(gca,'XTick',x_axis);
% set(gca,'XTickLabel',{num2str(x_axis(1)),num2str(x_axis(2))});

% A = failure_reached_array;
%
% figure
% imshow(A)
% figure
% imagesc(A)
% figure
% pcolor(A) ;
% shading interp;
% colorbar

% figure;
% heatmap(x_axis,y_axis,A,'colormap',parula)

% %% Plot tests
% close all;
% 
% time_span = 1:Ntimeloop;
% lonlat_plot = {[1 1];[1 4];[1 7];[1 10];[1 13];[5 4];[5 7];[5 10];[9 4];[9 7];[9 10];...
%     [25 4];[25 7];[25 10];[29 4];[29 7];[29 10]};
% radius_plot = 1;
% 
% xmin = min(min(tau_m,[],'all'),min(sigma_m,[],'all'));
% xmax = max(max(tau_m,[],'all'),max(sigma_m,[],'all'));
% 
% 
% for ii = 1:numel(lonlat_plot)
% 
%     lon_i = lonlat_plot{ii}(1);
%     lat_i = lonlat_plot{ii}(2);
% 
%     % Save as video
%     aviobj = VideoWriter(sprintf([nom_file '_Stress_lon_' num2str(lon(lon_i)*180/pi) '_lat_' num2str(lat(lat_i)*180/pi) '.avi']));
%     open(aviobj);
% 
%     f1 = figure;
%     set(f1,'Units','Normalized','OuterPosition',[0 0.5 0.5 0.5]); % pos x (px) pos y (px) width x (%) width y (%)
%     set(f1,'Units','Inches');
% 
%     for tt = time_span
% 
%         plot(squeeze(tau_m(tt,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
%         hold on;
%         plot(squeeze(sigma_m(tt,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
%         plot(squeeze(criterion(tt,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
%         plot(squeeze(-criterion(tt,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
%         hold off;
%         xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
%         ylabel(sprintf(['Radius (' depth_unit ')']));
%         legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
%             sprintf(['Failure Criterion C_{mc} (' stress_unit ')']),sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
%         title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
%         xlim([xmin xmax])
% 
%         % Freeze frame for longer
%         for j = 1:25
%             F = getframe(f1);
%             writeVideo(aviobj,F)
%         end
%     end
% 
%     close(f1)
%     close(aviobj);
% 
% end

%% Failure map correction
depth_interval = r_s;
depth_interval(2:end) = r_s(2:end) - r_s(1:end-1);

% Normalization to 1
depth_interval = depth_interval/sum(depth_interval);

for i_c = 1:numel(cohe)
    for i_f = 1:numel(friction)

        failure_radius_tot_bis(i_c,i_f) = sum(failure_radius_rel{i_c,i_f}.*depth_interval);

    end
end
