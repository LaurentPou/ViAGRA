%% Heatmap plots
size_font = 14;%22;
set(0,'DefaultTextFontSize',size_font);
set(0,'DefaultAxesFontSize',size_font);
%bool_save = 0; %
%bool_data_save = 0; %

y_axis = cohe; % First argument in matrix is column
x_axis = friction;

% Global colormap
h=figure;
heatmap(x_axis*180/pi,y_axis,failure_radius_tot_bis,'colormap',parula);
%heatmap(tan(x_axis),y_axis,failure_radius_tot,'colormap',parula);
xlabel('Angle of internal friction (°)');
ylabel('Cohesion (Pa)');
title(sprintf(['Global relative rupture ' legend_pressure]))