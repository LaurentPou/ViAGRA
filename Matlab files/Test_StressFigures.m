
% %% Max shear stress for A1
% close all;
% 
% xmin = 0;%min(min(tau_m,[],'all'),min(sigma_m,[],'all'));
% xmax = max(max(tau_m,[],'all'),max(sigma_m,[],'all'));
% 
% lon_i = lon_plot;
% lat_i = lat_plot;
% 
% 
% figure;
% plot(squeeze(tau_m(2,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
% hold on;
% plot(squeeze(criterion(2,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% hold off;
% xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
% ylabel(sprintf(['Radius (' depth_unit ')']));
% legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),...
%     sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
% %title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
% xlim([xmin xmax]);
% 
% %% Max shear stress for M = pi/2
% close all;
% 
% xmin = 0;%min(min(tau_m,[],'all'),min(sigma_m,[],'all'));
% xmax = max(max(tau_m,[],'all'),max(sigma_m,[],'all'));
% 
% lon_i = 1;%lon_plot;
% lat_i = 1;%lat_plot;
% 
% 
% figure;
% plot(squeeze(tau_m(3,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
% hold on;
% plot(squeeze(criterion(3,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% hold off;
% xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
% ylabel(sprintf(['Radius (' depth_unit ')']));
% legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),...
%     sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
% %title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
% xlim([xmin xmax]);


%% Stress radius plots
close all;

xmin = min(min(tau_m,[],'all'),min(sigma_m,[],'all'));
xmax = max(max(tau_m,[],'all'),max(sigma_m,[],'all'));

lon_i = 1;%lon_plot;
lat_i = 1;%lat_plot;

figure;
plot(squeeze(tau_m(1,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
hold on;
plot(squeeze(sigma_m(1,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
plot(squeeze(criterion(1,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
plot(squeeze(-criterion(1,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
hold off;
xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
ylabel(sprintf(['Radius (' depth_unit ')']));
legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
    sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
%title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
xlim([xmin xmax]);

figure;
plot(squeeze(tau_m(49,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
hold on;
plot(squeeze(sigma_m(49,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
plot(squeeze(criterion(49,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
plot(squeeze(-criterion(49,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
hold off;
xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
ylabel(sprintf(['Radius (' depth_unit ')']));
legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
    sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
%title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
xlim([xmin xmax]);

figure;
plot(squeeze(tau_m(97,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
hold on;
plot(squeeze(sigma_m(97,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
plot(squeeze(criterion(97,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
plot(squeeze(-criterion(97,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
hold off;
xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
ylabel(sprintf(['Radius (' depth_unit ')']));
legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
    sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
%title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
xlim([xmin xmax]);

figure;
plot(squeeze(tau_m(145,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
hold on;
plot(squeeze(sigma_m(145,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
plot(squeeze(criterion(145,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
plot(squeeze(-criterion(145,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
hold off;
xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
ylabel(sprintf(['Radius (' depth_unit ')']));
legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
    sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
%title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
xlim([xmin xmax]);

figure;
plot(squeeze(tau_m(193,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
hold on;
plot(squeeze(sigma_m(193,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
plot(squeeze(criterion(193,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
plot(squeeze(-criterion(193,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
hold off;
xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
ylabel(sprintf(['Radius (' depth_unit ')']));
legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
    sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
%title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
xlim([xmin xmax]);

%% Stress radius plots
% close all;
% 
% xmin = min(min(tau_m,[],'all'),min(sigma_m,[],'all'));
% xmax = max(max(tau_m,[],'all'),max(sigma_m,[],'all'));
% 
% lon_i = 1;%lon_plot;
% lat_i = 1;%lat_plot;
% 
% figure;
% plot(squeeze(tau_m(1,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
% hold on;
% plot(squeeze(sigma_m(1,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
% plot(squeeze(criterion(1,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% plot(squeeze(-criterion(1,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% hold off;
% xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
% ylabel(sprintf(['Radius (' depth_unit ')']));
% legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
%     sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
% %title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
% xlim([xmin xmax]);
% 
% figure;
% plot(squeeze(tau_m(2,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
% hold on;
% plot(squeeze(sigma_m(2,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
% plot(squeeze(criterion(2,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% plot(squeeze(-criterion(2,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% hold off;
% xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
% ylabel(sprintf(['Radius (' depth_unit ')']));
% legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
%     sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
% %title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
% xlim([xmin xmax]);
% 
% figure;
% plot(squeeze(tau_m(3,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
% hold on;
% plot(squeeze(sigma_m(3,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
% plot(squeeze(criterion(3,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% plot(squeeze(-criterion(3,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% hold off;
% xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
% ylabel(sprintf(['Radius (' depth_unit ')']));
% legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
%     sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
% %title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
% xlim([xmin xmax]);
% 
% figure;
% plot(squeeze(tau_m(4,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
% hold on;
% plot(squeeze(sigma_m(4,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
% plot(squeeze(criterion(4,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% plot(squeeze(-criterion(4,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% hold off;
% xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
% ylabel(sprintf(['Radius (' depth_unit ')']));
% legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
%     sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
% %title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
% xlim([xmin xmax]);
% 
% figure;
% plot(squeeze(tau_m(5,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'b','LineWidth',3);
% hold on;
% plot(squeeze(sigma_m(5,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'r','LineWidth',3);
% plot(squeeze(criterion(5,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% plot(squeeze(-criterion(5,lon_i,lat_i,:))*stress_factor,r_s*depth_factor,'k--','LineWidth',3);
% hold off;
% xlabel(sprintf(['\\tau_m and \\sigma_m(' stress_unit ')']));
% ylabel(sprintf(['Radius (' depth_unit ')']));
% legend(sprintf(['Shear stress \\tau_m (' stress_unit ')']),sprintf(['Normal stress \\sigma_m (' stress_unit ')']),...
%     sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));%,sprintf(['Failure Criterion C_{mc} (' stress_unit ')']));
% %title(sprintf(['Stress lon ' num2str(lon(lon_i)*180/pi) '°, colat ' num2str(lat(lat_i)*180/pi) '°, time ' num2str(tt) ' out of ' num2str(Ntimeloop)]));
% xlim([xmin xmax]);
