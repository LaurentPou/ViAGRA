close all;

time_plotA = 1;
lon_plotA = lon_plot;%round((360-31.1)/Dlon);
lat_plotA = lat_plot;%round((90-13.2)/Dlat);
nb_rs = length(r_s);
array_rs = 1:round(nb_rs/4):nb_rs;

for time_plotA = [52]%[4 12 44 52 53]%1:Nt %[2 21 22]

    h=figure;
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rr(time_plotA,lon_plotA,lat_plotA,:)))*stress_factor,'-^','Color','r','MarkerIndices',array_rs,'LineWidth',3);
    hold on;
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_tt(time_plotA,lon_plotA,lat_plotA,:)))*stress_factor,'-s','Color','g','MarkerIndices',array_rs,'LineWidth',3);
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_pp(time_plotA,lon_plotA,lat_plotA,:)))*stress_factor,'-o','Color','b','MarkerIndices',array_rs,'LineWidth',3);
    %plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:)))*stress_factor,'magenta');
    plot((ref(1)-r_s).*depth_factor,zeros(1,numel(r_s)),'black');
    hold off;
    xlabel(sprintf(['Depth (' depth_unit ')']));
    if strcmp(stress_unit,stress_r_unit)
        ylabel(sprintf([' Stress (' stress_unit ')']));
    else
        ylabel(sprintf([' Stress (' stress_r_unit ' or ' stress_unit ')']));
    end
    legend(sprintf(['\\sigma_{rr} in ' stress_r_unit]),...
        sprintf(['\\sigma_{\\theta\\theta} in ' stress_r_unit]),...
        sprintf(['\\sigma_{\\phi\\phi} in ' stress_r_unit]));
    title(sprintf(['Timestep ' num2str(time_plotA) ' out of ' num2str(Nt)]))

    h=figure;
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rt(time_plotA,lon_plotA,lat_plotA,:)))*stress_factor,'-^','Color','r','MarkerIndices',array_rs,'LineWidth',3);
    hold on;
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_tp(time_plotA,lon_plotA,lat_plotA,:)))*stress_factor,'-o','Color','g','MarkerIndices',array_rs,'LineWidth',3);
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_pr(time_plotA,lon_plotA,lat_plotA,:)))*stress_factor,'-s','Color','b','MarkerIndices',array_rs,'LineWidth',3);
    plot((ref(1)-r_s).*depth_factor,zeros(1,numel(r_s)),'black');
    hold off;
    xlabel(sprintf(['Depth (' depth_unit ')']));
    if strcmp(stress_unit,stress_r_unit)
        ylabel(sprintf([' Stress (' stress_unit ')']));
    else
        ylabel(sprintf([' Stress (' stress_r_unit ' or ' stress_unit ')']));
    end
    legend(sprintf(['\\sigma_{r\\theta} in ' stress_r_unit]),...
        sprintf(['\\sigma_{\\theta\\phi} in ' stress_r_unit]),...
        sprintf(['\\sigma_{\\phi r} in ' stress_r_unit]));
    title(sprintf(['Timestep ' num2str(time_plotA) ' out of ' num2str(Nt)]))

end
