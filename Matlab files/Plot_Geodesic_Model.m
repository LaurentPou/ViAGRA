function Plot_Geodesic_Model(r,rho,g,P,bool_fus,depth_factor,depth_unit)
% Plot the loaded model with the added layers for calculations

rkm = r*depth_factor;

if bool_fus
    
    figure;
    hold on;
    plot(rho*depth_factor,rkm,'r');
    plot(g,rkm,'b');
    plot(P/10^9,rkm,'g');
    legend('Density (kg/cm^3)','Gravity (m/s^2)','Pressure (Gpa)');
    xlabel('Pressure, Density or Gravitational acceleration');
    ylabel(sprintf(['Radius (' depth_unit ')']));
    hold off;

else
    
    figure;
    subplot(1,3,1);
    plot(rho,rkm);
    axis([min(0,min(rho)) max(rho) min(rkm) max(rkm)]);
    xlabel('Density in kg/m^3');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    subplot(1,3,2);
    plot(g,rkm);
    xlabel('Gravity in m/s^2');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    subplot(1,3,3);
    plot(P/10^9,rkm);
    xlabel('Pressure in GPa');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    
end
