function Plot_Mechanical_Model(r,rho,lambda,mu,V,bool_fus,depth_factor,depth_unit)
% Plot the loaded model with the added layers for calculations

rkm = r*depth_factor;

if bool_fus
    
    figure;
    subplot(1,2,1)
    hold on;
    plot(lambda/10^10,rkm,'r');
    plot(mu/10^10,rkm,'b');
    plot(rho*depth_factor,rkm,'g');
    legend('\lambda','\mu','density');
    xlabel('\lambda or \mu in 10*GPa, or density in kg/cm^3');
    ylabel(sprintf(['Radius (' depth_unit ')']));
    hold off;
    
    subplot(1,2,2);
    plot(V,rkm);
    xlabel('\nu Dynamic viscosity in GPa.s');
    ylabel(sprintf(['Radius (' depth_unit ')']));
    
else

    figure;
    subplot(1,4,1);
    plot(rho,rkm);
    axis([min(0,min(rho)) max(rho) min(rkm) max(rkm)]);
    xlabel('Density in kg/m^3');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    subplot(1,4,2);
    plot(lambda/10^9,rkm);
    xlabel('\lambda in GPa');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    subplot(1,4,3);
    plot(mu/10^9,rkm);
    xlabel('\mu in GPa');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    subplot(1,4,4);
    semilogx(V/10^9,rkm);
    xlabel('\nu Dynamic viscosity in GPa.s');
    ylabel(sprintf(['Radius (' depth_unit ')']));

end



