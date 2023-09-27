function Plot_Seismic_Model(r,rho,lambda,mu,V,per,bool_fus,depth_factor,depth_unit)
% Plot the loaded model with the added layers for calculations

% Radius in kilometers for lisibility
rkm = r*depth_factor;

% Conversion of lambda and mu into Vs and Vp, in km/s
Vp = sqrt((lambda+2*mu)./rho)/1000;
Vs = sqrt(mu./rho)/1000;

% Conversion of lambda, mu and per into Q
K = lambda+2/3*mu;
Q = V*(2*pi/per)./K;

if bool_fus
    
    figure;
    subplot(1,3,1)
    hold on;
    plot(Vs,rkm,'r');
    plot(Vp,rkm,'b');
    legend('Vs','Vp');
    xlabel('Seismic velocity in km/s');
    ylabel(sprintf(['Radius (' depth_unit ')']));
    hold off;

    subplot(1,3,2);
    plot(rho*depth_factor,rkm,'g');
    xlabel('Density in kg/cm^3');
    ylabel(sprintf(['Radius (' depth_unit ')']));
    
    subplot(1,3,3);
    semilogx(Q,rkm);
    xlabel('Quality factor Q');
    ylabel(sprintf(['Radius (' depth_unit ')']));
    
else
    
    figure;
    subplot(1,4,1);
    plot(rho,rkm);
    axis([min(0,min(rho)) max(rho) min(rkm) max(rkm)]);
    xlabel('Density in kg/m^3');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    subplot(1,4,2);
    plot(Vp,rkm);
    xlabel('Vp in km/s');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    subplot(1,4,3);
    plot(Vs,rkm);
    xlabel('Vs in km/s');
    ylabel(sprintf(['Radius (' depth_unit ')']));

    subplot(1,4,4);
    plot(Q,rkm);
    xlabel('Quality factor Q');
    ylabel(sprintf(['Radius (' depth_unit ')']));
    
end


