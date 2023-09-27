x_values = linspace(c_m(tt,ix,j,radius)-tau_m(tt,ix,j,radius),c_m(tt,ix,j,radius)+tau_m(tt,ix,j,radius),1000);
t = 0:pi/360:2*pi;


figure;
hold on;

scatter(c_m(tt,ix,j,radius),0); % Center of MC circle

plot(tau_m(tt,ix,j,radius)*cos(t) + c_m(tt,ix,j,radius),tau_m(tt,ix,j,radius)*sin(t),':'); % MC circle

plot(x_values,cohe+x_values*tan(friction)); % MC criterion

xline(0,'LineWidth',3,'Color','k');
yline(0,'LineWidth',3,'Color','k');