% Save Mesh
time_save = Time;
lon_save = lon;
lat_save = lat;

var_name01 = 'time_save';
var_name02 = 'lon_save';
var_name03 = 'lat_save';


file_name01 = sprintf([var_name01 '.txt']);
file_name02 = sprintf([var_name02 '.txt']);
file_name03 = sprintf([var_name03 '.txt']);


save(file_name01, var_name01, '-ASCII','-append');
save(file_name02, var_name02, '-ASCII','-append');
save(file_name03, var_name03, '-ASCII','-append');

% Save stress
% Reminder : dimensions are time, lon, lat, radius
sigma_r_save = srr(:,:,:,radius_plot)*stress_r_factor;
sigma_t_save = stt(:,:,:,radius_plot)*stress_r_factor;
sigma_p_save = spp(:,:,:,radius_plot)*stress_r_factor;
sigma_1_save = s1(:,:,:,radius_plot)*stress_r_factor;
sigma_3_save = s3(:,:,:,radius_plot)*stress_r_factor;

criterion_save = (- criterion(:,:,:,radius_plot) + tau_m(:,:,:,radius_plot))*stress_factor;


var_name1 = 'sigma_r_save';
var_name2 = 'sigma_t_save';
var_name3 = 'sigma_p_save';
var_name4 = 'sigma_1_save';
var_name5 = 'sigma_3_save';
var_name6 = 'criterion_save';

file_name1 = sprintf([var_name1 '.txt']);
file_name2 = sprintf([var_name2 '.txt']);
file_name3 = sprintf([var_name3 '.txt']);
file_name4 = sprintf([var_name4 '.txt']);
file_name5 = sprintf([var_name5 '.txt']);
file_name6 = sprintf([var_name6 '.txt']);

save(file_name1, var_name1, '-ASCII','-append');
save(file_name2, var_name2, '-ASCII','-append');
save(file_name3, var_name3, '-ASCII','-append');
save(file_name4, var_name4, '-ASCII','-append');
save(file_name5, var_name5, '-ASCII','-append');
save(file_name6, var_name6, '-ASCII','-append');

% Save depth of failure
if exist('frictionbis','var');
    max_deep_failure_save = (r_s(end)-r_s(min_failure))*depth_factor;
    min_deep_failure_save = (r_s(end)-r_s(max_failure2))*depth_factor;
    max_surface_failure_save = (r_s(end)-r_s(min_failure2))*depth_factor;
    min_surface_failure_save = (r_s(end)-r_s(max_failure))*depth_factor;
    friction_angle_save = frictionbis*180/pi;

    var_name10 = 'max_deep_failure_save';
    var_name20 = 'min_deep_failure_save';
    var_name30 = 'max_surface_failure_save';
    var_name40 = 'min_surface_failure_save';
    var_name50 = 'friction_angle_save';

    file_name10 = sprintf([var_name10 '.txt']);
    file_name20 = sprintf([var_name20 '.txt']);
    file_name30 = sprintf([var_name30 '.txt']);
    file_name40 = sprintf([var_name40 '.txt']);
    file_name50 = sprintf([var_name50 '.txt']);

    save(file_name10, var_name10, '-ASCII','-append');
    save(file_name20, var_name20, '-ASCII','-append');
    save(file_name30, var_name30, '-ASCII','-append');
    save(file_name40, var_name40, '-ASCII','-append');
    save(file_name50, var_name50, '-ASCII','-append');
    
    
end
