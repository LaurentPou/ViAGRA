function retour = InterpCriterion(criterion,dim_values,bool_uplow)
% This function smoothes calculations such as criterion which are heavily
% discretized and thus has a stair-like aspect
% DIM_VALUES is the values of the parameters on which CRITERION must be
% smoothed
% BOOL_UPLOW : 1 if upper_interp wanted, -1 if lower_interp wanted,
% otherwise average of both is given

n = numel(criterion);

% All data come by pair, except the last point

% Upper interpolation
n_upper = floor((n+1)/2);
criterion_upper_interp = zeros(1,n_upper);
upper_values = zeros(1,n_upper);
k = 1;

for i = 1:2:n
    criterion_upper_interp(k) = criterion(i);
    upper_values(k) = dim_values(i);
    k = k+1;
end

% Lower interpolation
n_lower = floor((n)/2);
criterion_lower_interp = zeros(1,n_lower);
lower_values = zeros(1,n_lower);
k = 1;

for i = 2:2:n
    criterion_lower_interp(k) = criterion(i);
    lower_values(k) = dim_values(i);
    k = k+1;
end

upper_interp = interp1(upper_values,criterion_upper_interp,dim_values);
lower_interp = interp1(lower_values,criterion_lower_interp,dim_values);

if bool_uplow == 1
    retour = upper_interp;
elseif bool_uplow == -1
    retour = lower_interp;
else
    retour = 0.5*(upper_interp + lower_interp);
end