function param_fit = PolyInterpolation(param,var_ref,var_eval,degree)


p = polyfit(var_ref,param,degree);
param_fit = polyval(p,var_eval);
