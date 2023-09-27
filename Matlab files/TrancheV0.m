function [r1,rho1,lambda1,mu1,V1]=TrancheV0(r,rho,lambda,mu,V,nbtranche,bool_interp_tranche)
% TrancheV0 needs a point at the center
% Bool is true for linear interpolation, and wrong for step functions

n=numel(r);
ii=1;

for i=1:n-1

%     ii_start = ii;
    
    % Linear interpolation
    for j=1:nbtranche
        
        if bool_interp_tranche
            
            r1(ii)=j*(r(i+1)-r(i))/nbtranche+r(i);
            rho1(ii)=j*(rho(i+1)-rho(i))/nbtranche+rho(i);
            lambda1(ii)=j*(lambda(i+1)-lambda(i))/nbtranche+lambda(i);
            mu1(ii)=j*(mu(i+1)-mu(i))/nbtranche+mu(i);
            V1(ii)=j*(V(i+1)-V(i))/nbtranche+V(i);
        
        else
            
            r1(ii)=j*(r(i+1)-r(i))/nbtranche+r(i);
            rho1(ii)=rho(i);
            lambda1(ii)=lambda(i);
            mu1(ii)=mu(i);
            V1(ii)=V(i);
            
        end
        
        
        ii = ii+1;
        
    end
    
%     ii_end = ii-1;
    
%     % Smoothen discontinuity
%     if i > 1 && r(i)-r(i-1) <= 1000
%         
%         r_value = [r1(i) r1(floor((ii_end-ii_start)/3+ii_start)) r1(floor(2*(ii_end-ii_start)/3+ii_start)) r1(ii_end)];
%         rho_value = [rho1(i) rho1(floor((ii_end-ii_start)/3+ii_start)) rho1(floor(2*(ii_end-ii_start)/3+ii_start)) rho1(ii_end)];
%         lambda_value = [lambda1(i) lambda1(floor((ii_end-ii_start)/3+ii_start)) lambda1(floor(2*(ii_end-ii_start)/3+ii_start)) lambda1(ii_end)];
%         mu_value = [mu1(i) mu1(floor((ii_end-ii_start)/3+ii_start)) mu1(floor(2*(ii_end-ii_start)/3+ii_start)) mu1(ii_end)];
%         V_value = [V1(i) V1(floor((ii_end-ii_start)/3+ii_start)) V1(floor(2*(ii_end-ii_start)/3+ii_start)) V1(ii_end)];
%         
%         r_eval = r1(ii_start:ii_end);
%         
% %         rho_p = polyfit(r_value,rho_value,2);
% %         lambda_p = polyfit(r_value,lambda_value,2);
% %         mu_p = polyfit(r_value,mu_value,2);
% %         V_p = polyfit(r_value,V_value,2);
%         
%         rho1(ii_start:ii_end) = interp1(r_value,rho_value,r_eval,'spline');
%         lambda1(ii_start:ii_end) = interp1(r_value,lambda_value,r_eval,'spline');
%         mu1(ii_start:ii_end) = interp1(r_value,mu_value,r_eval,'spline');
%         V1(ii_start:ii_end) = interp1(r_value,V_value,r_eval,'spline');
%         
%         
%     end
    
end
