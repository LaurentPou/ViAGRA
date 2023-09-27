function g=gravity(rho,r,nbcouches,G)
% --------------------------------------------------------------------------
% Calcul de g pour chaques couches (et de D pour la subroutine propag)
% --------------------------------------------------------------------------
g(1)=(4/3)*pi*G*rho(1)*r(1);
for i=2:nbcouches
    g(i)=g(i-1)*(r(i-1)/r(i))^2+(4/3)*pi*G*rho(i)*(r(i)-r(i-1)*((r(i-1)/r(i))^2));
end
%********************************************
