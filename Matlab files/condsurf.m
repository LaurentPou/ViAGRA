function [coeff,determinant]=condsurf(Y,n,cst3,appli_pression,appli_potentiel)
% --------------------------------------------------------------------------
% application des conditions de surfaces (Saiko) pour le calcul des coeff de la combi lineaire de Yi
%determiné par la distribution de masse sigme=
% --------------------------------------------------------------------------


%IMPORTANT
%**************************************
%pression: (k2') conditions de surface sur y2
if appli_pression==1
cond1=-(2*n+1)*cst3;
else
cond1=0;
end
%**************************************
%charge= pression et potentiel
%**************************************
%potentiel: (k2) conditons de surface sur y5 et y6
if appli_potentiel==1
cond3=(2*n+1);
else
cond3=0;
end
%**************************************


l=size(Y,1);
a=[l-15,l-12,l-9,l-6,l-3,l];%ligne de y1:6 pour la derniere couche
SURF1=[Y(a(2),1),Y(a(2),2),Y(a(2),3)];
SURF2=[Y(a(4),1),Y(a(4),2),Y(a(4),3)];
SURF3=[Y(a(6),1)+(n+1)*Y(a(5),1),   Y(a(6),2)+(n+1)*Y(a(5),2),  Y(a(6),3)+(n+1)*Y(a(5),3)];
SURF=[SURF1;SURF2;SURF3];
COND=[cond1;0;cond3];
coeff=SURF\COND;
determinant=det(inv(SURF));


for i=1:6
Ysurf(i)=coeff(1)*Y(a(i),1)+coeff(2)*Y(a(i),2)+coeff(3)*Y(a(i),3);%la solution est une combinaison linéaire des 3 sol independante
end
%********************************************

