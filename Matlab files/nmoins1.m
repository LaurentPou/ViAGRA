function C=nmoins1(y0,sc,I)
% --------------------------------------------------------------------------
% conditionnement des yi,0 pour chaque couche de la super couche
% --------------------------------------------------------------------------
C1=[];
C=[];
 for i=I(1,sc)+1:I(1,sc+1)
  if i==I(1,sc)+1 %initialisation pour recupérer les valeur de la super couche précédente via y0  
     C0=[y0(1);0;0;y0(2);0;0;y0(3);0;0;y0(4);0;0;y0(5);0;0;y0(6);0;0];% amodifier si ppc ~=3
  else 
     C0=zeros(18,1);
  end
  C=[C;C0];
 end 
%********************************************
