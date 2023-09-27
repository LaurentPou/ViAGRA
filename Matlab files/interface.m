function y0=interface(Y,I,mu,sc)
% --------------------------------------------------------------------------
% calcul les yi,0 de la super couche suivante a partir des yi,1 calculé de la couche précédente
% --------------------------------------------------------------------------
l=size(Y,1);
a=[l-15,l-12,l-9,l-6,l-3,l];%ligne de y1:6 pour la derniere couche
if mu(I(1,sc+1))==0 %interface liquide/solide
    %interface='liquide/solide'
    for k=1:2
        y0(1,k)=Y(a(1),k); %y1,y2,y5,y6 sont continus a l'interface
        y0(2,k)=Y(a(2),k);
        y0(3,k)=0;
        y0(5,k)=Y(a(5),k);
        y0(6,k)=Y(a(6),k);
        
    end
    y0(1,3)=0; % on redéfini une solution indep dans la couche solide
    y0(2,3)=0; % on redéfini une solution indep dans la couche solide
    y0(3,3)=1e-3; % on redéfini une solution indep dans la couche solide
    y0(4,3)=0; % on redéfini une solution indep dans la couche solide
    y0(5,3)=0; % on redéfini une solution indep dans la couche solide
    y0(6,3)=0; % on redéfini une solution indep dans la couche solide
    
    
else   %interface  solide/liquide, on recombine les trois solutions indep en deux solution qui satisassent y4 nul
    %interface='solide/liquide'
    K1=Y(a(4),1)/Y(a(4),3);
    K2=Y(a(4),2)/Y(a(4),3);
    
    for i=1:6
        y0(i,1)=Y(a(i),1)-K1*Y(a(i),3);%cette combinaison satisfait y4=0
        y0(i,2)=Y(a(i),2)-K2*Y(a(i),3);%cette combinaison satisfait y4=0
        y0(i,3)=0;  %la troisieme sol indep n'existe plus dans la couche liquide
    end
end
%********************************************

