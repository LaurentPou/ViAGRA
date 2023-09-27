function L=liquid(G,lambda,rho,g,r,D,n,w,cst1,cst2,ref)
% --------------------------------------------------------------------------
% Remplissage de la matrice du problème pour une couche liquide
% --------------------------------------------------------------------------
clear L;
%soit p les positions des yij, j=1:3 fonction de r et D
p(3)=r;
p(2)=p(3)-D/2;
Zero=[0,0,0];

%remplissage de la matrice du système ODE (yi)
%Y1
l101=[1,0,0];
l10=[l101,Zero,Zero,Zero,Zero,Zero];

l111=[1/D,-2/p(2),-1/D];
l112=[0,1/lambda,0];
l113=[0,(n*(n+1))/p(2),0];
l11=[l111,l112,l113,Zero,Zero,Zero];

l121=[0,2/D,-2/D-2/p(3)];
l122=[0,0,1/lambda];
l123=[0,0,(n*(n+1))/p(3)];
l12=[l121,l122,l123,Zero,Zero,Zero];
%Y2
l202=[1,0,0];
l20=[Zero,l202,Zero,Zero,Zero,Zero];

l211=[0,-cst1*w^2*rho-4*rho*g*cst1/p(2),0];
l212=[1/D,0,-1/D];
l213=[0,cst1*n*(n+1)*rho*g/p(2),0];
l216=[0,-cst1*rho,0];
l21=[l211,l212,l213,Zero,Zero,l216];

l221=[0,0,-cst1*w^2*rho-4*rho*g*cst1/p(2)];
l222=[0,2/D,-2/D];
l223=[0,0,cst1*n*(n+1)*rho*g/p(2)];
l226=[0,0,-cst1*rho];
l22=[l221,l222,l223,Zero,Zero,l226];

%Y3
l303=[1,0,0];
l30=[Zero,Zero,l303,Zero,Zero,Zero];

if w==0 %y3 indeterminé, fixé a 0
    l313=[0,1,0];
    l31=[Zero,Zero,l313,Zero,Zero,Zero];

    l323=[0,0,1];
    l32=[Zero,Zero,l323,Zero,Zero,Zero];
else
    %l311=[0,g/(p(2)*w^2),0];% on ne calcul pas de point intermediaire, car le syst deviens instable
    %l312=[0,-1/(p(2)*w^2*rho*cst1),0];
    l313=[0,-1,1];
    %l315=[0,-1/(p(2)*w^2),0];
    l31=[Zero,Zero,l313,Zero,Zero,Zero];

    l321=[0,0,g/p(3)];
    l322=[0,0,-1/(p(3)*rho*cst1)];
    l323=[0,0,-w^2];
    l325=[0,0,-1/p(3)];
    l32=[l321,l322,l323,Zero,l325,Zero];
end
   
%Y4 %nul dans la couche liquide
l40=[Zero,Zero,Zero,1,0,0,Zero,Zero];
l41=[Zero,Zero,Zero,0,1,0,Zero,Zero];
l42=[Zero,Zero,Zero,0,0,1,Zero,Zero];

%Y5
l505=[1,0,0];
l50=[Zero,Zero,Zero,Zero,l505,Zero];

l511=[0,cst2*rho,0];
l515=[1/D,0,-1/D];
l516=[0,1,0];
l51=[l511,Zero,Zero,Zero,l515,l516];

l521=[0,0,cst2*rho];
l525=[0,2/D,-2/D];
l526=[0,0,1];
l52=[l521,Zero,Zero,Zero,l525,l526];

%Y6
l606=[1,0,0];
l60=[Zero,Zero,Zero,Zero,Zero,l606];

l613=[0,-(cst2*rho*n*(n+1))/p(2),0];
l615=[0,n*(n+1)/(p(2)^2),0];
l616=[1/D,-2/p(2),-1/D];
l61=[Zero,Zero,l613,Zero,l615,l616];

l623=[0,0,-(cst2*rho*n*(n+1))/p(3)];
l625=[0,0,n*(n+1)/(p(3)^2)];
l626=[0,2/D,-2/D-2/p(3)];
l62=[Zero,Zero,l623,Zero,l625,l626];


L=[l10;l11;l12;l20;l21;l22;l30;l31;l32;l40;l41;l42;l50;l51;l52;l60;l61;l62];
L=sparse(L);
