function S=solid(G,lambda,mu,rho,g,r,D,n,w,cst1,cst2,ref)
% --------------------------------------------------------------------------
% Matrice du problème pour une couche solide
% --------------------------------------------------------------------------
clear S;
%soit p les positions des yij, j=1:3 fonction de r et D
p(3)=r;
p(2)=p(3)-D/2;
p(1)=p(2)-D/2;

Zero=[0,0,0];
%remplissage de la matrice
%Y1
s101=[1,0,0];
s10=[s101, Zero,  Zero,  Zero,  Zero,  Zero];

s111=[1/D,-2*lambda/((lambda+2*mu)*p(2)),-1/D];
s112=[0,1/(lambda+2*mu),0];
s113=[0,(lambda*n*(n+1))/((lambda+2*mu)*p(2)),0];
s11=[s111,s112,s113,Zero,Zero,Zero];

s121=[0,1/(D/2),-1/(D/2)-2*lambda/((lambda+2*mu)*p(3))];
s122=[0,0,1/(lambda+2*mu)];
s123=[0,0,(lambda*n*(n+1))/((lambda+2*mu)*p(3))];
s12=[s121,s122,s123,Zero,Zero,Zero];

%Y2
s202=[1,0,0];
s20=[Zero,s202,Zero,Zero,Zero,Zero];

s211=[0,(1/(p(2)^2))*(-(cst1*w^2*rho*p(2)^2)-cst1*4*rho*g*p(2)+(4*mu*(3*lambda+2*mu))/(lambda+2*mu)),0];
s212=[1/D,-4*mu/((lambda+2*mu)*p(2)),-1/D];
s213=[0,(n*(n+1)/(p(2)^2))*(cst1*rho*g*p(2)-(2*mu*(3*lambda+2*mu))/(lambda+2*mu)),0];
s214=[0,(n*(n+1))/p(2),0];
s216=[0,-cst1*rho,0];
s21=[s211,s212,s213,s214,Zero,s216];

s221=[0,0,(1/(p(3)^2))*(-(cst1*w^2*rho*p(3)^2)-cst1*4*rho*g*p(3)+(4*mu*(3*lambda+2*mu))/(lambda+2*mu))];
s222=[0,1/(D/2),-1/(D/2)-4*mu/((lambda+2*mu)*p(3))];
s223=[0,0,(1/(p(3)^2))*(cst1*n*(n+1)*rho*g*p(3)-(2*mu*(3*lambda+2*mu)*(n+1)*n)/(lambda+2*mu))];
s224=[0,0,(n*(n+1))/p(3)];
s226=[0,0,-cst1*rho];
s22=[s221,s222,s223,s224,Zero,s226];

%Y3
s303=[1,0,0];
s30=[Zero,Zero,s303,Zero,Zero,Zero];

s311=[0,-1/p(2),0];
s313=[1/D,1/p(2),-1/D];
s314=[0,1/mu,0];
s31=[s311,Zero,s313,s314,Zero,Zero];

s321=[0,0,-1/p(3)];
s323=[0,1/(D/2),-1/(D/2)+1/p(3)];
s324=[0,0,1/mu];
s32=[s321,Zero,s323,s324,Zero,Zero];

%Y4
s404=[1,0,0];
s40=[Zero,Zero,Zero,s404,Zero,Zero];

s411=[0,(1/(p(2)^2))*(cst1*rho*g*p(2)-(2*mu*(3*lambda+2*mu))/(lambda+2*mu)),0];
s412=[0,-lambda/(p(2)*(lambda+2*mu)),0];
s413=[0,(1/(p(2)^2))*(-cst1*rho*w^2*p(2)^2+(2*mu/(lambda+2*mu)*(lambda*(2*n^2+2*n-1)+2*mu*(n^2+n-1)))),0];
s414=[1/D,-3/(p(2)),-1/D];
s415=[0,-cst1*rho/p(2),0];
s41=[s411,s412,s413,s414,s415,Zero];

s421=[0,0,(1/(p(3)^2))*(cst1*rho*g*p(3)-(2*mu*(3*lambda+2*mu))/(lambda+2*mu))];
s422=[0,0,-lambda/(p(3)*(lambda+2*mu))];
s423=[0,0,(1/(p(3)^2))*(-cst1*rho*w^2*p(3)^2+(2*mu/(lambda+2*mu))*(lambda*(2*n^2+2*n-1)+2*mu*(n^2+n-1)))];
s424=[0,1/(D/2),-1/(D/2)-3/(p(3))];
s425=[0,0,-cst1*rho/p(3)];
s42=[s421,s422,s423,s424,s425,Zero];

%Y5

s505=[1,0,0];
s50=[Zero,Zero,Zero,Zero,s505,Zero];

s511=[0,cst2*rho,0];
s515=[1/D,0,-1/D];
s516=[0,1,0];
s51=[s511,Zero,Zero,Zero,s515,s516];

s521=[0,0,cst2*rho];
s525=[0,1/(D/2),-1/(D/2)];
s526=[0,0,1];
s52=[s521,Zero,Zero,Zero,s525,s526];

%Y6

s606=[1,0,0];
s60=[Zero,Zero,Zero,Zero,Zero,s606];

s613=[0,-(cst2*rho*n*(n+1))/p(2),0];
s615=[0,n*(n+1)/(p(2)^2),0];
s616=[1/D,-2/p(2),-1/D];
s61=[Zero,Zero,s613,Zero,s615,s616];

s623=[0,0,-(cst2*rho*n*(n+1))/p(3)];
s625=[0,0,n*(n+1)/(p(3)^2)];
s626=[0,1/(D/2),-1/(D/2)-2/p(3)];
s62=[Zero,Zero,s623,Zero,s625,s626];

%concatenation
S=[s10;s11;s12;s20;s21;s22;s30;s31;s32;s40;s41;s42;s50;s51;s52;s60;s61;s62];
S=sparse(S);
%********************************************

