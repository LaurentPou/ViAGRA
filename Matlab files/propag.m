function M=propag(G,lambda,mu,rho,g,r,D,n,I,sc,w,cst1,cst2,ref)
% --------------------------------------------------------------------------
% Inversion de la matrice de la "supercouche"C jusqu'a la prochaine interface S/L,L/S,ou surface
% --------------------------------------------------------------------------
T=[];
T1=[];
T2=[];


for i=I(1,sc)+1:I(1,sc+1)%on assemble chaque sous couche de la super couche
    if mu(i)~=0                                                        %couche solide
        T=solid(G,lambda(i),mu(i),rho(i),g(i),r(i),D(i),n,w,cst1,cst2,ref);
    else                                                               %couche liquide
        T=liquid(G,lambda(i),rho(i),g(i),r(i),D(i),n,w,cst1,cst2,ref);
    end
    T1=blkdiag(T1,T);%association des systems pour une super couche(=concaténation des block solide ou liquide selon la diagonale)
    %if i/50==floor(i/50)
    % nombre_couches=i%suivi consol
    %end
end

w1=size(T,1);
if (I(1,sc)+1-I(1,sc+1))~=0 %si il y a plus d'une couche, le yi0=y(i-1),1
    V=Vgen;
    for i=I(1,sc)+1:I(1,sc+1)
        T2=blkdiag(T2,V);
    end
    w=zeros(w1,size(T2,2));
    T2=[w;T2];
    T2((size(T2,1)-17):size(T2,1),:)=[];
    T1=T1+T2;
end

M=sparse(T1); % libère de l'espace mémoire pour les termes non nuls de T1--> accélère nettenment l'inversion

function V=Vgen
v=[-1,0,0,-1,0,0,-1,0,0,-1,0,0,-1,0,0,-1];%génération d'une matrice 18*18 pour yi0=yi-1,1
   V=diag(v);
   v=zeros(18-2,2); 
   V=[v,V];
   v=zeros(2,18);
   V=[V;v];
   V=sparse(V);
 