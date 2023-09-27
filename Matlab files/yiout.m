function out=yiout(Y,COEFF,I,mu)%réordonne les yi le long de r, yi1 et yi2 pour chaque couche et yi3 pour la derniere couche
A=COEFF(1);
B=COEFF(2);
C1=COEFF(3);%calculé en surface
for sc=1:I(2,1)
    if mu(I(1,sc+1))~=0%calculé sous une couche liquide
    C2(sc)=-Y(I(1,sc+1)*18-6,1)/Y(I(1,sc+1)*18-6,3)*A-Y(I(1,sc+1)*18-6,2)/Y(I(1,sc+1)*18-6,3)*B;
    end
end

for sc=1:I(2,1)
    for i=I(1,sc)+1:I(1,sc+1)
        i=i-1;
        for j=1:2
            if mu(i+1)~=0%couche solide
                if sc==I(2,1)  %il n'existe plus de couche liquide au dessus
                    C=C1;
                else           %il existe une couche liquide au dessus
                    C=C2(sc);
                end
                
                y1(j+i*2)=A*Y(j+i*18,1)   +B*Y(j+i*18,2)   +C*Y(j+i*18,3);
                y2(j+i*2)=A*Y(3+j+i*18,1) +B*Y(3+j+i*18,2) +C*Y(3+j+i*18,3);
                y3(j+i*2)=A*Y(6+j+i*18,1) +B*Y(6+j+i*18,2) +C*Y(6+j+i*18,3);
                y4(j+i*2)=A*Y(9+j+i*18,1) +B*Y(9+j+i*18,2) +C*Y(9+j+i*18,3);
                y5(j+i*2)=A*Y(12+j+i*18,1)+B*Y(12+j+i*18,2)+C*Y(12+j+i*18,3);
                y6(j+i*2)=A*Y(15+j+i*18,1)+B*Y(15+j+i*18,2)+C*Y(15+j+i*18,3);
                
            else %couche liquide
                
                y1(j+i*2)=A*Y(j+i*18,1)   +B*Y(j+i*18,2)   ;
                y2(j+i*2)=A*Y(3+j+i*18,1) +B*Y(3+j+i*18,2) ;
                y3(j+i*2)=A*Y(6+j+i*18,1) +B*Y(6+j+i*18,2) ;
                y4(j+i*2)=A*Y(9+j+i*18,1) +B*Y(9+j+i*18,2) ;
                y5(j+i*2)=A*Y(12+j+i*18,1)+B*Y(12+j+i*18,2);
                y6(j+i*2)=A*Y(15+j+i*18,1)+B*Y(15+j+i*18,2);
                
            end
        end
    end
end

r=I(1,I(2,1)+1)*18;
t=size(y6,2);
y1(t+1)=A*Y(r-15,1)+B*Y(r-15,2)+C1*Y(r-15,3);
y2(t+1)=A*Y(r-12,1)+B*Y(r-12,2)+C1*Y(r-12,3);
y3(t+1)=A*Y(r-9,1) +B*Y(r-9,2) +C1*Y(r-9,3);
y4(t+1)=A*Y(r-6,1) +B*Y(r-6,2) +C1*Y(r-6,3);
y5(t+1)=A*Y(r-3,1) +B*Y(r-3,2) +C1*Y(r-3,3);
y6(t+1)=A*Y(r,1)   +B*Y(r,2)   +C1*Y(r,3);


out=[y1;y2;y3;y4;y5;y6];
