function I=locinterface(mu,nbcouches)
% --------------------------------------------------------------------------
%I revoye les valeurs des couches aux interfaces L/S et S/L (localisation)
% --------------------------------------------------------------------------
temp(1,:)=(mu==0); % 1 si mu(i)=0, 0 si mu(i)~=0
j=2;
I(1,1)=1;
for i=2:nbcouches
    if temp(1,i)~=temp(1,i-1); % interface
        I(1,j)=i;
        j=j+1;
    end
end
I=I-1;
I(1,j)=nbcouches;
I(2,1)=j-1; %nb de "super couche"
%********************************************
