function [yi,rp]=decond(out,ref,r,kpot,vpot)

yi(1,:)=ref(1)*out(1,:);
yi(2,:)=ref(3)*out(2,:);
yi(3,:)=ref(1)*out(3,:);
yi(4,:)=ref(3)*out(4,:);
yi(5,:)=ref(1)*ref(4)*out(5,:);
yi(6,:)=ref(4)*out(6,:);

yi=yi*kpot/vpot; %dimenssionement pour le potentiel de marré kpot

rp(1)=0; %rayon corresondant aux yi
rp(2)=((r(2)-r(1))/2);
j=3;
for i=1:size(r,2)-1
    rp(j)=r(i);
    j=j+1;
    rp(j)=r(i)+((r(i+1)-r(i))/2);
    j=j+1;
end
rp(j)=r(i+1);
