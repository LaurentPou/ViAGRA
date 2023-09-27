function [plo]=Cizdir1(lon,lat,disp,Nt,Time,per,SN,mi,ma)

A1=squeeze(disp(SN,:,:));

lonn=lon*180/pi;%-0.25;
latt=(90-lat*180/pi);%lat*180/pi-90;%+0.25;
[Plg,Plt]=meshgrid(lonn,latt);

%
ma=max([max(A1)]);
mi=min([min(A1)]);

%figure
MAP3D_HAM(Plg,Plt,A1,mi,ma)
