function [plo]=Cizdir1mod(lon,lat,disp,Nt,Time,per)
NNv=Nt;
lm=int8(1:NNv/8:NNv+1)
A1=squeeze(disp(lm(1),:,:));
A2=squeeze(disp(lm(2),:,:));
A3=squeeze(disp(lm(3),:,:));
A4=squeeze(disp(lm(4),:,:));
A5=squeeze(disp(lm(5),:,:));
A6=squeeze(disp(lm(6),:,:));
A7=squeeze(disp(lm(7),:,:));
A8=squeeze(disp(lm(8),:,:));
%A9=squeeze(disp(lm(9),:,:));

ma=max([max(A1),max(A2),max(A3),max(A4),max(A5),max(A6),max(A7),max(A8)]);
mi=min([min(A1),min(A2),min(A3),min(A4),min(A5),min(A6),min(A7),min(A8)]);


%%%%%%%%%%%%%%%%3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lonn=lon*180/pi;%-0.25;
latt=lat*180/pi-90;%+0.25;
[Plg,Plt]=meshgrid(lonn,latt);


% MAP3D_HAM(Plg,Plt,A1,mi,ma)
% MAP3D_MILLER(Plg,Plt,A1,mi,ma)
% MAP3D_ORTHO(Plg,Plt,A1,mi,ma)
% MAP3D_STEREOGRAPHIC(Plg,Plt,A1,mi,ma)

%figure %(3 components, 6 plots)

%subplot(221)

%MAP3D_HAM(Plg,Plt,A1,mi,ma)
%subplot(222)
%MAP3D_MILLER(Plg,Plt,A1,mi,ma)
%subplot(223)
% MAP3D_ORTHO(Plg,Plt,A1,mi,ma)
%subplot(224)
%MAP3D_STEREOGRAPHIC(Plg,Plt,A1,mi,ma)
%set(get(h,'title'),'string','Potential');


 
 figure
 MAP3D_ORTHO(Plg,Plt,A7,mi,ma)
 %text(0,0,['\lambda=  ',num2str(Time(lm(1))/per*360),'^0'])
 colorbar('h')
 
  
 %figure
 %MAP3D_ORTHO(Plg,Plt,A2,mi,ma)
 %text(0,0,['\lambda=  ',num2str(Time(lm(2))/per*360),'^0'])
 %colorbar('h')
  
 %figure
 %MAP3D_ORTHO(Plg,Plt,A3,mi,ma)
 %text(0,0,['\lambda=  ',num2str(Time(lm(3))/per*360),'^0'])
 %colorbar('h')
  
 %figure
 %MAP3D_ORTHO(Plg,Plt,A4,mi,ma)
 %text(0,0,['\lambda=  ',num2str(Time(lm(4))/per*360),'^0'])
 %colorbar('h')
  
 %figure
 %MAP3D_ORTHO(Plg,Plt,A5,mi,ma)
 %text(0,0,['\lambda=  ',num2str(Time(lm(5))/per*360),'^0'])
 %colorbar('h')
  
 %figure
 %MAP3D_ORTHO(Plg,Plt,A6,mi,ma)
 %text(0,0,['\lambda=  ',num2str(Time(lm(6))/per*360),'^0'])
 %colorbar('h')
 
 %figure
 %MAP3D_ORTHO(Plg,Plt,A7,mi,ma)
 %text(0,0,['\lambda=  ',num2str(Time(lm(7))/per*360),'^0'])
 %colorbar('h')
 
 %figure
 %MAP3D_ORTHO(Plg,Plt,A8,mi,ma)
 %text(0,0,['\lambda=  ',num2str(Time(lm(8))/per*360),'^0'])
 %colorbar('h')
 
 
 
 
 %set(get(h,'title'),'string','Potential');
