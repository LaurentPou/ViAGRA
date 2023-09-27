function [AZ2]=MAP3D_ORTHO(Plg,Plt,A1,mi,ma)

%lonn=lon*180/pi;%-0.25;
%latt=lat*180/pi-90;%+0.25;
%[Plg,Plt]=meshgrid(lonn,latt);

m_proj('ortho','lat',40','long',-180');
m_pcolor(Plg,Plt,A1');shading interp; 
hold on;
m_pcolor(Plg-360,Plt,A1');shading interp;
%m_coast('patch','r'); %
m_grid('linest','-','xticklabels',[],'yticklabels',[]);
%patch(.55*[-1 1 1 -1],.25*[-1 -1 1 1]-.55,'w');
%text(0,-.55,'M\_Map','fontsize',25,'color','b',...%
%   'vertical','middle','horizontal','center');
%set(gcf,'units','inches','position',[2 2 3 3]);
%set(gcf,'paperposition',[3 3 3 3]);

caxis([mi,ma])