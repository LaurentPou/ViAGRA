function [AZ2]=MAP3D_HAM(Plg,Plt,A1,mi,ma)

%lonn=lon*180/pi;%-0.25;
%latt=lat*180/pi-90;%+0.25;
%[Plg,Plt]=meshgrid(lonn,latt);

%figure
m_proj('hammer-aitoff','clongitude',-180);
m_pcolor(Plg,Plt,A1');shading interp; 
hold on;
m_pcolor(Plg-360,Plt,A1');shading interp;
%m_coast('patch',[.6 1 .6]);
m_grid('xaxis','middle');
%h=colorbar('h');
%shading interp
%set(get(h,'title'),'string','Surface displacement');
caxis([mi,ma])
%m_grid('linestyle','none','box','fancy','tickdir','out');
%m_grid('box','fancy','tickdir','in','xaxis','middle')%,'xlabeldir','middle','yticklabels',[-40; -20; 0; 20; 40] );
