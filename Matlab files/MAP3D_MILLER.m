function [AZ2]=MAP3D_STEREOGRAPHIC(Plg,Plt,A1,mi,ma)



m_proj('miller','lat',82);
m_pcolor(Plg,Plt,A1');shading interp; 
hold on;
m_pcolor(Plg-360,Plt,A1');shading interp;
%m_coast('color',[0 .6 0]);
%m_line(lon,lat,'linewi',3,'color','r');
m_grid('linestyle','none','box','fancy','tickdir','out');
caxis([mi,ma])
