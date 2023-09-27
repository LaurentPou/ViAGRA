function [AZ2]=MAP3D_STEREOGRAPHIC(Plg,Plt,A1,mi,ma)


m_proj('stereographic','lat',90,'long',0,'radius',25); %90,30,25

m_pcolor(Plg,Plt,A1');shading interp; 
hold on;
m_pcolor(Plg-360,Plt,A1');shading interp;
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
caxis([mi,ma])

