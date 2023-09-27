function y0=init(nbdeg)
% --------------------------------------------------------------------------
%initialisation au centre des yi0 avec potentiel et deplacement nuls (y1,3,5)
% --------------------------------------------------------------------------
y0(:,1)=[0;1e-3;0;0;0;0];
y0(:,3)=[0;0;0;1e-3;0;0];
y0(:,2)=[0;0;0;0;0;1e-3]; % cette solution indep (k=3) n'existe plus ds une couche liquide
%car y4=0. a l'interface liquide solide, une nouvelle sol indep est définie. cf interface
%********************************************
