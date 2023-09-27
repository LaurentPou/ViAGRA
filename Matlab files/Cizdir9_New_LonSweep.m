function [plo]=Cizdir9_New_LonSweep(lon,lat,dat_i,Nt,time_plot,per,name,bool_color_num,bool_save,model_name,bool_maps_type)
% plot at a given time, longitude sweep
% bool_maps_type = bool array for maps. Order is: ham (trinkel proj),
% miller, ortho, ortho 3D, stereographic (north pole)

bool_ham = bool_maps_type(1);
bool_miller = bool_maps_type(2);
bool_ortho = bool_maps_type(3);
bool_ortho_3D = bool_maps_type(4);
bool_stereo = bool_maps_type(5);

dat_mid = squeeze(dat_i(time_plot,:,:)); % Size: Nlon, Nlat
Nlon = numel(lon);
lon_idx = 1:Nlon;
lon_2idx = cat(2,lon_idx,lon_idx);


NNv=floor((Nlon)/2); % Nt/2
%lm=int8(1:NNv/8:NNv+1); % % NNv/8
lm=(1:floor(NNv/8):NNv+1); % % NNv/8
A1=Cizdir9_Newmod(dat_mid,lm(1),lon_2idx,Nlon);
A2=Cizdir9_Newmod(dat_mid,lm(2),lon_2idx,Nlon);
A3=Cizdir9_Newmod(dat_mid,lm(3),lon_2idx,Nlon);
A4=Cizdir9_Newmod(dat_mid,lm(4),lon_2idx,Nlon);
A5=Cizdir9_Newmod(dat_mid,lm(5),lon_2idx,Nlon);
A6=Cizdir9_Newmod(dat_mid,lm(6),lon_2idx,Nlon);
A7=Cizdir9_Newmod(dat_mid,lm(7),lon_2idx,Nlon);
A8=Cizdir9_Newmod(dat_mid,lm(8),lon_2idx,Nlon);
A9=Cizdir9_Newmod(dat_mid,lm(9),lon_2idx,Nlon);

% NNv=floor(numel(lon)/2);
% lm=int8(1:NNv/8:NNv+1);
% A1=squeeze(dat(:,lm(1),:));
% A2=squeeze(dat(:,lm(2),:));
% A3=squeeze(dat(:,lm(3),:));
% A4=squeeze(dat(:,lm(4),:));
% A5=squeeze(dat(:,lm(5),:));
% A6=squeeze(dat(:,lm(6),:));
% A7=squeeze(dat(:,lm(7),:));
% A8=squeeze(dat(:,lm(8),:));
% A9=squeeze(dat(:,lm(9),:));

ma=max([max(A1),max(A2),max(A3),max(A4),max(A5),max(A6),max(A7),max(A8),max(A9)]);
ma1 = max(max(A1));
mi=min([min(A1),min(A2),min(A3),min(A4),min(A5),min(A6),min(A7),min(A8),min(A9)]);
mi1 = min(min(A2));

%%%%%%%%%%%%%%%%3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lonn=lon*180/pi;%-0.25;
latt=(90-lat*180/pi);%lat*180/pi-90;%+0.25;
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

% Legend to indicates if criterion is reached or not
bool_criterion_find = ~isempty(strfind(name,'Failure criterion'));
if ~bool_color_num
    if ma > 0 && mi > 0
        array_bar = [0,mi,ma];
        label_bar = {'Reached','Not reached','Not reached at all'};
    elseif  ma < 0 && mi < 0
        array_bar = [mi,ma,0];
        label_bar = {'Well reached','Reached','Barely reached'};
    else
        array_bar = [mi,0,ma];
        label_bar = {'Well reached','Reached','Not reached'};
    end
end


% Before : Time(lm(1))/per*360 instead of lon(lm(1))*180/pi

%%
if bool_ham

    figure %17
    subplot(331)
    MAP3D_HAM(Plg,Plt,A1,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(1))*180/pi),'^o'])
    subplot(332)
    MAP3D_HAM(Plg,Plt,A2,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(2))*180/pi),'^o'])
    subplot(333)
    MAP3D_HAM(Plg,Plt,A3,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(3))*180/pi),'^o'])
    subplot(334)
    MAP3D_HAM(Plg,Plt,A4,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(4))*180/pi),'^o'])
    subplot(335)
    MAP3D_HAM(Plg,Plt,A5,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(5))*180/pi),'^o'])
    subplot(336)
    MAP3D_HAM(Plg,Plt,A6,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(6))*180/pi),'^o'])
    subplot(337)
    MAP3D_HAM(Plg,Plt,A7,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(7))*180/pi),'^o'])
    subplot(338)
    MAP3D_HAM(Plg,Plt,A8,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(8))*180/pi),'^o'])
    subplot(339)
    MAP3D_HAM(Plg,Plt,A9,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(9))*180/pi),'^o'])
    c = colorbar;
    if bool_criterion_find && ~bool_color_num
        c = colorbar('Ticks',array_bar,'TickLabels',label_bar);
    end
    set(c,'position',[0.9 0.08 0.025 0.9]); % x pos, y pos, largeur, ?
    ylabel(c,sprintf(name),'FontSize',15);
    suptitle(sprintf(['Map of \n' name ' over half its surface']));

end

% figure % 18
% subplot(331)
%  MAP3D_MILLER(Plg,Plt,A1,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(1))*180/pi),'^o'])
% subplot(332)
%  MAP3D_MILLER(Plg,Plt,A2,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(2))*180/pi),'^o'])
% subplot(333)
%  MAP3D_MILLER(Plg,Plt,A3,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(3))*180/pi),'^o'])
% subplot(334)
%  MAP3D_MILLER(Plg,Plt,A4,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(4))*180/pi),'^o'])
% subplot(335)
%  MAP3D_MILLER(Plg,Plt,A5,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(5))*180/pi),'^o'])
% subplot(336)
%  MAP3D_ORTHO(Plg,Plt,A6,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(6))*180/pi),'^o'])
% subplot(337)
%  MAP3D_ORTHO(Plg,Plt,A7,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(7))*180/pi),'^o'])
%  subplot(338)
%  MAP3D_ORTHO(Plg,Plt,A8,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(8))*180/pi),'^o'])
%  subplot(339)
%  MAP3D_STEREOGRAPHIC(Plg,Plt,A9,mi,ma)
%  text(0,0,['\lambda=  ',num2str(lon(lm(9))*180/pi),'^o'])
%  colorbar('h')
%

%%
if bool_miller

    figure % 19
    subplot(331)
    MAP3D_MILLER(Plg,Plt,A1,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(1))*180/pi),'^o'])
    subplot(332)
    MAP3D_MILLER(Plg,Plt,A2,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(2))*180/pi),'^o'])
    subplot(333)
    MAP3D_MILLER(Plg,Plt,A3,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(3))*180/pi),'^0'])
    subplot(334)
    MAP3D_MILLER(Plg,Plt,A4,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(4))*180/pi),'^0'])
    subplot(335)
    MAP3D_MILLER(Plg,Plt,A5,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(5))*180/pi),'^o'])
    subplot(336)
    MAP3D_MILLER(Plg,Plt,A6,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(6))*180/pi),'^o'])
    subplot(337)
    MAP3D_MILLER(Plg,Plt,A7,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(7))*180/pi),'^o'])
    subplot(338)
    MAP3D_MILLER(Plg,Plt,A8,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(8))*180/pi),'^o'])
    subplot(339)
    MAP3D_MILLER(Plg,Plt,A9,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(9))*180/pi),'^o'])
    c = colorbar;
    if bool_criterion_find && ~bool_color_num
        c = colorbar('Ticks',array_bar,'TickLabels',label_bar);
    end
    set(c,'position',[0.9 0.08 0.025 0.9]); % x pos, y pos, largeur, ?
    ylabel(c,sprintf(name),'FontSize',15);
    suptitle(sprintf(['Miller representation of \n' name ' over half its surface']));

end

%%
if bool_ortho

    figure % 20
    subplot(331)
    MAP3D_ORTHO(Plg,Plt,A1,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(1))*180/pi),'^o'])
    subplot(332)
    MAP3D_ORTHO(Plg,Plt,A2,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(2))*180/pi),'^o'])
    subplot(333)
    MAP3D_ORTHO(Plg,Plt,A3,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(3))*180/pi),'^o'])
    subplot(334)
    MAP3D_ORTHO(Plg,Plt,A4,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(4))*180/pi),'^o'])
    subplot(335)
    MAP3D_ORTHO(Plg,Plt,A5,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(5))*180/pi),'^o'])
    subplot(336)
    MAP3D_ORTHO(Plg,Plt,A6,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(6))*180/pi),'^o'])
    subplot(337)
    MAP3D_ORTHO(Plg,Plt,A7,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(7))*180/pi),'^o'])
    subplot(338)
    MAP3D_ORTHO(Plg,Plt,A8,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(8))*180/pi),'^o'])
    subplot(339)
    MAP3D_ORTHO(Plg,Plt,A9,mi,ma)
    text(0,0,['\lambda=  ',num2str(lon(lm(9))*180/pi),'^0'])
    c = colorbar;
    if bool_criterion_find && ~bool_color_num
        c = colorbar('Ticks',array_bar,'TickLabels',label_bar);
    end
    set(c,'position',[0.9 0.08 0.025 0.9]); % x pos, y pos, largeur, ?
    ylabel(c,sprintf(name),'FontSize',15);
    suptitle(sprintf(['Orthographic representation of \n' name ' over half its surface']));

end

%%
if bool_ortho_3D

    h = figure; % 21
    MAP3D_ORTHO(Plg,Plt,A1,mi,ma);
    text(-0.2,0,'\lambda= ')
    text(0.0,0,num2str(lon(lm(1))*180/pi))
    text(0.1,0,'  ^o')
    c = colorbar;
    if bool_criterion_find && ~bool_color_num
        c = colorbar('Ticks',array_bar,'TickLabels',label_bar);
    end
    set(c,'position',[0.9 0.08 0.025 0.9]); % x pos, y pos, largeur, ?
    ylabel(c,sprintf(name),'FontSize',15);
    suptitle(sprintf(['Orthographic representation of \n' name ' ']));

    A_SavePlot(bool_save,h,sprintf(['Orthographic representation of \n' name ' ' model_name]));

end

%%
if bool_stereo

    figure; % 22
    subplot(3,3,1);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A1,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(1))*180/pi),'^o']);
    subplot(3,3,2);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A2,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(2))*180/pi),'^o']);
    subplot(3,3,3);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A3,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(3))*180/pi),'^o']);
    subplot(3,3,4);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A4,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(4))*180/pi),'^o']);
    subplot(3,3,5);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A5,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(5))*180/pi),'^o']);
    subplot(3,3,6);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A6,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(6))*180/pi),'^0']);
    subplot(3,3,7);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A7,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(7))*180/pi),'^o']);
    subplot(3,3,8);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A8,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(8))*180/pi),'^o']);
    subplot(3,3,9);
    MAP3D_STEREOGRAPHIC(Plg,Plt,A9,mi,ma);
    text(0,0,['\lambda=  ',num2str(lon(lm(9))*180/pi),'^o']);
    c = colorbar;
    if bool_criterion_find && ~bool_color_num
        c = colorbar('Ticks',array_bar,'TickLabels',label_bar);
    end
    set(c,'position',[0.9 0.08 0.025 0.9]); % x pos, y pos, largeur, ?
    ylabel(c,sprintf(name),'FontSize',15);
    suptitle(sprintf(['Stereographic representation of \n' name ' over half its surface']));

end


%set(get(h,'title'),'string','Potential');
