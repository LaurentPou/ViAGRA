%close all;


if bool_plot_disp
    h=figure;
    subplot(1,3,1);
    plot(squeeze(urr(time_plot,lon_plot,lat_plot,:))*disp_factor,r_s*depth_factor);
    xlabel(sprintf(['u_{rr} in ' disp_unit]));
    subplot(1,3,2);
    plot(squeeze(utt(time_plot,lon_plot,lat_plot,:))*disp_factor,r_s*depth_factor);
    xlabel(sprintf(['u_{\\theta\\theta} in ' disp_unit]));
    subplot(1,3,3);
    plot(squeeze(upp(time_plot,lon_plot,lat_plot,:))*disp_factor,r_s*depth_factor);
    xlabel(sprintf(['u_{\\phi\\phi} in ' disp_unit]));
    
    A_SavePlot(bool_save,h,sprintf(['Displacements ' model_name]));
    
end

if bool_plot_strain
    h=figure;
    subplot(2,3,1);
    plot(squeeze(e_rr(time_plot,lon_plot,lat_plot,:)),r_s*depth_factor);
    xlabel('e rr');
    subplot(2,3,2);
    plot(squeeze(e_tt(time_plot,lon_plot,lat_plot,:)),r_s*depth_factor);
    xlabel('e tt');
    subplot(2,3,3);
    plot(squeeze(e_pp(time_plot,lon_plot,lat_plot,:)),r_s*depth_factor);
    xlabel('e pp');
    
    subplot(2,3,4);
    plot(squeeze(e_tp(time_plot,lon_plot,lat_plot,:)),r_s*depth_factor);
    xlabel('e tp');
    subplot(2,3,5);
    plot(squeeze(e_pr(time_plot,lon_plot,lat_plot,:)),r_s*depth_factor);
    xlabel('e pr');
    subplot(2,3,6);
    plot(squeeze(e_rt(time_plot,lon_plot,lat_plot,:)),r_s*depth_factor);
    xlabel('e rt');
    
    A_SavePlot(bool_save,h,sprintf(['Strains ' model_name]));
end


if bool_plot_stress_spherical
    
    h=figure;
    subplot(2,3,1);
    plot(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:))*stress_r_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{rr} in ' stress_r_unit]));
    subplot(2,3,2);
    plot(squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{\\theta\\theta} in ' stress_r_unit]));
    subplot(2,3,3);
    plot(squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{\\phi\\phi} in ' stress_r_unit]));
    
    subplot(2,3,4);
    plot(squeeze(sigma_tp(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{\\theta\\phi} in ' stress_r_unit]));
    subplot(2,3,5);
    plot(squeeze(sigma_pr(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{\\phi r} in ' stress_r_unit]));
    subplot(2,3,6);
    plot(squeeze(sigma_rt(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{r\\theta} in ' stress_r_unit]));
    
    A_SavePlot(bool_save,h,sprintf(['Stress ' model_name]));
    
%     h=figure;
%     plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:)))*stress_factor,'r');
%     hold on;
%     plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:)))*stress_factor,'g');
%     plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:)))*stress_factor,'b');
%     %plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:)))*stress_factor,'magenta');
%     plot((ref(1)-r_s).*depth_factor,zeros(1,numel(r_s)),'black');
%     hold off;
%     xlabel(sprintf(['Depth (' depth_unit ')']));
%     legend(sprintf(['\\sigma_{rr} in ' stress_r_unit]),...
%         sprintf(['\\sigma_{\\theta\\theta} in ' stress_r_unit]),...
%         sprintf(['\\sigma_{\\phi\\phi} in ' stress_r_unit]));

nb_rs = length(r_s);
array_rs = 1:round(nb_rs/4):nb_rs;

    h=figure;
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:)))*stress_factor,'-^','Color','r','MarkerIndices',array_rs);
    hold on;
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:)))*stress_factor,'-s','Color','g','MarkerIndices',array_rs);
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:)))*stress_factor,'-o','Color','b','MarkerIndices',array_rs);
    plot((ref(1)-r_s).*depth_factor,zeros(1,numel(r_s)),'black');
    hold off;
    xlabel(sprintf(['Depth (' depth_unit ')']));
    legend(sprintf(['\\sigma_{rr} in ' stress_r_unit]),...
        sprintf(['\\sigma_{\\theta\\theta} in ' stress_r_unit]),...
        sprintf(['\\sigma_{\\phi\\phi} in ' stress_r_unit]));
    title(sprintf(['Timestep ' num2str(time_plot) ' out of ' num2str(Nt)]))
    
    A_SavePlot(bool_save,h,sprintf(['All stress vs Depth at Location ' model_name]));

%     h=figure;
%     plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rt(time_plot,lon_plot,lat_plot,:)))*stress_factor,'r');
%     hold on;
%     plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_tp(time_plot,lon_plot,lat_plot,:)))*stress_factor,'g');
%     plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_pr(time_plot,lon_plot,lat_plot,:)))*stress_factor,'b');
%     %plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:)))*stress_factor,'magenta');
%     plot((ref(1)-r_s).*depth_factor,zeros(1,numel(r_s)),'black');
%     hold off;
%     xlabel(sprintf(['Depth (' depth_unit ')']));
%     legend(sprintf(['\\sigma_{r\\theta} in ' stress_r_unit]),...
%         sprintf(['\\sigma_{\\theta\\phi} in ' stress_r_unit]),...
%         sprintf(['\\sigma_{\\phi r} in ' stress_r_unit]));
%     

    h=figure;
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_rt(time_plot,lon_plot,lat_plot,:)))*stress_factor,'-^','Color','r','MarkerIndices',array_rs);
    hold on;
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_tp(time_plot,lon_plot,lat_plot,:)))*stress_factor,'-o','Color','g','MarkerIndices',array_rs);
    plot((ref(1)-r_s).*depth_factor,fliplr(squeeze(sigma_pr(time_plot,lon_plot,lat_plot,:)))*stress_factor,'-s','Color','b','MarkerIndices',array_rs);
    plot((ref(1)-r_s).*depth_factor,zeros(1,numel(r_s)),'black');
    hold off;
    xlabel(sprintf(['Depth (' depth_unit ')']));
    legend(sprintf(['\\sigma_{r\\theta} in ' stress_r_unit]),...
        sprintf(['\\sigma_{\\theta\\phi} in ' stress_r_unit]),...
        sprintf(['\\sigma_{\\phi r} in ' stress_r_unit]));
    title(sprintf(['Timestep ' num2str(time_plot) ' out of ' num2str(Nt)]))

    h=figure;
    plot(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['Trace of regular stress tensor in ' stress_unit]));
    
    
end

if bool_plot_stress_cartesian && bool_stress_cartesian
    
    h=figure;
    subplot(2,3,1);
    plot(squeeze(sigma_xx(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{xx} in ' stress_unit]));
    subplot(2,3,2);
    plot(squeeze(sigma_yy(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{yy} in ' stress_unit]));
    subplot(2,3,3);
    plot(squeeze(sigma_zz(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{zz} in ' stress_unit]));
    
    subplot(2,3,4);
    plot(squeeze(sigma_xy(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{xy} in ' stress_unit]));
    subplot(2,3,5);
    plot(squeeze(sigma_yz(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{yz} in ' stress_unit]));
    subplot(2,3,6);
    plot(squeeze(sigma_xz(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_{xz} in ' stress_unit]));
    
    A_SavePlot(bool_save,h,sprintf(['Cartesian stress ' model_name]));
    
    h=figure;
    plot(squeeze(sigma_xx(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_yy(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_zz(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['Trace of regular stress tensor in ' stress_unit]));
    
end
    




if bool_plot_stress_criteria && bool_stress_criteria
   
    h=figure;
    subplot(3,2,1);
    plot(squeeze(s_e(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_ep in ' stress_unit]));
    subplot(3,2,2);
    plot(squeeze(sigma_e(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['sigma_ep in ' stress_unit]));
    subplot(3,2,3);
    plot(squeeze(s_e21(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_e21p in ' stress_unit]));
    subplot(3,2,4);
    plot(squeeze(s_e22(1,1,1,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_e22p in ' stress_unit]));
    subplot(3,2,5);
    plot(squeeze(s_m21(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_m21p in ' stress_unit]));
    subplot(3,2,6);
    plot(squeeze(s_m22(1,1,1,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_m22p in ' stress_unit]));
    
    A_SavePlot(bool_save,h,sprintf(['Stress elastic criteria ' model_name]));
    
    h=figure;
    subplot(1,2,1);
    hold on;
    plot(squeeze(s_e21(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor,'r');
    plot(squeeze(s_m21(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor,'g');
    legend(sprintf(['s_e21p in ' stress_unit]),sprintf(['s_m21p in ' stress_unit]));
    hold off;
    subplot(1,2,2);
    hold on;
    plot(squeeze(s_e22(1,1,1,:))*stress_factor,r_s*depth_factor,'r');
    plot(squeeze(s_m22(1,1,1,:))*stress_factor,r_s*depth_factor,'g');
    legend(sprintf(['s_e22p in ' stress_unit]),sprintf(['s_m22p in ' stress_unit]));
    hold off;
    
    A_SavePlot(bool_save,h,sprintf(['Stress elastic criteria evaluation ' model_name]));
    
    disp('s_ep maximal at :');
    disp(ref(1)*depth_factor*(rp(find(squeeze(s_e(time_plot,lon_plot,lat_plot,:)) == max(squeeze(s_e(time_plot,lon_plot,lat_plot,:)))))));
    disp(sprintf([depth_unit ' radius']));
    disp('');
    
    disp('sigma_ep maximal at :');
    disp(ref(1)*depth_factor*(rp(find(squeeze(sigma_e(time_plot,lon_plot,lat_plot,:)) == max(squeeze(sigma_e(time_plot,lon_plot,lat_plot,:)))))));
    disp(sprintf([depth_unit ' radius']));
    disp('');
    
    disp('s_e21p maximal at :');
    disp(ref(1)*depth_factor*(rp(find(squeeze(s_e21(time_plot,lon_plot,lat_plot,:)) == max(squeeze(s_e21(time_plot,lon_plot,lat_plot,:)))))));
    disp(sprintf([depth_unit ' radius']));
    disp('');
    
    disp('s_e22p maximal at :');
    disp(ref(1)*depth_factor*(rp(find(squeeze(s_e22(time_plot,lon_plot,lat_plot,:)) == max(squeeze(s_e22(time_plot,lon_plot,lat_plot,:)))))));
    disp(sprintf([depth_unit ' radius']));
    disp('');
    
    %structure strainstress de sortie.
    strainstress.s_ep=abs(s_e(:,:,:,end));
    strainstress.sigma_ep=abs(sigma_e(:,:,:,end));
    strainstress.s_e21p=abs(s_e21(:,:,:,end));
    strainstress.s_e22p=abs(s_e22(:,:,:,end));
    
end

if bool_plot_deviatoric_stress
    
    h=figure;
    subplot(2,3,1);
    plot(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:))*stress_r_factor,r_s*depth_factor);
    xlabel(sprintf(['s_{rr} in ' stress_r_unit]));
    subplot(2,3,2);
    plot(squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_{\\theta\\theta} in ' stress_unit]));
    subplot(2,3,3);
    plot(squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_{\\phi\\phi} in ' stress_unit]));
    
    subplot(2,3,4);
    plot(squeeze(sigma_tp(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_{\\theta\\phi} in ' stress_unit]));
    subplot(2,3,5);
    plot(squeeze(sigma_pr(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_{\\phi r} in ' stress_unit]));
    subplot(2,3,6);
    plot(squeeze(sigma_rt(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['s_{r\\theta} in ' stress_unit]));
    
    A_SavePlot(bool_save,h,sprintf(['Deviatoric stress ' model_name]));
    
    h=figure;
    plot(squeeze(sigma_rr(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_tt(time_plot,lon_plot,lat_plot,:))+squeeze(sigma_pp(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['Trace of deviatoric stress tensor in ' stress_unit]));
    
end

if bool_plot_main_stress
    
    h=figure;
    subplot(1,3,1);
    plot(squeeze(s1(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_1 in ' stress_unit]));
    subplot(1,3,2);
    plot(squeeze(s2(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_2 in ' stress_unit]));
    subplot(1,3,3);
    plot(squeeze(s3(time_plot,lon_plot,lat_plot,:))*stress_factor,r_s*depth_factor);
    xlabel(sprintf(['\\sigma_3 in ' stress_unit]));
    
    A_SavePlot(bool_save,h,sprintf(['Principal stress ' model_name]));
    
end

 
%%
if bool_surface_stress_maps && bool_stress_criteria
 
%     maa=max([max(max(max(s_ep))),max(max(max(sigma_ep)))]);
%     mii=min([min(min(min(s_ep))),min(min(min(sigma_ep)))]);
% 
%     h=figure;
%     title('Effective stress','FontSize',14)
%     subplot(221)
%     Cizdir1(lon,lat,s_ep,Nt,time_plot,per,1,mii,maa);
%     title('(a)','FontSize',18)
%     h=colorbar('h');
%     subplot(222)
%     Cizdir1(lon,lat,sigma_ep,Nt,time_plot,1,1,mii,maa);
%     title('(b)','FontSize',18)
%     subplot(223)
%     Cizdir1(lon,lat,s_e21p,Nt,time_plot,per,1,mii,maa);
%     title('(b)','FontSize',18)
%     h=colorbar('h');
%     subplot(224)
%     Cizdir1(lon,lat,s_e22p,Nt,time_plot,per,1,mii,maa);
%     title('(b)','FontSize',18)
%     h=colorbar('h');

    maa=max([max(max(max(s_e(:,:,:,end))))*stress_factor,max(max(max(sigma_e(:,:,:,end))))*stress_factor]);
    mii=min([min(min(min(s_e(:,:,:,end)))*stress_factor),min(min(min(sigma_e(:,:,:,end))))*stress_factor]);
    h=figure;
    subplot(221)
    Cizdir1(lon,lat,s_e(:,:,:,end)*stress_factor,Nt,time_plot,per,1,mii,maa);
    h=colorbar('h');
    title(sprintf(['Effective stress in' stress_unit]),'FontSize',14)
    subplot(222)

    Y1=abs(sigma_rr(:,:,:,end)*stress_r_factor);
    maa=max([max(max(max(Y1)))]);
    mii=min([min(min(min(Y1)))]);

    Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
    title(sprintf(['\\sigma_{rr} in ' stress_r_unit]),'FontSize',18)
    h=colorbar('h');

    subplot(223)
    Y1=abs(sigma_pp(:,:,:,end)*stress_factor);
    maa=max([max(max(max(Y1)))]);
    mii=min([min(min(min(Y1)))]);

    Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
    title(sprintf(['\\sigma_{\\phi\\phi} in ' stress_unit]),'FontSize',18)
    h=colorbar('h');



    subplot(224)
    Y1=abs(sigma_tt(:,:,:,end)*stress_factor);
    maa=max([max(max(max(Y1)))]);
    mii=min([min(min(min(Y1)))]);

    Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
    title(sprintf(['\\sigma_{\\theta\\theta} in ' stress_unit]),'FontSize',18)
    h=colorbar('h');

    
%     h=figure;
% 
%     Y1=abs(sigma_rr(:,:,:,end)*stress_r_factor);
%     Y2=abs(sigma_rr(:,:,:,end)*stress_r_factor);
%     maa=max([max(max(max(Y1))),max(max(max(Y2)))]);
%     mii=min([min(min(min(Y1))),min(min(min(Y2)))]);
% 
%     subplot(321)
%     title(sprintf(['from \\sigma_{rr} in ' stress_r_unit]),'FontSize',14)
%     Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
%     title(sprintf(['from \\sigma_{rr} in ' stress_r_unit]),'FontSize',18)
% 
%     subplot(322)
%     Cizdir1(lon,lat,Y2,Nt,time_plot,per,1,mii,maa);
%     title(sprintf(['from strain_{rr} in ' stress_r_unit]),'FontSize',18)
% 
%     Y1=abs(sigma_pp(:,:,:,end));
%     Y2=abs(sigma_pp(:,:,:,end));
%     maa=max([max(max(max(Y1))),max(max(max(Y2)))]);
%     mii=min([min(min(min(Y1))),min(min(min(Y2)))]);
% 
%     subplot(323)
%     title('from sigma pp','FontSize',14)
%     Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
%     title('from stress pp','FontSize',18)
% 
%     subplot(324)
%     Cizdir1(lon,lat,Y2,Nt,time_plot,per,1,mii,maa);
%     title('from strain pp','FontSize',18)
% 
%     Y1=abs(sigma_tt(:,:,:,end));
%     Y2=abs(sigma_tt(:,:,:,end));
%     maa=max([max(max(max(Y1))),max(max(max(Y2)))]);
%     mii=min([min(min(min(Y1))),min(min(min(Y2)))]);
% 
%     subplot(325)
%     title('from sigma pp','FontSize',14)
%     Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
%     title('from stress tt','FontSize',18)
% 
%     subplot(326)
%     Cizdir1(lon,lat,Y2,Nt,time_plot,per,1,mii,maa);
%     title('from strain tt','FontSize',18)
% 
% 
% 
% 
%     h=figure;
% 
%     Y1=abs(sigma_tp(:,:,:,end));
%     Y2=abs(sigma_tp(:,:,:,end));
%     maa=max([max(max(max(Y1))),max(max(max(Y2)))]);
%     mii=min([min(min(min(Y1))),min(min(min(Y2)))]);
% 
%     subplot(321)
%     title('from sigma rr','FontSize',14)
%     Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
%     title('from stress rr','FontSize',18)
% 
%     subplot(322)
%     Cizdir1(lon,lat,Y2,Nt,time_plot,per,1,mii,maa);
%     title('from strain rr','FontSize',18)
% 
%     Y1=abs(sigma_pr(:,:,:,end));
%     Y2=abs(sigma_pr(:,:,:,end));
%     maa=max([max(max(max(Y1))),max(max(max(Y2)))]);
%     mii=min([min(min(min(Y1))),min(min(min(Y2)))]);
% 
%     subplot(323)
%     title('from sigma pp','FontSize',14)
%     Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
%     title('from stress pp','FontSize',18)
% 
%     subplot(324)
%     Cizdir1(lon,lat,Y2,Nt,time_plot,per,1,mii,maa);
%     title('from strain pp','FontSize',18)
% 
%     Y1=abs(sigma_rt(:,:,:,end));
%     Y2=abs(sigma_rt(:,:,:,end));
%     maa=max([max(max(max(Y1))),max(max(max(Y2)))]);
%     mii=min([min(min(min(Y1))),min(min(min(Y2)))]);
% 
%     subplot(325)
%     title('from sigma pp','FontSize',14)
%     Cizdir1(lon,lat,Y1,Nt,time_plot,per,1,mii,maa);
%     title('from stress tt','FontSize',18)
% 
%     subplot(326)
%     Cizdir1(lon,lat,Y2,Nt,time_plot,per,1,mii,maa);
%     title('from strain tt','FontSize',18)

end

%%

if bool_surface_stress_3D
    if bool_sweep_lon
    %     Cizdir9_New_LonSweep(lon,lat,real(ur(:,:,:,radius_plot))*disp_factor,Nt,time_plot,per,sprintf[('Radial displacement (' disp_unit ')')],bool_color_num,bool_save,model_name);
  
        Cizdir9_New_LonSweep(lon,lat,sigma_rr(:,:,:,radius_plot)*stress_r_factor,Nt,time_plot,per,sprintf(['\\\\sigma_{r} (' stress_r_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_LonSweep(lon,lat,sigma_tt(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['\\\\sigma_{\\\\theta} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_LonSweep(lon,lat,sigma_pp(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['\\\\sigma_{\\\\phi} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_LonSweep(lon,lat,sigma_rt(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['\\\\sigma_{r\\\\theta} (' stress_r_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_LonSweep(lon,lat,sigma_tp(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['\\\\sigma_{\\\\theta\\\\phi} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_LonSweep(lon,lat,sigma_pr(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['\\\\sigma_{\\\\phi r} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_LonSweep(lon,lat,s1(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['Principal max stress \\\\sigma_{1} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_LonSweep(lon,lat,s3(:,:,:,radius_plot)*stress_factor,Nt,time_plot,per,sprintf(['Principal min stress \\\\sigma_{3} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);
    
    elseif bool_sweep_time
        
        %     Cizdir9_New_TimeSweep(Time,lat,real(ur(:,:,:,radius_plot)),Nt,time_plot,per,'Radial displacement (m)',bool_color_num,model_name);
    
        Cizdir9_New_TimeSweep(Time,lat,sigma_rr(:,:,:,radius_plot)*stress_r_factor,Nt,lon_plot,per,sprintf(['Stress \\\\sigma_r (' stress_r_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_TimeSweep(Time,lat,sigma_tt(:,:,:,radius_plot)*stress_factor,Nt,lon_plot,per,sprintf(['\\\\sigma_{\\\\theta} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_TimeSweep(Time,lat,sigma_pp(:,:,:,radius_plot)*stress_factor,Nt,lon_plot,per,sprintf(['\\\\sigma_{\\\\phi} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_TimeSweep(Time,lat,s1(:,:,:,radius_plot)*stress_factor,Nt,lon_plot,per,sprintf(['Principal max stress \\\\sigma_{1} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);

        Cizdir9_New_TimeSweep(Time,lat,s3(:,:,:,radius_plot)*stress_factor,Nt,lon_plot,per,sprintf(['Principal min stress \\\\sigma_{3} (' stress_unit ')']),bool_color_num,bool_save,model_name,bool_maps_type);
        
        
    end
end