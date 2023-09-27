%function [disp,lat,lon,time]=PLOTS(A1,per,kpot,Vpot)
% This module uses rp (actual radius) in the spherical harmonics instead
% of ref(1) (r_ref)

%% Constants for stress calculations
%lambda_s=data(:,3);% lambda coeff de Lam� de la couche
%mu_s=data(:,4);%

% Size : nbcouches
mu_s0=real(mu*ref(3));%/ref(3);
lambda_s0=lambda*ref(3);%/ref(3);
%A0=(lambda_s0+2*mu_s0);
rr=ref(1); 

% Size : nbcouches*2+1 (central point and medians added)
y1=real(yi(1,:));%/ref(1);
y2=real(yi(2,:));%/ref(3);%/(ref(1)*ref(3)*sqrt(ref(4)/ref(1)));
y3=real(yi(3,:));%/ref(1);
y4=real(yi(4,:));%/ref(3);%/(ref(1)*ref(3)*sqrt(ref(4)/ref(1)));
y5=real(yi(5,:));%/(ref(1)*ref(4));%/(ref(1)*sqrt(ref(4)/ref(1)));
y6=real(yi(6,:));%/ref(4);%/ref(4);

if bool_greff_lefftz %y3 and y7 almost equal
    y7=real(yi(7,:));
end

% Size harmonisation
% r_s will be rewritten
[r_s,rho_s,lambda_s,mu_s,V_s]=TrancheV0(r_data,rho_data,lambda_data,mu_data,V_data,2*nbtranche,bool_interp_tranche);
rho_s = [rho_s(1) rho_s];
%lambda_s = [lambda_s(1) lambda_s]; 
%mu_s = [mu_s(1) mu_s]; % This mu_s can create discontinuities for sigma_tt and sigma_pp, only on Europa_Marusiak2021 (probably due to ocean layer)
%V_s = [V_s(1) V_s];


lambda_s=[lambda_s0(1) reshape([lambda_s0;lambda_s0],1,nbcouches*2)];
mu_s=[mu_s0(1) reshape([mu_s0;mu_s0],1,nbcouches*2)];

%A_s = (lambda_s+2*mu_s);
A0=(lambda_s0+2*mu_s0);
A_s=[A0(1) reshape([A0;A0],1,nbcouches*2)];

r_s = rp*ref(1);
r_s(1) = r_s(2)/10;

dy1=1./(A_s).*( y2-lambda_s./r_s.*(2*y1-n_harmo*(n_harmo+1)*y3)); %y3 and y7 almost equal
dy3=-y1./r_s+y3./r_s+y4./mu_s;

%% Stress and Strain initialisation (Order t p r)
Nradius = nbcouches*2+1;
Ntimeloop = Nt; % 1 Or equal to Nt

if ~bool_direct_stress && ~bool_stress_from_strain
    error('Please select a stress method calculation');
end

% Matrix initialisation - not defined because sparse thus faster

% Time = zeros(1,Ntimeloop);
% 
% urr = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% utt = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% upp = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% 
% sigma_rr = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_tt = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_pp = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_tp = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_pr = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_rt = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% 
% sigma_rr2 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_tt2 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_pp2 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_tp2 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_pr2 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_rt2 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% 
% sigma_xx = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_yy = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_zz = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_xy = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_yz = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_xz = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% 
% srr = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% stt = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% spp = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% stp = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% spr = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% srt = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% 
% e_rr = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% e_tt = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% e_pp = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% e_tp = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% e_pr = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% e_rt = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% 
% s_e = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% s_e21 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% s_m21 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% s_e22 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% s_m22 = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% sigma_e = zeros(Ntimeloop,Nlon,Nlat,Nradius);
% 
% phi = zeros(Ntimeloop,Nlon,Nlat,Nradius);

%% Loops
% Time loop
t=0;

simu_time = cputime;

for tt=1:1:Ntimeloop+1
    Time(tt)=t;
    cos_nt=cos(nn*t);%pulsation
    sin_nt=sin(nn*t);
    
    
    A20=-3*sqrt(pi/5)*cos_nt*epot;
    A22=epot*1/2*sqrt(6*pi/5)*(3*cos_nt+4/1i*sin_nt);
    A22n=epot*1/2*sqrt(6*pi/5)*(3*cos_nt-4/1i*sin_nt);
    
    for radius=1:1:Nradius
        for ix=1:1:Nlon
            for j=1:Nlat
                sinlat=sin(lat(j));
                coslat=cos(lat(j));
                cotlat=cot(lat(j));
                sinlon=sin(lon(ix));
                coslon=cos(lon(ix));
                
                %HARMONICS
                
                
                SY=(A20*Y20(ix,j)+A22*Y22(ix,j)+A22n*Y22n(ix,j));
                SdYdt=(A20*dY20dt(ix,j)+A22*dY22dt(ix,j)+A22n*dY22ndt(ix,j));
                SdYdp=(A20*dY20dp(ix,j)+A22*dY22dp(ix,j)+A22n*dY22ndp(ix,j));
                Sd2Ydt=(A20*d2Y20dt(ix,j)+A22*d2Y22dt(ix,j)+A22n*d2Y22ndt(ix,j));
                Sd2Ydp=(A20*d2Y20dp(ix,j)+A22*d2Y22dp(ix,j)+A22n*d2Y22ndp(ix,j));
                Sd2Ydpt=(A20*d2Y20dpt(ix,j)+A22*d2Y22dpt(ix,j)+A22n*d2Y22ndpt(ix,j));
                
                
                %T=(Sd2Ydt+cotlat*SdYdt+1/sinlat/sinlat*Sd2Ydp);
                %displacement
%                 urr(tt,ix,j,radius)=y1(radius)*SY;
%                 utt(tt,ix,j,radius)=y3(radius)*SdYdt;
%                 upp(tt,ix,j,radius)=y3(radius)/sinlat*SdYdp;

                if bool_displacement
                    urr(tt,ix,j,radius)=(yi(1,radius)*(A20*Y20(ix,j)+A22*Y22(ix,j)+A22n*Y22n(ix,j)));
                    %ur_21(tt,i,j)=real(yi(1,end))*(A20*Y20(i,j)+A21*Y21(i,j)+A22*Y22(i,j)+A22n*Y22n(i,j));
                    utt(tt,ix,j,radius)=(yi(3,radius)*(A20*dY20dt(ix,j)+A22*dY22dt(ix,j)+A22n*dY22ndt(ix,j)));
                    upp(tt,ix,j,radius)=(yi(3,radius)/sinlat*(A20*dY20dp(ix,j)+A22*dY22dp(ix,j)+A22n*dY22ndp(ix,j)));
                end
                
                
                if bool_direct_stress
                    % STRESSES
                    sigma_rr(tt,ix,j,radius)=y2(radius)*SY;
                    sigma_tt(tt,ix,j,radius)=(lambda_s(radius)*dy1(radius)+A_s(radius)/r_s(radius)*(2*y1(radius)-6*y3(radius))-2*mu_s(radius)/r_s(radius)*y1(radius))*SY - 2*mu_s(radius)/r_s(radius)*y3(radius)*(cotlat*SdYdt+1/sinlat^2*Sd2Ydp);
                    %sigma_tt1(tt,i,j,radius)=(lambda_s(radius)*dy1(radius)+A(radius)/r_s(radius)*(2*y1(radius)-6*y3(radius))-2*mu_s(radius)/r_s(radius)*y1(radius))*SY;
                    %sigma_tt2(tt,i,j,radius)=- 2*mu_s(radius)/r_s(radius)*y3(radius)*(coslat/sinlat*SdYdt+1/sinlat^2*Sd2Ydp);
                    sigma_pp(tt,ix,j,radius)=(lambda_s(radius)*dy1(radius)+A_s(radius)/r_s(radius)*(2*y1(radius)-6*y3(radius))-2*mu_s(radius)/r_s(radius)*y1(radius))*SY - 2*mu_s(radius)/r_s(radius)*y3(radius)*Sd2Ydt;
                    %sigma_pp1(tt,i,j,radius)=(lambda_s(radius)*dy1(radius)+A(radius)/r_s(radius)*(2*y1(radius)-6*y3(radius))-2*mu_s(radius)/r_s(radius)*y1(radius))*SY;
                    %sigma_pp2(tt,i,j,radius)= - 2*mu_s(radius)/r_s(radius)*y3(radius)*Sd2Ydt;
                    sigma_tp(tt,ix,j,radius)=2*mu_s(radius)/r_s(radius)*y3(radius)*(Sd2Ydpt/sinlat-coslat/sinlat^2*SdYdp);
                    sigma_pr(tt,ix,j,radius)=y4(radius)*SdYdp/sinlat;
                    sigma_rt(tt,ix,j,radius)=y4(radius)*SdYdt;
                end
                
                
                if bool_strain
                    if bool_greff_lefftz
                        % STRAIN
                        e_rr(tt,ix,j,radius)=dy1(radius)*SY;
                        e_tt(tt,ix,j,radius)=y1(radius)*SY/r_s(radius)+y3(radius)/r_s(radius)*Sd2Ydt;
                        e_pp(tt,ix,j,radius)=y7(radius)/r_s(radius)/sinlat/sinlat*Sd2Ydp+y7(radius)*cotlat/r_s(radius)*SdYdt+y1(radius)*SY/r_s(radius); % diff
                        %e_tt(tt,ix,j,radius)=e_pp(tt,ix,j,radius);
                        e_tp(tt,ix,j,radius)=y3(radius)/r_s(radius)/sinlat*Sd2Ydpt-y7(radius)/r_s(radius)/sinlat*cotlat*SdYdp;
                        
                        if mu_s(radius)~= 0
                            e_rt(tt,ix,j,radius)=y4(radius)/(2*mu_s(radius))*SdYdt;
                            e_pr(tt,ix,j,radius)=y4(radius)/(2*mu_s(radius)*sinlat)*SdYdp;
                        else
                            e_pr(tt,ix,j,radius)=0;
                            e_rt(tt,ix,j,radius)=0;
                        end
                    else
                        % STRAIN
                        e_rr(tt,ix,j,radius)=dy1(radius)*SY;
                        e_tt(tt,ix,j,radius)=y1(radius)*SY/r_s(radius)+y3(radius)/r_s(radius)*Sd2Ydt;
                        e_pp(tt,ix,j,radius)=y3(radius)/r_s(radius)/sinlat/sinlat*Sd2Ydp+y3(radius)*cotlat/r_s(radius)*SdYdt+y1(radius)*SY/r_s(radius); % diff
                        %e_tt(tt,ix,j,radius)=e_pp(tt,ix,j,radius);

                        %disp('Change in strain formula');
%                         e_tp(tt,ix,j,radius)=(1/(2*r_s(radius)))*((1/sinlat)*(y3(radius)*Sd2Ydpt)+y3(radius)/sinlat*Sd2Ydpt-(y3(radius)/sinlat)*cotlat);
% 
%                         if mu_s(radius)~= 0
%                             e_rt(tt,ix,j,radius)=(1/2)*((1/sinlat)*(y1(radius)*SdYdt)+dy3(radius)-y3(radius)/r_s(radius));
%                             e_pr(tt,ix,j,radius)=(1/2)*((1/(r_s(radius)*sinlat))*y1(radius)*SdYdp + dy3(radius)/sinlat - y3(radius)/sinlat/r_s(radius));
%                         else
%                             e_pr(tt,ix,j,radius)=0;
%                             e_rt(tt,ix,j,radius)=0;
%                         end

                        e_tp(tt,ix,j,radius)=y3(radius)/r_s(radius)/sinlat*(Sd2Ydpt-cotlat*SdYdp);
                        
                        if mu_s(radius)~= 0
                            e_rt(tt,ix,j,radius)=y4(radius)/(2*mu_s(radius))*SdYdt;
                            e_pr(tt,ix,j,radius)=y4(radius)/(2*mu_s(radius)*sinlat)*SdYdp;
                        else
                            e_pr(tt,ix,j,radius)=0;
                            e_rt(tt,ix,j,radius)=0;
                        end


                    end
                end

                if bool_stress_from_strain && bool_strain
                    % STRESS FROM STRAIN
                    sigma_rr(tt,ix,j,radius)=2*mu_s(radius)*e_rr(tt,ix,j,radius)+lambda_s(radius)*(e_rr(tt,ix,j,radius)+e_tt(tt,ix,j,radius)+e_pp(tt,ix,j,radius));
                    sigma_tt(tt,ix,j,radius)=2*mu_s(radius)*e_tt(tt,ix,j,radius)+lambda_s(radius)*(e_rr(tt,ix,j,radius)+e_tt(tt,ix,j,radius)+e_pp(tt,ix,j,radius));
                    sigma_pp(tt,ix,j,radius)=2*mu_s(radius)*e_pp(tt,ix,j,radius)+lambda_s(radius)*(e_rr(tt,ix,j,radius)+e_tt(tt,ix,j,radius)+e_pp(tt,ix,j,radius));
                    sigma_tp(tt,ix,j,radius)=2*mu_s(radius)*e_tp(tt,ix,j,radius);
                    sigma_pr(tt,ix,j,radius)=2*mu_s(radius)*e_pr(tt,ix,j,radius);
                    sigma_rt(tt,ix,j,radius)=2*mu_s(radius)*e_rt(tt,ix,j,radius);
                elseif bool_stress_from_strain && ~bool_strain
                    error('If stress is calculated from strain, then bool_strain must be true !');
                end
                
%                 % Test, to match exact figure from Minshull and Goulty
%                 % Not needed, as we use spherical coordinates (r,t,p)
%                 like Alterman whereas M&G use spherical polar coordinates
%                 (r,p,t)
%                 a = sigma_tt(tt,ix,j,radius);
%                 b = sigma_pp(tt,ix,j,radius);
%                 c = sigma_tp(tt,ix,j,radius);
%                 d = sigma_pr(tt,ix,j,radius);
%                 e = sigma_rt(tt,ix,j,radius);
%                 
%                 sigma_tt(tt,ix,j,radius)=b;
%                 sigma_pp(tt,ix,j,radius)=a;
%                 sigma_tp(tt,ix,j,radius)=c;
%                 sigma_pr(tt,ix,j,radius)=e;
%                 sigma_rt(tt,ix,j,radius)=d;

                %% Converting stress signs so that compression is positive for Mohr Coulomb failure criterion
                sigma_rr(tt,ix,j,radius) = -sigma_rr(tt,ix,j,radius);
                sigma_tt(tt,ix,j,radius) = -sigma_tt(tt,ix,j,radius);
                sigma_pp(tt,ix,j,radius) = -sigma_pp(tt,ix,j,radius);
                sigma_rt(tt,ix,j,radius) = -sigma_rt(tt,ix,j,radius);
                sigma_tp(tt,ix,j,radius) = -sigma_tp(tt,ix,j,radius);
                sigma_pr(tt,ix,j,radius) = -sigma_pr(tt,ix,j,radius);

                
                if bool_stress_cartesian
                    % STRESS CONVERSION - From spherical coordinates to
                    % cartesian
                    sigma_xx(tt,ix,j,radius)=coslat*coslon*sigma_rr(tt,ix,j,radius) - sinlat*coslon*sigma_tt(tt,ix,j,radius) + sinlon*sigma_pp(tt,ix,j,radius);
                    sigma_yy(tt,ix,j,radius)=coslat*sinlon*sigma_rr(tt,ix,j,radius) - sinlat*sinlon*sigma_tt(tt,ix,j,radius) - coslon*sigma_pp(tt,ix,j,radius);
                    sigma_zz(tt,ix,j,radius)=sinlat*sigma_rr(tt,ix,j,radius) + coslat*sinlon*sigma_tt(tt,ix,j,radius);
                    sigma_xy(tt,ix,j,radius)=1; %?
                    sigma_yz(tt,ix,j,radius)=1; %?
                    sigma_xz(tt,ix,j,radius)=1; %?
                end
                
%                 srr(tt,ix,j,radius)=(sigma_rr(tt,ix,j,radius));
%                 stt(tt,ix,j,radius)=(sigma_tt(tt,ix,j,radius));
%                 spp(tt,ix,j,radius)=(sigma_pp(tt,ix,j,radius));
%                 stp(tt,ix,j,radius)=(sigma_tp(tt,ix,j,radius));
%                 spr(tt,ix,j,radius)=(sigma_pr(tt,ix,j,radius));
%                 srt(tt,ix,j,radius)=(sigma_rt(tt,ix,j,radius));
% 
%                 
%                 if bool_deviatoric_stress
%                     % EFFECTIVE STRESS - Stress deviator tensor calculation
%                     sm=1/3*(srr(tt,ix,j,radius)+stt(tt,ix,j,radius)+spp(tt,ix,j,radius));
% 
%                     srr(tt,ix,j,radius)=srr(tt,ix,j,radius)-sm;
%                     stt(tt,ix,j,radius)=stt(tt,ix,j,radius)-sm;
%                     spp(tt,ix,j,radius)=spp(tt,ix,j,radius)-sm;
%                 end
                
                % No need to add pressure yet, as it is a diagonal tensor
                % it does not change diagonalisation
                M_sigma = [sigma_rr(tt,ix,j,radius) sigma_rt(tt,ix,j,radius) sigma_pr(tt,ix,j,radius);...
                    sigma_rt(tt,ix,j,radius) sigma_tt(tt,ix,j,radius) sigma_tp(tt,ix,j,radius);...
                    sigma_pr(tt,ix,j,radius) sigma_tp(tt,ix,j,radius) sigma_pp(tt,ix,j,radius)];
                
                if bool_eigen
                    
                    [eig_vect,eig_val] = eig(M_sigma); % Slow...
                    
                    % Orthogonality verification
                    prod_scal1 = dot(eig_vect(:,1),eig_vect(:,2));
                    prod_scal2 = dot(eig_vect(:,2),eig_vect(:,3));
                    prod_scal3 = dot(eig_vect(:,3),eig_vect(:,1));
                    
                    if max([abs(prod_scal1) abs(prod_scal2) abs(prod_scal3)]) > tolerance
                        fprintf(['\nThere is a non orthogonal diagonal matrix here at time '... 
                            num2str(tt) 's, lat ' num2str(ix*180/pi) '�, lon ' num2str(j*180/pi)...
                            '� and radius ' num2str(r_s(radius)) 'm\nMax value of scalar product is '...
                            num2str(max([abs(prod_scal1) abs(prod_scal2) abs(prod_scal3)])) '\n']);
                    end
                    
                    sigma_val = diag(eig_val); % eig return a diag matrix of eigen values
                    val_p = sort(sigma_val); % Ascending order
                    s3(tt,ix,j,radius) = val_p(1);
                    s2(tt,ix,j,radius) = val_p(2);
                    s1(tt,ix,j,radius) = val_p(3);

                    P{tt,ix,j,radius} = eig_vect; % diagonalization matrix such as P-1 M_sigma P = eig_val

%                     val_p = eig(M_sigma);
%                     smin = val_p(1);
%                     smax = val_p(3);

                    % s_0=(srr+stt+spp)/3; % Toujours nul
                    % s_e(tt,ix,j,radius)=sqrt(1/2*((srr-s_0)^2+(stt-s_0)^2+(spp-s_0)^2)+(stp^2+spr^2+srt^2));
                    % s_e21(tt,ix,j,radius)=sqrt(1/2*(stt-srr)^2); % von Mises maximum shear stress (equatorial plane) ?
                    
                else
                    
                    s_simple = sort([srr(tt,ix,j,radius);stt(tt,ix,j,radius);spp(tt,ix,j,radius)]);
                    s3(tt,ix,j,radius) = s_simple(1);
                    s2(tt,ix,j,radius) = s_simple(2);
                    s1(tt,ix,j,radius) = s_simple(3);

                end
                
                if bool_stress_criteria
                    s_e(tt,ix,j,radius)=sqrt(1/2*((s1(tt,ix,j,radius)-s2(tt,ix,j,radius))^2+(s1(tt,ix,j,radius)-s3(tt,ix,j,radius))^2+(s2(tt,ix,j,radius)-s3(tt,ix,j,radius))^2)); % von Mises stress
                    s_e21(tt,ix,j,radius)=sqrt(1/2*(s3(tt,ix,j,radius)-s1(tt,ix,j,radius))^2); % Max shear stress
                    s_m21(tt,ix,j,radius)=sqrt(1/2*(s3(tt,ix,j,radius)+s1(tt,ix,j,radius))^2); % Max normal stress
                    
                    %STRESS IN RADIAL PLANE ; mb delete this ?
                    if sqrt(sigma_tp(tt,ix,j,radius)^2) ==0

                        s_radial_1= stt(tt,ix,j,radius);
                        s_radial_2= spp(tt,ix,j,radius);
                        s_radial_3= 0;
                    else

                        s_radial_1=(stt(tt,ix,j,radius)+spp(tt,ix,j,radius)+sqrt((stt(tt,ix,j,radius)-spp(tt,ix,j,radius))^2+4*stp(tt,ix,j,radius)^2))/2;
                        s_radial_2=(stt(tt,ix,j,radius)+spp(tt,ix,j,radius)-sqrt((stt(tt,ix,j,radius)-spp(tt,ix,j,radius))^2+4*stp(tt,ix,j,radius)^2))/2;
                        s_radial_3=0;

                    end

                    sigma_e(tt,ix,j,radius)=sqrt(1/2*((s_radial_1-s_radial_2)^2+(s_radial_2-s_radial_3)^2+(s_radial_3-s_radial_1)^2)); % von Mises stress in the radial plane
                    s_e22(tt,ix,j,radius)=1/2*(s_radial_1-s_radial_2); % von Mises maximum shear stress (radial plane)
                    s_m22(tt,ix,j,radius)=1/2*(s_radial_1+s_radial_2); % von Mises maximum normal stress (radial plane)
                    
                end
                
            end
        end
    end
    t=t+Tf/Ntimeloop;
end

% simu_time = cputime - simu_time;
% fprintf(['\nElapsed time : ' num2str(simu_time) 's\n']);

Stress_Plots;