% Code calculating the displacement and stress on a satellite in
% synchronous rotation around its main body ; those outputs are due to the
% excentricity of the satellite orbit
% This function for cohesion and friction sweep at once

%% 1-Parametering:


%% Temporary parameters fixed
obl = 0;
G = 6.67259e-11;% G constante gravi

%% Rheology
% Choice of rheology : %maxwell (1) ; %Kelvin-Voight (2) ; %Standard liner Solid (SLS) (3)
%Bourger's Body (4) ; %Caputo (5) (non lineraire)
type=1;

%% Plots
size_font = 14;%22;
set(0,'DefaultTextFontSize',size_font);
set(0,'DefaultAxesFontSize',size_font);
bool_save = 0; %
bool_data_save = 0; % Only data for Failure as of now
bool_fus = 0;

disp_factor = 1;%10^9;
disp_unit = 'm';%'nm';
gravity_factor = 1;%10^12;
gravity_unit = 'm';%'pm';
gravity_ref_factor = 1;%10^6;
gravity_ref_unit = 'm';%'�m';
depth_factor = 10^(-3);
depth_unit = 'km';
stress_r_factor = 1;%10^(-5);%10^6;
stress_r_unit = 'Pa';%'bar';%'�Pa';
stress_factor = 1;%10^(-5);%10^3;
stress_unit = 'Pa';%'bar';%'mPa';

%% Boolean for plots and calculations
bool_surface_disp_maps = 1; % Surface displacement maps
bool_surface_stress_maps = 1; % Surface stress maps
bool_surface_Mohr_maps = 1;
bool_surface_stress_3D = 1; %
bool_plot_yi = 1; %

bool_ham = 1;
bool_miller = 1;
bool_ortho = 1;
bool_ortho_3D = 1;
bool_stereo = 1;
bool_maps_type = [bool_ham bool_miller bool_ortho bool_ortho_3D bool_stereo];

bool_displacement = 1; %
bool_plot_displacement = 1;
bool_plot_displacement_time = 1;
bool_plot_all_displacement = 1; % Plots ur but also ut and up

bool_strain = 1; %
bool_plot_disp = 1; % wrt to depth
bool_plot_strain = 1; %

bool_rotation_speed = 0;
bool_plot_rotation_speed = 0;

bool_plot_acc = 0;
bool_plot_tilt = 0;
bool_plot_gravity = 0; %

bool_direct_stress = 0; %
bool_plot_stress_spherical = 1;

if ~bool_direct_stress
    bool_strain = 1;
    bool_stress_from_strain = 1;
else
    bool_stress_from_strain = 0;
end

bool_stress_cartesian = 0;
bool_plot_stress_cartesian = 0;

bool_deviatoric_stress = 0;
bool_plot_deviatoric_stress = 0;

bool_stress_criteria = 0;
bool_plot_stress_criteria = 0;

bool_eigen = 1; % Diagonalisation (1) or Approximation (0)
bool_plot_main_stress = 1; %

bool_greff_lefftz = 0; % For Love numbers spheroidal
bool_disp_LN_GL = 1 && bool_greff_lefftz; % Display Greff Lefftz analytical Love numbers

if bool_greff_lefftz
    bool_strain = 1;
    bool_stress_from_strain = 1;
    bool_direct_stress = 0;
end

% Relative to failure criteria
bool_pressure_failure = 1; % Add pressure in criterion
bool_plot_failure = 1; %
bool_plot_failure_log = 0; %

bool_loaded_model = 0;
bool_final_model = 0; %

bool_seismic_model = 1;
bool_geodesic_model = 1;
bool_mechanical_model = 1;

% For Cizdir9 plots
bool_color_num = 1;
bool_sweep_lon = 1;
bool_sweep_time = 0;

c_displacement=1; % Surface Displacement Calculations ;
c_stress=1; % Stress Calculations
c_failure=1; % Failure Criterion Calculations
c_dis_int=0; % Internal Dissipation Calculations

bool_print_tidal_outputs = 1;

if bool_pressure_failure
    disp('Lithostatic pressure is in failure criterion');
    legend_pressure = 'with pressure in failure criterion';
else
    disp('Lithostatic pressure is neglected in failure criterion');
    legend_pressure = 'without pressure in failure criterion';
end

%% Constants
N_lat=12; % N_lat + 1 pts ; is in fact colatitude
N_lon=32; % N_lon + 1 pts
Frac_lat=100; % lat(1) = 1/fraclat*lat(2)to avoid 0 in denominator
Nt=16;

time_plot=5;%floor(Nt/4)+1;
lat_plot=1;%1;%round(N_lat*((90+13.2)/180)); % lat -13.2°
lon_plot=1;%round(N_lon*((360-31.1)/360)); % longitude -31.1°
radius_plot=0.9999;%0.8;% 0 for center, 1 for surface

% Inversion parameters
n_harmo=2; % degre d'harmonique calcul?s% ATTENTION 2 est la valeur minimum, le degre 1 correspondant au deplacement du CM est ignore
n_poly=8; % degre d'interpolation polynomiale
%nbtranche=10;%nb de tranche par couche (DidymoonC : 1000, 500 otherwise) , 10 Moon (20 pour Bills), max 50 for Europa (to fit with Tobie et al 2005, otherwise numerical errors appear)
appli_pression=0; % oui si == 1
appli_potentiel=~appli_pression; % oui si ==1
transfo_laplace=appli_pression; % Laplace <-> pression, Fourier <-> potentiel
bool_interp_tranche = 1; % Linear interpolation between interfaces with different physical parameters or not
% TBD: Currently must be activated

flag=0;%initialisation du flag de conditionnement matriciel

% Limite valeur nulle
tolerance = 1e-10;

% Model name
model_name = sprintf([nom_file '_' num2str(nbtranche)]);


%% 2a-Initialisation
if numel(friction) > 1
    if numel(cohe) > 1
        bool_failure_phi = 0; % Multiple values of phi
        bool_failure_sweep = 1; % Multiple values of phi and cohesion
    else
        bool_failure_phi = 1;
        bool_failure_sweep = 0;
    end
else
    bool_failure_phi = 0;
    bool_failure_sweep = 0;
end
% -------------------------------------------------------------------------
data=load(strcat(path,nom_file,'.txt'));
r_data=data(:,1);% limites des couches
rho_data=data(:,2);% rho densite de la couche
lambda_data=data(:,3);% lambda coeff de Lam? de la couche
mu_data=data(:,4);% mu rigidite de la couche
V_data=data(:,5);%Viscosite

[r,rho,lambda,mu,V]=TrancheV0(r_data,rho_data,lambda_data,mu_data,V_data,nbtranche,bool_interp_tranche);%d?coupe chaque couche en nbtranche sscouches

%% Model refining
nbcouches=size(mu,2);% nbcouches le nombre couches du mod?le
I=locinterface(mu,nbcouches);%recherche les interfaces solide/liquide,liquide/solide,solide/solide,liquide/liquide
y0=init(n_harmo);%initialisation des trois solutions independante
g=gravity(rho,r,nbcouches,G);% calcul de g pour chaque couche
g_data=gravity(rho_data,r_data,numel(r_data),G);
m2=masse(G,r,g);
press = Pressure(r,rho,g);
press_data = Pressure(r_data,rho_data,g_data);

% if bool_seismic_model
%     Plot_Seismic_Model(r_data,rho_data,lambda_data,mu_data,V_data,per,1);
% end
%
% if bool_geodesic_model
%     Plot_Geodesic_Model(r_data,rho_data,g_data,press_data,1)
% end
%
% if bool_mechanical_model
%     Plot_Mechanical_Model(r_data,rho_data,lambda_data,mu_data,V_data,1);
% end


if per==0
    w=0;
else
    w=1i*2*pi/per; %fr?quence associ? ? per
end


if transfo_laplace==1,w=-w*1i;end%'fr?quence' Laplace

lambda_bck=lambda;
mu_bck=mu;
% --------------------------------------------------------------------------
%solution visco?lastique par principe de correspondance
% -------------------------------------------------------------------------

[mu,lambda]=rheologie_cpx(mu,lambda,V,w,r,type);

% --------------------------------------------------------------------------
%% 2b-Adimensionnement
% --------------------------------------------------------------------------
%ref=[r(nbcouches),rho(1),1e11,g(nbcouches)];%
ref=[r(nbcouches),m2/(4*pi*r(nbcouches)^3),G*m2^2/(4*pi*r(nbcouches)^4),g(nbcouches)];%
r=r/ref(1);
rho=rho/ref(2);
lambda=lambda/ref(3);
mu=mu/ref(3);
g=g/ref(4);
w=w/(sqrt(ref(4)/ref(1)));

%figure(1),subplot(2,2,1),plot(rho,r),title('rho'),subplot(2,2,2),plot(lambda,r),title('lambda'),subplot(2,2,3),plot(mu,r),title('mu'),subplot(2,2,4),plot(V,r),title('V');

%constantes d'adimensionnement
cst1=(ref(2)*ref(4)*ref(1))/ref(3);   %A
cst2=4*pi*ref(2)*G*ref(1)/ref(4)  ;  %B
cst3=ref(4)^2/(4*pi*ref(3)*G)     ;  %C

D=thickness(r,nbcouches);% D ?paisseur de la couche (km)

% --------------------------------------------------------------------------
%% 3-propagation jusqu'a la surface
% --------------------------------------------------------------------------
sol=[];
for sc=1:I(2,1) % super couche (=couches jusqu'a la prochaine interface)
    clear Y;
    M=propag(G,lambda,mu,rho,g,r,D,n_harmo,I,sc,w,cst1,cst2,ref);%inversion de la matrice de la super couche
    for k=1:3 % k solutions ind?pendantes
        C=nmoins1(y0(:,k),sc,I);% conditionement des yi,0 des d?buts de super couche
        condi=condest(M);
        if condi > 1e20
            flag=1;
        end
        Y(:,k)=M\C; %solution pour la super couche sc et le degr? n pour le couple y0k % Where conditioning matters
    end
    sol=[sol;Y];
    clear M, clear C, clear y0;
    if sc~=I(2,1)
        y0=interface(Y,I,mu,sc); %recalcul de yO pour chaque interface.rmq:la derniere interface (sc=I(2,1))correspond a la surface, y0 est calcul? mais inutile et inutilis?
    else

        % --------------------------------------------------------------------------
        %% 4-surface
        % --------------------------------------------------------------------------
        [coeff,determinant]=condsurf(Y,n_harmo,cst3,appli_pression,appli_potentiel); % Where conditioning matters
        out=yiout(sol,coeff,I,mu);%r?arangement des yi le long de r
    end
end
% --------------------------------------------------------------------------
%% 5-sorties
% --------------------------------------------------------------------------
Vpot=m2*G/ref(1);
kpot=1;

if ~bool_greff_lefftz
    % Spherical code from Alterman and Karatekin

    [yi,rp]=decond(out,ref,r,kpot,Vpot); %output for an external potentail=1
    LN=Love(out(:,size(out,2)));%calcul les nombres de love en surface a la fin de l'?tude
    % out are yi adimensionned but yi have the good dimension with a potential
    % correction
else
    % Spheroidal correction from Greff Lefftz
    disp('Warning : yi2, 4 and 6 not usable');
    disp('Warning : yi7 is yi3 wrt phi');
    [dh,dl_theta,dk,dl_phi] = Greff_Lefftz_mod(choice_planet,n_harmo,N_lat,bool_disp_LN_GL);

    % First approx : homothetic transformation
    out(1,:) = out(1,:)*((out(1,end)+dh(N_lat/2+1))/out(1,end));
    out(3,:) = out(3,:)*((out(1,end)+dl_theta(N_lat/2+1))/out(1,end));
    out(5,:) = (out(5,:)-1)*(((out(1,end)-1)+dk(N_lat/2+1))/(out(1,end)-1))+1;

    out = cat(1,out,out(3,:)*((out(1,end)+dl_phi(N_lat/2+1))/out(1,end)));

    disp(strcat('dh=',num2str(dh(N_lat/2+1))))
    disp(strcat('dl_theta=',num2str(dl_theta(N_lat/2+1))))
    disp(strcat('dl_phi=',num2str(dl_phi(N_lat/2+1))))
    disp(strcat('dk=',num2str(dh(N_lat/2+1))))

    [yi,rp]=decond(out,ref,r,kpot,Vpot); %output for an external potentail=1

    yi = cat(1,yi,yi(3,:)*((yi(1,end)+dl_phi(N_lat/2+1))/yi(1,end)));

    LN=Love(out(:,size(out,2)));%calcul les nombres de love en surface a la fin de l'?tude
    % out are yi adimensionned but yi have the good dimension with a potential
    % correction

end

disp(strcat('h =',num2str(LN(1))))
disp(strcat('l =',num2str(LN(2))))
if bool_greff_lefftz
    disp(strcat('l_phi =',num2str(out(7,end))));
end
disp(strcat('k =',num2str(LN(3))))


disp(strcat('y1 =',num2str(yi(1,end))))

% error('Stooop');

if bool_plot_yi
    %figure,for i=1:6,subplot(3,2,i),plot(out(i,:),rp*ref(1)*depth_factor),xlabel(strcat('y',int2str(i))),end
    figure,for i=1:6,subplot(1,6,i),plot(out(i,:),rp*ref(1)*depth_factor),xlabel(strcat('y',int2str(i))),end
end

nn=sqrt((G*(M1+m2)/(a^3)));

% [Ed,Q]=dissip(LN,nn,ref(1),e,G);
% disp(strcat('Q =',num2str(Q)))
% disp(strcat('dE (GW)=',num2str(Ed/1E9)))


% EXTERNAL POTENTIAL and MESH
epot=(nn*ref(1))^2*e; %potentiel subit par le coprs ?tudi? (eccentricity)
opot=(nn*ref(1))^2*obl; % (obliquity)

HarmonicsMesh;
radius_plot = A_ArrayValueToIndex(radius_plot,rp);


% Failure radius simplified in terms of supercouches
couches_idx = 1:nbtranche*2:2*nbcouches+1;

% Omitting transition regions
bot_layers_idx = cat(2,1,couches_idx(2:2:end)) - 1;
top_layers_idx = cat(2,couches_idx(1:2:end-1),couches_idx(end)) + 1;

bot_layers_idx(1) = bot_layers_idx(1) + 1;
top_layers_idx(end) = top_layers_idx(end) - 1;


% surface displacement
Tf=per;
if c_displacement
    Displacement_Locked;
end

%%
% error('Stooop');

if bool_seismic_model
    %Plot_Seismic_Model(ref(1)*r,ref(2)*real(rho),ref(3)*real(lambda),ref(3)*real(mu),real(V),real(per),1);
    Plot_Seismic_Model(ref(1)*r,ref(2)*real(rho),ref(3)*real(lambda),ref(3)*real(mu),real(V),real(per),bool_fus,depth_factor,depth_unit);
end

if bool_geodesic_model
    %Plot_Geodesic_Model(ref(1)*r,ref(2)*real(rho),ref(4)*real(g),press,1);
    Plot_Geodesic_Model(ref(1)*r,ref(2)*real(rho),ref(4)*real(g),press,bool_fus,depth_factor,depth_unit);
end

if bool_mechanical_model
    %Plot_Mechanical_Model(ref(1)*r,ref(2)*real(rho),ref(3)*real(lambda),ref(3)*real(mu),real(V),1);
    Plot_Mechanical_Model(ref(1)*r,ref(2)*real(rho),ref(3)*real(lambda),ref(3)*real(mu),real(V),bool_fus,depth_factor,depth_unit);
end

%% Stress and failure calculations

if c_stress
    Stress;
    if c_failure
        if bool_failure_phi
            Failure_phi;
        elseif bool_failure_sweep
            Failure_Sweep;
        else
            Failure;
            %Failure_Sweep;
        end
    end
end

%% Save files

if bool_data_save && c_stress
    %A_WriteIntoFile;
    tau_m = abs(s1-s3)/2; % Should be positive
    sigma_m = (s1+s3)/2; % Can be either sign
    save(sprintf([nom_file '_MainStress.mat']),'Ntimeloop','Nlon','Nlat','Nradius',...
        'Time','lon','lat','r_s','s1','s3','tau_m','sigma_m');
    save(sprintf([nom_file '_SphericalStress.mat']),'Ntimeloop','Nlon','Nlat','Nradius',...
        'Time','lon','lat','r_s','sigma_rr','sigma_tt','sigma_pp','sigma_rt','sigma_tp','sigma_pr');
    save(sprintf([nom_file '_Others.mat']),'rho_s','r_s','nbcouches','G','cohe','friction',...
        'nom_file','stress_factor','depth_factor','stress_unit','depth_unit',...
        'bool_pressure_failure','legend_pressure','bool_data_save');
end

%g_s =gravity(rho_s,r_s,2*nbcouches+1,G);% calcul de g pour chaque couche
%press_s = Pressure(r_s,rho_s,g_s);
