clear all

planete='Prem5c_II';
nbcouches=30;%pas d'intergration
pres=1; %pression
poten=1;%potentiel
Lapla=1;%transfo de laplace
yr=60*60*24*365.25;%année
degre=[2]%,6,15,30,100]%,20,25,30,50,70,90,100];

calc_t=1;
visu_s=1;
visu_det=0;

%***************************************************************************************
%spectre de viscoéla (en seconde)
%***************************************************************************************
s=1;
for i=1:8
   for j=[1,3] 
spectre(s)=j*1000*yr*str2num(strcat('1e',num2str(i-4)));
s=s+1;
   end
end
svalues=fliplr(1./spectre);
%***************************************************************************************


disp(strcat('modèle:',planete))

for ind=1:length(degre)%boucle sur n
deg=degre(ind);    
disp(strcat('degré ',int2str(deg),'...'))

%***************************************************************************************
disp('réponse élastique (s)...')           %réponse dirac du système
%***************************************************************************************
[LNs_e,yis_e,determinant,flag_e]=elastique(planete,deg,nbcouches,pres,poten,Lapla);
LNs_e(:,deg)=LNs_e;
disp('    h         k         l ')
disp(LNs_e(:,deg)')



%***************************************************************************************
disp('réponse viscoélastique (s)...')      %calcul de la ln_e-ln_ve pour i svalues
%***************************************************************************************
for i=1:size(spectre,2) 
    if i/5==floor(i/5)
    disp(strcat(num2str(floor(i/size(spectre,2)*100)),'%'))
    end
    [LNs_ve,yis_ve,determinant,flag_ve]=viscoelastique(planete,spectre(i),deg,nbcouches,pres,poten,Lapla);
 LNs_v_temp(i,1)=LNs_ve(1)-LNs_e(1,deg);
 LNs_v_temp(i,2)=LNs_ve(2)-LNs_e(2,deg);
 LNs_v_temp(i,3)=LNs_ve(3)-LNs_e(3,deg);
 deter(i)=determinant;
    for j=1:6
    yis_v_temp(i,j)=yis_ve(j)-yis_e(j);
    end
end
D(:,deg)=deter';
LNs_v(:,:,deg)=(LNs_v_temp);
yis_v(:,:,deg)=(yis_v_temp);




if calc_t==1;
%***************************************************************************************
disp('résidues (s)...')    %calcul des résidus sur la base quadratique des svalues
%***************************************************************************************


imax=size(spectre,2);
for i=1:imax
for j=1:imax
m(i,j)=((1/svalues(i))*(1/svalues(j)))/((1/svalues(i))+(1/svalues(j)));%generation matrice quadratique
end
end

for i=1:imax %conditionnement
    m(i,:)=m(i,:)/max(m(i,:));
    LNs_v(i,:,deg)=LNs_v(i,:,deg)/max(m(i,:));
    yis_v(i,:,deg)=yis_v(i,:,deg)/max(m(i,:));
end

hresidues(:,deg)=flipud(-m\LNs_v(:,1,deg));
kresidues(:,deg)=flipud(-m\LNs_v(:,2,deg));
lresidues(:,deg)=flipud(-m\LNs_v(:,3,deg));
for i=1:6
yresidues(:,i,deg)=flipud(-m\yis_v(:,i,deg));
end


%***************************************************************************************
disp('réponse viscoélastique (t)...')    %calcul de la réponse viscoélastique f(t)
%***************************************************************************************

t=1;
for i=1:10
   for j=[1,1.2,1.4,1.8,2,2.5,3,4,6,8] 
tps(t)=j*yr*str2num(strcat('1e',num2str(i)));
t=t+1;
   end
end
    
for t=1:size(tps,2) %Laplace inverse
        loveh(t,deg)=hresidues(1,deg)*exp(-tps(t)*svalues(1));
        lovek(t,deg)=kresidues(1,deg)*exp(-tps(t)*svalues(1));
        lovel(t,deg)=lresidues(1,deg)*exp(-tps(t)*svalues(1));
        for k=1:6
        y(t,k,deg)=yresidues(1,k,deg)*exp(-tps(t)*svalues(1));
        end
    for i=2:imax
        loveh(t,deg)=loveh(t,deg)+hresidues(i,deg)*exp(-tps(t)*svalues(i));
        lovek(t,deg)=lovek(t,deg)+kresidues(i,deg)*exp(-tps(t)*svalues(i));
        lovel(t,deg)=lovel(t,deg)+lresidues(i,deg)*exp(-tps(t)*svalues(i));
        for k=1:6
        y(t,k,deg)=y(t,k,deg)+yresidues(i,k,deg)*exp(-tps(t)*svalues(i));
        end
    end
end


end % if calc_t
end % boucle sur deg


%***************************************************************************************
disp('plot')
%***************************************************************************************
%r Red
%g Green
%b Blue
%c Cyan
%m Magenta
%y Yellow
%k Black 
%w White 
couleur=char('b','c','g','y','r','m','k','b','c','g','y','m','k','r','m','k','b','c','g','y','m','k','r');
Tscale=1000*yr;

%---------------------------------------------------------------------------------------
if visu_s==1;

figure %spectre des nombres de love
love=1;
titre=strcat('spectre ln_e_e(',int2str(love),')-ln_e(',int2str(love),')',planete);
subplot(3,1,love),plot(spectre/Tscale,LNs_v(:,love,2)/(2*2+1),'k*-')
title(titre),xlabel('spectre (ka)'),ylabel('h/(2n+1)')
set(gca,'xscale','log')%set(gca,'yscale','log')
hold on
 for deg=2:length(degre)
    plot(spectre/Tscale,LNs_v(:,love,degre(deg))/(2*degre(deg)+1),couleur(deg-1))
 end
love=3;
titre=strcat('spectre ln_e_e(',int2str(love),')-ln_e(',int2str(love),')',planete);
subplot(3,1,2),plot(spectre/Tscale,LNs_v(:,love,2),'k*-')
title(titre),xlabel('spectre (ka)'),ylabel('l')
set(gca,'xscale','log')%set(gca,'yscale','log')
hold on
 for deg=2:length(degre)
    plot(spectre/Tscale,LNs_v(:,love,degre(deg)),couleur(deg-1))
end
love=2;
titre=strcat('spectre ln_e_e(',int2str(love),')-ln_e(',int2str(love),')',planete);
subplot(3,1,3),plot(spectre/Tscale,LNs_v(:,love,2)/(2*2+1),'k*-')
title(titre),xlabel('spectre (ka)'),ylabel('k')
set(gca,'xscale','log')%set(gca,'yscale','log')
hold on
 for deg=2:length(degre)
    plot(spectre/Tscale,LNs_v(:,love,degre(deg))/(2*degre(deg)+1),couleur(deg-1))
 end


if visu_det==1; 
 
figure %spectre du determinant
for deg=1:length(degre)
plot(spectre/(Tscale),D(:,degre(deg)),couleur(deg))
title('spectre du determinant'),xlabel('spectre (ka)'),ylabel('det')
set(gca,'xscale','log')%set(gca,'yscale','log')
hold on
end


figure %plot des svalues (det==0)
flag2=0;
for deg=1:length(degre)
    clear zero;
    ind=1;
    for i=2:size(D,1)
        if D(i-1,degre(deg))*D(i,degre(deg))<0 
        zero(ind)=i;
        ind=ind+1;
        flag2(deg)=1;
        end
    end
    if flag2(deg)==1;
    for i=1:length(zero)
    svalu(i,degre(deg))=spectre(zero(i));
    end
    plot(degre(deg),svalu(:,degre(deg))/Tscale,'*')
    title('spectre des zero du determinant'),    xlabel('degre d harmonique'),    ylabel('modes propres')
    set(gca,'yscale','log')
    hold on
    end
end
hold off

end %if visu_det

end %if visu_s
%---------------------------------------------------------------------------------------


%---------------------------------------------------------------------------------------
if calc_t==1;
    
figure
for i=1:6
subplot(4,3,i),plot(spectre/(1000*yr),yis_v(:,i,2)),set(gca,'xscale','log')
title('spectres yi'),xlabel('spectre (ka)'),ylabel(strcat('y',int2str(i)))
set(gca,'xscale','log')
subplot(4,3,i+6),plot(svalues,yresidues(:,i,2)),set(gca,'xscale','log'),ylabel(strcat('y_residus',int2str(i)))%,set(gca,'yscale','log')
title('spectres residus yi')
end
 
 
%-----------          --------------          ------------            ----------
figure
subplot(3,1,1)
plot(tps/(yr*1000),(loveh(:,2)+LNs_e(1,2))/(2*2+1),'k')
title(strcat('evotution temporelle pour  -',planete)),xlabel('temps (kya)'),ylabel('h/(2n+1)')
set(gca,'xscale','log');
hold on
 for deg=2:length(degre)
    plot(tps/(yr*1000),(loveh(:,degre(deg))+LNs_e(1,degre(deg)))/(2*degre(deg)+1),couleur(deg-1))
end
subplot(3,1,2)
plot(tps/(yr*1000),lovel(:,2)+LNs_e(3,2),'k')
title(strcat('evotution temporelle pour  -',planete)),xlabel('temps (kya)'),ylabel('l')
set(gca,'xscale','log');
hold on
 for deg=2:length(degre)
    plot(tps/(yr*1000),lovel(:,degre(deg))+LNs_e(3,degre(deg)),couleur(deg-1))
 end

subplot(3,1,3)
plot(tps/(yr*1000),lovek(:,2)+LNs_e(2,2),'k')
title(strcat('evotution temporelle pour  -',planete)),xlabel('temps (kya)'),ylabel('k')
set(gca,'xscale','log');
hold on
 for deg=2:length(degre)
    plot(tps/(yr*1000),lovek(:,degre(deg))+LNs_e(2,degre(deg)),couleur(deg-1))
 end

 
figure
for k=1:6
subplot(3,2,k),plot(tps/(yr*1000),y(:,k,2))
title(strcat('evotution temporelle pour  -',planete)),xlabel('temps (kya)'),ylabel(strcat('y',int2str(k)))
set(gca,'xscale','log')%,set(gca,'yscale','log');
hold on
 for deg=2:length(degre)
    plot(tps/(yr*1000),y(:,k,degre(deg)),couleur(deg-1))
 end
end

end %calc_t
%---------------------------------------------------------------------------------------

%fin
%***************************************************************************************
