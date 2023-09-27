%routine de conversion des valeur Vs, Vp en paramètre de Lamé lambda et mu

path='C:\Documents and Settings\visiteur\Desktop\Nico\ext\';
nom_file='PREM40';
data=load(strcat(path,nom_file,'.txt'));

r=data(:,1);
rho=data(:,2);
vp=data(:,3);
vs=data(:,4);

for i=1:size(vp,1)
mu(i)=rho(i)*vs(i)^2;
lambda(i)=rho(i)*vp(i)^2-2*mu(i);
end

fid=fopen('Prem40.txt','w');
for i=1:size(vp,1)
fprintf(fid,'%f \t %f \t %f \t %f \r',r(i),rho(i),lambda(i),mu(i));
end
fclose(fid);
'conversion ok'
fichier=strcat(path,nom_file,'.txt')