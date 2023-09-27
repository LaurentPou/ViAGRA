function m=masse(G,r,g)
s=size(g,2);
m=g(s)*r(s)^2/G; %calcul de la masse a partir de g=GM/a^2 (g précalculé par gravity)

