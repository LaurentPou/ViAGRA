function P = Pressure(r,rho,g)

P = zeros(1,numel(rho));

for i=numel(rho)-1:-1:1
   P(i) = P(i+1) + (r(i+1)-r(i))*rho(i)*g(i);
end

end

