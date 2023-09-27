function [mu,lambda]=rheologie_cpx(mu,lambda,V,w,r,type)

%rheologie  lineaire
%------------------------------------------------------------------
% Maxwell (1)
if type==1
    %disp ('rheologie de Maxwell')
    for i=1:size(r,2)
        if V(i)~=0
            mu(i)=w*mu(i)/(w+mu(i)/V(i));
            %lambda(i)=lambda(i)*((w+mu(i)/V(i)*(2*mu(i)/3*lambda(i)+1))/(w+mu(i)/V(i)));% pour solution compressible visco-ela...
        else
            mu(i)=mu(i);
            %lambda(i)=lambda(i);
        end
    end
end
% Kelvin-Voight (2)
if type==2
    %disp ('rheologie de Kelvin-Voight')
    for i=1:size(r,2)
        if V(i)~=0
            mu(i)=mu(i)*(1+(V(i)*w)/mu(i));
            %lambda(i)=lambda(i)*((w+mu(i)/V(i)*(2*mu(i)/3*lambda(i)+1))/(w+mu(i)/V(i)));% pour solution compressible visco-ela...
        else
            mu(i)=mu(i);
            %lambda(i)=lambda(i);
        end
    end
end

%Standard liner Solid (SLS) (3)
if type==3
    %disp ('rheologie Standard liner Solid (SLS)')
    for i=1:size(r,2)
        if V(i)~=0
            mu1(i)=mu(i);
            mu2(i)=0.15*mu(i);
            mu(i)=(mu1(i)*(w+(mu2(i)/V(i))))/(w+(mu1(i)+mu2(i))/V(i));
            %lambda(i)=lambda(i)*((w+mu(i)/V(i)*(2*mu(i)/3*lambda(i)+1))/(w+mu(i)/V(i)));% pour solution compressible visco-ela...
        else
            mu(i)=mu(i);
            %lambda(i)=lambda(i);
        end
    end
end

%Bourger's Body (4)
if type == 4
    %disp ('rheologie Bourger s Body')
    for i=1:size(r,2)
        if V(i)~=0
            mu1(i)=mu(i);
            mu2(i)=0.15*mu(i);
            V1(i)=V(i);
            V2(i)=0.55*V(i);
            mu(i)=(mu1(i)*w^2+w*(mu1(i)*mu2(i))/V2(i))/(w^2+(mu1(i)/V1(i)+mu2(i)/V2(i)+mu1(i)/V2(i))*w+(mu1(i)*mu2(i))/(V1(i)*V2(i)));
            %lambda(i)=lambda(i)*((w+mu(i)/V(i)*(2*mu(i)/3*lambda(i)+1))/(w+mu(i)/V(i)));% pour solution compressible visco-ela...
        else
            mu(i)=mu(i);
            %lambda(i)=lambda(i);
        end
    end
end
%rheologie non lineaire
%------------------------------------------------------------------
%Caputo (5)
if type ==5
    %disp ('rheologie Caputo')
    for i=1:size(r,2)
        if V(i)~=0
            z=0.7;
            mu(i)=w^z*mu(i)/(w^z+mu(i)/V(i));
            %lambda(i)=lambda(i)*((w+mu(i)/V(i)*(2*mu(i)/3*lambda(i)+1))/(w+mu(i)/V(i)));% pour solution compressible visco-ela...
        else
            mu(i)=mu(i);
            %lambda(i)=lambda(i);
        end
    end
end