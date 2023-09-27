function [Ed,Q]=dissip(LN,n,R,e,G)

imk2=imag(LN(3));
Ed=-21/2*n^5*R^5*e^2/G*imk2;
Q=-(abs(LN(3))/imk2);
