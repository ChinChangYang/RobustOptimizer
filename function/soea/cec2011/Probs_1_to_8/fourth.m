function f=fourth(x,u) 
[ps d]=size(x)
f=sum(Y.^2,2)+0.1*repmat(u*u,ps,1);
% f(i)=x(1)^2+x(2)^2+0.1*u^2;
