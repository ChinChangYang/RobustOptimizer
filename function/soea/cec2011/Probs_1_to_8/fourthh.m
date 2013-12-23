tol=1.0e-08;
tspan=[0 0.78];
yo =[ 0.09 0.09 ]';
u=2;
options = odeset('RelTol',tol);
[T,Y] = ode45(@(t,y) rigid(t,y,u),tspan,yo,options);
j=sum(sum(Y.^2,2)+(0.1)*u*u)