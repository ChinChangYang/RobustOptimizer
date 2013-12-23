function dy = intgrl(t,x,u)

exp(25*pinv(x+2)*x)
exp(25*x*pinv(x+2))
pause
dy = zeros(2,1);    % a column vector
dy(1) = -(2.+u)*(x+0.25)+(x+0.5)*exp(25*pinv(x+2)*x);
dy(2) = 0.5-x-(x+0.5)*exp(25*x/(x+2));
