function dy = intgrl(t,x,u)
dy = zeros(2,1);    % a column vector
dy(1) = -(2+u)*(x(1)+0.25)+(x(2)+0.5)*exp(25*x(1)/(x(1)+2));
dy(2) = 0.5-x(2)-(x(2)+0.5)*exp(25*x(1)/(x(1)+2));
