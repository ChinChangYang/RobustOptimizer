function dy = diffsolv(t,x,u)
load c_bifunc_data;% c(i,j) is saved here.
ml=[1 u u.^2 u.^3];
mlt=repmat(ml,10,1);
k=sum(c.*mlt,2);
dy = zeros(7,1);    % a column vector
dy(1) = -k(1)*x(1);
dy(2) = k(1)*x(1)-(k(2)+k(3))*x(2)+k(4)*x(5);
dy(3) = k(2)*x(2);
dy(4) = -k(6)*x(4)+k(5)*x(5);
dy(5) = k(3)*x(2)+k(6)*x(4)-(k(4)+k(5)+k(8)+k(9))*x(5)+k(7)*x(6)+k(10)*x(7);
dy(6) = k(8)*x(5)-k(7)*x(6);
dy(7) = k(9)*x(5)-k(10)*x(7);
