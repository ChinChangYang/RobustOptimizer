function th = angle3d(x,j,i,k)
%
%
a =x(i,:)-x(j,:);
b =x(k,:)-x(j,:);
A =sqrt(a(1)^2+a(2)^2+a(3)^2);
B =sqrt(b(1)^2+b(2)^2+b(3)^2);
c =dot(a,b);
if A*B == 0
    th=pi;
else
    th =acos(c/(A*B));
end