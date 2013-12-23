clear all
global initial_flag
initial_flag=0;
xmax=10;
xmin=-xmax;
vmax=0.4*xmax;
w=0.9
n=20;
c1=2;
c2=2;
m=6;
n_p=n;
I_fno=1;
Fm=0.9;Fn=0.5;
cr=0.9;
itr=30000/n;
for i=1:n
    for j=1:m
      x(i,j)=xmin +rand*(xmax-xmin);
      xbest(i,j)=xmin +rand*(xmax-xmin);
      xgbest(i,j)=xmin +rand*(xmax-xmin);
      v(i,j)=0;
    end
    fitbx(i)=11111111;% local best positions
%     fitb(i)=inf;
end% initialization
gbest=inf;
fitx=[];
for l=1:itr
for i=1:n    
fitx(i)=benchmark_func(x(i,:),I_fno);
end% fitx calculation

kk=fitx<fitbx;
kkc=fitx>fitbx;
xbest=xbest.*repmat(kkc',1,m)+ x.*repmat(kk',1,m);
fitbx=fitbx.*kkc+fitx.*kk;
[p q]=min(fitbx)
if gbest>p
    fitbx=p;
xgbest=repmat(xbest(q,:),n,1);
end

for i=1:n
    for j=1:m
      v(i,j)=w*v(i,j) + c1* rand(1)*(xbest(i,j)-x(i,j)) + c2*rand(1)*(xgbest(1,j)- x(i,j));
      x(i,j)=x(i,j) + v(i,j);
      v(i,j)=max(-vmax,min(v(i,j),vmax));
      v(i,j)=max(-xmax,min(x(i,j),xmax));
    end 
end
p
l
end
theta=2*pi/100;
t=1:0.1:25;
y_t=1+x(1)*sin(x(2)*t*theta+x(3)*sin(x(4)*t*theta+x(5)*sin(x(6)*t*theta)));
    y_0_t=1*sin(5*t*theta-1.5*sin(4.8*t*theta+2*sin(4.9*t*theta)));
  plot(t,y_0_t)
  hold on
    plot(t,y_t,'k')