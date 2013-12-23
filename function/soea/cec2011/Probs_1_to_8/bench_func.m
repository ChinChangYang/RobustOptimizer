function f=bench_func(x,fun_num)
[ps,d]=size(x);

%fun_num=1   Parameter Estimation for Frequency-Modulated (FM) Sound Waves,initialization range=[0,6.35], bound=[-6.4,6.35] , length of x=6.                                                                                   
%fun_num=2   Lennard Jones potential value. x passed to this function must be n dimentional array where, n is perfectly divisible by 3.
%fun_num=3   The bifunctional catalyst blend optimal control problem.this is a one dimensional problem with bound=[0.6,0.9]. relative tolerance is set to 1e-001.
%fun_num=4   The optimal control of a non-linear stirred tank reactor. this is a one dimensional problem with no bound and initialization range=[0,5.0].relative tolerance is set to 1e-001.
%fun_num=5   Tersoff potential for model Si(B).
%fun_num=6   Tersoff potential for model Si(C).
%fun_num=7   spread spectrum radar polyphase problem. bound=[0,2*pi].MAX NO: OF DIMENSIONS=20.

%fun_num=8   Transmission Network Expansion Planning (TNEP) problem. The
%            objective function is to be minimized subject to the constaraints given in
%            the problem defination file. Maximum and minimum bound of the population
%            that to be generated is 0 and 15. The minimum dimension should
%            be 7.

switch fun_num
    case 1

           if d<6
              disp('dimension-size should be six.')
           else
           if d>6
              disp('dimension-size is more than 6.')
              disp('function has been evaluated on first six dimensions.')
           end
           end
         theta=2*pi/100;
         f=0;
           for t=0:100
               y_t=x(1)*sin(x(2)*t*theta+x(3)*sin(x(4)*t*theta+x(5)*sin(x(6)*t*theta)));
               y_0_t=1*sin(5*t*theta-1.5*sin(4.8*t*theta+2*sin(4.9*t*theta)));
               f=f+(y_t-y_0_t)^2;
           end


case 2
    % lennard jones potential problem.
% x passed to this function must be n dimentional array where, n is...
% perfectly divisible by 3.
          r=[];
          p=size(x);
          if rem(p(2),3)~=0
             disp('x passed to this function must be n dimentional array where, n is perfectly divisible by 3.')
          end
           n=p(2)/3;
           x=reshape(x,3,n)';
           v=0;
           a=ones(n,n);
           b=repmat(2,n,n);
         for i=1:(n-1)
           for j=(i+1):n
            r(i,j)=sqrt(sum((x(i,:)-x(j,:)).^2));
            v=v+(a(i,j)/r(i,j)^12-b(i,j)/r(i,j)^6);
           end
         end
         f=v;

case 3
    %The bifunctional catalyst blend optimal control problem
          tol=1.0e-01;% tol is a matter of concern. decreasing it make algo fast. 
          tspan=[0 0.78];% check for the range
          yo =[1 0 0 0 0 0 0];
          u=x;%u should be passed here.
          options = odeset('RelTol',tol);
       [T,Y] = ode45(@(t,y) diffsolv(t,y,u),tspan,yo,options);
       w=size(Y);
       f=Y(w(1),w(2))*1e+003;
 
case 4
% The optimal control of a non-linear stirred tank reactor.       
          tol=1.0e-01;
          tspan=[0 0.78];
          yo =[ 0.09 0.09]';
          u=x;%u should be passed here.
          options = odeset('RelTol',tol);
       [T,Y] = ode45(@(t,y) intgrl(t,y,u),tspan,yo,options);
          f=sum(sum(Y.^2,2)+(0.1)*(u).*(u));



case 5
    % tersoff potential problem, model Si(B)
    p=size(x);
if rem(p(2),3)~=0
    disp('x passed to this function must be n dimentional array where, n is perfectly divisible by 3.')
end
NP=p(2)/3;
x=reshape(x,3,NP)';
R1=3.0;
R2=0.2;
A=3.2647e+3;
B=9.5373e+1;
lemda1=3.2394;
lemda2=1.3258;
lemda3=1.3258;
c=4.8381;
d=2.0417;
n1=22.956;
gama=0.33675;
h=0;
E=zeros(1,NP);
r=zeros(NP);
for i=1:NP
    for j=1:NP
% %         if i==j
% %             continue
% %         end
         r(i,j)=sqrt(sum((x(i,:)-x(j,:)).^2));
if r(i,j)<(R1-R2)
    fcr(i,j)=1;
elseif  r(i,j)>(R1+R2)
    fcr(i,j)=0;
else
    fcr(i,j)=0.5-0.5*sin(pi/2*(r(i,j)-R1)/R2);
end

VRr(i,j)=A*exp(-lemda1*r(i,j));
VAr(i,j)=B*exp(-lemda2*r(i,j));
    end
end
for i=1:NP
    for j=1:NP
        if i==j
          continue
        end
jeta=zeros(NP,NP);
for k=1:NP
    if i==k || j==k        continue  
    end
     rd1=sqrt(sum((x(i,:)-x(k,:)).^2));
      rd3=sqrt(sum((x(k,:)-x(j,:)).^2));
       rd2=sqrt(sum((x(i,:)-x(j,:)).^2));
       ctheta_ijk=(rd1^2+rd2^2-rd3^3)/(2*rd1*rd2);
       G_th_ijk =1+(c^2)/(d^2)-(c^2)/(d^2+(h-ctheta_ijk)^2);
       jeta(i,j)=jeta(i,j)+fcr(i,k)*G_th_ijk*exp(lemda3^3*(r(i,j)-r(i,k))^3);
end

Bij=(1+(gama*jeta(i,j))^n1)^(-0.5/n1);
E(i)=E(i)+fcr(i,j)*(VRr(i,j)-Bij*VAr(i,j))/2;
    end    
end
f=sum(E);

case 6
    % tersoff potential for model Si(C)
         p=size(x);
if rem(p(2),3)~=0
    disp('x passed to this function must be n dimentional array where, n is perfectly divisible by 3.')
end
NP=p(2)/3;
x=reshape(x,3,NP)';
R1=2.85;
R2=0.15;
A=1.8308e+3;
B=4.7118e+2;
lemda1=2.4799;
lemda2=1.7322;
lemda3=1.7322;
c=1.0039e+05;
d=1.6218e+01;
n1=7.8734e-01;
gama=1.0999e-06;
h=-5.9826e-01;
E=zeros(1,NP);
r=zeros(NP);
for i=1:NP
    for j=1:NP
         r(i,j)=sqrt(sum((x(i,:)-x(j,:)).^2));
if r(i,j)<(R1-R2)
    fcr(i,j)=1;
elseif  r(i,j)>(R1+R2)
    fcr(i,j)=0;
else
    fcr(i,j)=0.5-0.5*sin(pi/2*(r(i,j)-R1)/R2);
end

VRr(i,j)=A*exp(-lemda1*r(i,j));
VAr(i,j)=B*exp(-lemda2*r(i,j));
    end
end
for i=1:NP
    for j=1:NP
        if i==j
          continue
        end
jeta=zeros(NP,NP);
for k=1:NP
    if i==k || j==k        continue  
    end
     rd1=sqrt(sum((x(i,:)-x(k,:)).^2));
      rd3=sqrt(sum((x(k,:)-x(j,:)).^2));
       rd2=sqrt(sum((x(i,:)-x(j,:)).^2));
       ctheta_ijk=(rd1^2+rd2^2-rd3^3)/(2*rd1*rd2);
       G_th_ijk =1+(c^2)/(d^2)-(c^2)/(d^2+(h-ctheta_ijk)^2);
       jeta(i,j)=jeta(i,j)+fcr(i,k)*G_th_ijk*exp(lemda3^3*(r(i,j)-r(i,k))^3);
end

Bij=(1+(gama*jeta(i,j))^n1)^(-0.5/n1);
E(i)=E(i)+fcr(i,j)*(VRr(i,j)-Bij*VAr(i,j))/2;
    end    
end
f=sum(E);

case 7
% sprd_spectrum_rad_pphase
          hsum=[];
          var=2*d-1;
            for kk=1:2*var
              if rem(kk,2)~=0
                 i=(kk+1)/2;
                 hsum(kk)=0;
              for j=i:d    %fi(2i-1)X
                 summ=0;
                  for i1=(abs(2*i-j-1)+1):j
                     summ=x(i1)+summ;
                  end
                 hsum(kk)= cos(summ)+ hsum(kk);
              end
              else 
                i=kk/2;
                hsum(kk)=0;
              for j=(i+1):d    %fi(2i)X
                summ=0;
                 for i1=(abs(2*i-j)+1):j
                   summ=x(i1)+summ;
                 end
                hsum(kk)= cos(summ)+ hsum(kk);
              end
               hsum(kk)=hsum(kk)+0.5;
              end
            end
                f=max(hsum) ;  

   
case 8
    %Transmission Network Expansion Planning 
    % Maximum and minimum bound of the population that to be generated is 1 and 16.  
    sw=ceil(x);
    data6Bus;
    n1=length(Linedata(:,1));
    sw1=sw;
    for k=1:length(sw1)
        Linedata(n1+k,:)=Candidate(sw1(k),:);
    end
    n_orginalLine=n1;
    n=length(Pgen);
    B=zeros(n,n);
    Nline=length(Linedata(:,1));
    %Ncand=max(Candidate(:,1));
    Xline=Linedata(:,4);
    pijmax=Linedata(:,6);
    Tap=ones(n);
    for C=1:Nline
        bline(C)=1/Xline(C);
        k=Linedata(C,2);
        m=Linedata(C,3);
        B(k,m)=B(k,m)-(bline(C));
        B(m,k)=B(k,m);
        B(k,k)=B(k,k)+(bline(C));
        B(m,m)=B(m,m)+(bline(C));
    end    
B(1,1)=10000000;
X=inv(B);
delP= Pgen-Pload;
delP=(delP');

delta=X*(delP);
pij=zeros(Nline,1);
for k=1:Nline
i=Linedata(k,2);
j=Linedata(k,3);
pij(k)=(delta(i)-delta(j))/Xline(k);
end
PIPbase=0.0;
f=sum(Linedata(n_orginalLine+1:end,7))+30;
pen=0;

for i=1:length(Linedata(:,1))
    pen=pen+5000*max((abs(pij(i))-Linedata(i,6)),0);
end

for i=1:length(Candidate(:,1))
    [a ]=find(sw==i);
    if length(a)>3
        pen=pen+1000;
    end
end
f=f+pen;
    
otherwise % in no valid case
        disp('not valid function no: is specified')
end   