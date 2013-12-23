function[y sllreturn bwfn]=antennafunccircular(x1,null,phi_desired,distance)
%This function calculates the fitness value of array 'x1' which is returned in 'y'
global directivity

pi=3.141592654;
dim=length(x1);
y=0;
[temp num_null]=size(null);
num1=300;
phi=linspace(0,360,num1);
phizero=0;
yax(1)=array_factorcir(x1,(pi/180)*phi(1),phi_desired,distance,dim);
maxi=yax;
phi_ref=1;
for i=2:num1%This loop finds out the maximum gain 
    yax(i)=array_factorcir(x1,(pi/180)*phi(i),phi_desired,distance,dim);
    if maxi<yax(i)
        maxi=yax(i);
        phizero=phi(i);
        phi_ref=i;
    end;
end;
maxtem=0;
count=0;
if yax(1)>yax(num1) && yax(1)>yax(2)
    count=count+1;
    sidelobes(count)=yax(1);
    sllphi(count)=phi(1);
end
if yax(num1)>yax(1) && yax(num1)>yax(num1-1)
    count=count+1;
    sidelobes(count)=yax(num1);
    sllphi(count)=phi(num1);
end
for i=2:num1-1
    if yax(i)>yax(i+1) && yax(i)>yax(i-1)
        count=count+1;
        sidelobes(count)=yax(i);
        sllphi(count)=phi(i);
    end
end
sidelobes=sort(sidelobes,'descend');
upper_bound=180;
lower_bound=180;
y=sidelobes(2)/maxi;
sllreturn=20*log10(y);
for i=1:num1/2
    if (phi_ref+i)>num1-1
        upper_bound=180;
        break;
    end
    tem=yax(phi_ref+i);
    if yax(phi_ref+i)<yax(phi_ref+i-1) && yax(phi_ref+i)<yax(phi_ref+i+1)
        upper_bound=phi(phi_ref+i)-phi(phi_ref);
        break;
    end;
end

for i=1:num1/2
    
    if (phi_ref-i<2)
        lower_bound=180;
        break;
    end
    tem=yax(phi_ref-i);
    if yax(phi_ref-i)<yax(phi_ref-i-1) && yax(phi_ref-i)<yax(phi_ref-i+1)
        lower_bound=phi(phi_ref)-phi(phi_ref-i);
    break;
    end;
end
bwfn=upper_bound+lower_bound;
% bwfn
%y=maxtem;
y1=0;
for i=1:num_null%The objective function for null control
    %is calculated here
    y1=y1+(array_factorcir(x1,null(i),phi_desired,distance,dim)/maxi);
end;
y2=0;
uavg=trapezoidalcir(x1,0,2*pi,50,phi_desired,distance,dim);
y2=abs(2*pi*maxi*maxi/uavg);
directivity=10*log10(y2);
y3=abs(phizero-phi_desired);
if y3<5
    y3=0;
end
y=0;
if bwfn>80
    y=y+abs(bwfn-80);
end;
% directivity
%y3=abs((phizero-phi_desired)*(pi/180));
%y=y+y1+1/y2+y3;
y=sllreturn+y+y1+y3;