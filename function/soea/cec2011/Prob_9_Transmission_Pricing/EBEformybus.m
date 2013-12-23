% Formation of Ybus.
function [YIbus] = formybus(linedata,n)
busa=linedata(:,1);busb=linedata(:,2);
L=length(busa);
Z=linedata(:,3);ZI=imag(Z);
YI=ones(L,1)./ZI;
YIbus=zeros(n,n);

    for k=1:L
    
    n=busa(k);m=busb(k);
    YIbus(n,n)=YIbus(n,n)+YI(k);
    YIbus(n,m)=-YI(k)+YIbus(n,m);
    YIbus(m,n)=-YI(k)+YIbus(m,n);
    YIbus(m,m)=YIbus(m,m)+YI(k);    
    end

return;