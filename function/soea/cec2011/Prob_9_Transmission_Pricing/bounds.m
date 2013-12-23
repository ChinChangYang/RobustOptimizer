

EBEinputfile; 

% n=length(bus_spec(:,1));
Pg=(bus_spec(:,7))/100;
Pd=(bus_spec(:,5))/100;

g=find(Pg>0);
d=find(Pd>0); 

%%%%%%%%% define BT   %%%%%%%%%%%%%%%%%%%%%
BT=zeros(length(g),length(d));

BT(1,4)=5;BT(1,5)=10;BT(1,6)=5;
BT(2,3)=5;
BT(3,21)=2.5;
BT(4,21)=2.5;BT(4,16)=15;
BT(5,12)=2.5;BT(6,8)=2.5;

BT=BT/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GD_max=zeros(length(g),length(d));
for i=1:length(g)
    for j=1:length(d)
        GD_max(i,j)=min(Pg(g(i))-BT(i,j),Pd(d(j))-BT(i,j));
    end
end



