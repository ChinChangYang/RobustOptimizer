function [fitness PENALTY rate_d]=cost_fn(x)

% global t iter_max;
% Kp = 100*(1-(1-.001)*t/iter_max);
Kp=100;

EBEinputfile;
 n=length(bus_spec(:,1));
Pg=(bus_spec(:,7))/100;
Pd=(bus_spec(:,5))/100;

na=linedata(:,1); nb=linedata(:,2);

g=find(Pg>0);
d=find(Pd>0); 
Pg = Pg(g); 
Pd = Pd(d);

BT=zeros(length(g),length(d));
%%%%%%%%%   100 MW BILATERAL TRANSACTION   %%%%%%%%%%%%%

% BT(1,4)=40;BT(1,5)=10;BT(1,6)=15;
% BT(2,3)=5;BT(2,14)=5;
% BT(3,21)=2.5;
% BT(4,21)=2.5;BT(4,16)=15;
% BT(5,12)=2.5;BT(6,8)=2.5;
% 
% %%%%%%%%%   50 MW BILATERAL TRANSACTION   %%%%%%%%%%%%%
BT(1,4)=5;BT(1,5)=10;BT(1,6)=5;
BT(2,3)=5;
BT(3,21)=2.5;
BT(4,21)=2.5;BT(4,16)=15;
BT(5,12)=2.5;BT(6,8)=2.5;

BT=BT/100;


%% 

% eps = 1;
%%%%%%%%%%%%%%%%%%%%%%%
% BT=zeros(6,21);
Pg2 = sum(BT,2);  %   generations involved in BT 
Pd2 = sum(BT,1);  %    loads involved in BT
%%%%%%%%%%%%%%%


% %% Calculation of GD
GD = zeros(length(g),length(d));
for i=1:length(g);
    GD(i,:)=x( ((i-1)*length(d)+1) : (i*length(d)) );
end

%% calculation of PTDF

% xi=imag(linedata(:,3));
% flows

line_data=linedata;
%line_data(:,3)=line_data(:,3).*ranline(:,ns);
[YIbus] = EBEformybus(line_data,n);
YIbus(1,:)=[];
YIbus(:,1)=[];
Xt=inv(YIbus);
X=zeros(n,n);
for i=0:(n-1)
    for j=0:(n-1)
    if (i~=0) && (j~=0)
        X(i+1,j+1)=Xt(i,j);
    else
        X(i+1,j+1)=0;
    end
    end
end

% PTDF 
xij=imag(line_data(:,3));


for i=1:length(g)
    for j=1:length(d)
        for k=1:length(na)
        PTDF(i,j,k)= [X(na(k),g(i))-X(nb(k),g(i))-X(na(k),d(j))+X(nb(k),d(j))]/xij(k);
        end
    end
end

% GD
% PTDF


%% calculation of objective fn.

% Rg=[4.1925    4.1135    3.8882    4.3101    8.2130    8.5662];

Rg = [32.7290   32.1122   30.3532   33.6474   64.1156   66.8729];

Rd = [7.5449 10.7964 10.9944 11.0402 11.7990 15.3803 42.6800 41.4551 ...
      73.1939 57.0430 45.5920 43.6553 61.8002 59.6409 57.0279 ...
      51.0749 67.1070 60.6623 198.6744  178.9956  199.9483];



% FC=100*xij;
FC=100*xij/sum(xij);

flows=zeros(length(na));
cost_line=zeros(length(na));

for k=1:length(na)
    for i=1:length(g)
        for j=1:length(d)
        
            flows(k)=flows(k) + abs(PTDF(i,j,k)*GD(i,j)) + abs(PTDF(i,j,k)*BT(i,j));
           
        end
    end
    cost_line(k)=FC(k)/flows(k);
end


pg = Pg - Pg2;
% Pd'
% Pd2'
pd = Pd - Pd2';

rate_ebe1=0;
cost_l=zeros(length(g),length(d));
cost_gen=zeros(length(g));

for i=1:length(g)
    for j=1:length(d)
        for k=1:length(na)
            
            cost_l(i,j) = cost_l(i,j) + abs(cost_line(k)*PTDF(i,j,k));
        end
        
        cost_gen(i) = cost_gen(i) + GD(i,j)*cost_l(i,j);
    end
    rate_ebe1=rate_ebe1 + (cost_gen(i)/pg(i) - Rg(i))^2;
end
  
rate_ebe2=0;
% cost_l=zeros(length(g),length(d));
cost_load=zeros(length(d));

for j=1:length(d)
    for i=1:length(g)
%         for k=1:length(na)
%             
%             cost_l(i,j) = cost_l(i,j) + cost_line(k)*PTDF(i,j,k);
%         end
%         
        cost_load(j) = cost_load(j) + GD(i,j)*cost_l(i,j);
    end
    rate_ebe2=rate_ebe2 + (cost_load(j)/pd(j) - Rd(j))^2;
end
        
rate_d = rate_ebe1 + rate_ebe2;

%%    CONSTRAINT  VIOLATIONS  

Pg_x = sum(GD,2)+ sum(BT,2);
Pd_x = sum(GD,1)+ sum(BT,1);

% 100*sum(sum(BT,2))

Gpen = 0;
for i=1:length(g)
%     if abs(Pg_x(i)-Pg(i)) >=eps
%         Gpen = Gpen + MM;
%     else
        Gpen = Gpen + 100*abs(Pg_x(i)-Pg(i));
%     end
end
% [Pg_x Pg]


LDpen = 0;
for i=1:length(d)
%     if abs(Pd_x(i)-Pd(i)) >=eps
%         LDpen = LDpen + MM;
%     else
        LDpen = LDpen + 100*abs(Pd_x(i)-Pd(i));
%     end
end

PENALTY  = Gpen + LDpen;

% Gpen
% LDpen
fitness = rate_d + 50*Kp*PENALTY;
% save('test_c.mat');

