clear all
global initial_flag
initial_flag=0;
xmax=10;
xmin=-xmax;
w=0.9
n=20;
m=6;
n_p=n;
I_fno=1;
Fm=0.9;Fn=0.5;
cr=0.9;
rep=30000/n;
for i=1:n
    for j=1:m
      x(i,j)=xmin +rand*(xmax-xmin)
    end
%     fitb(i)=inf;
end% initialization
fitx=[];
for i=1:n    
fitx(i)=bench_func(x(i,:),I_fno);
end% fitx calculation
pop_in=randint(n_p,m,[1,n]);
for i=1:n_p
pvar=x([pop_in(i,:)],1:m);
p(i,:)=diag(pvar)';
end % [population generation based on permu&combn

for itr=1:rep
    F=0.8;
% for i=1:n_p    
% fitb(i)=benchmark_func(p(i,:),I_fno);
% end% fitb calculation
% [pp qq]=min(fitb);
% mut(1,:)=p(qq,:);

%     [pp qq]=min(fitx);
%     mutt(1,:)=p(qq,:);
% mut(1,:)=p(qq,:);
for i=1:n
    k1=ceil(rand*n_p);
    while k1==i 
        k1=ceil(rand*n_p);
    end 
    k2=ceil(rand*n_p);
        while k2==k1 | k2==i 
            k2=ceil(rand*n_p);
        end  
        k3=ceil(rand*n_p);
        while k3==k1 | k3==k2 | k3==i k3=ceil(rand*n_p); end
    for j=1:m
% rn=randint(1,3,[1,n_p]);
% x_new(i,j)=mut(j)+ F*(p(rn(2),j)-p(rn(3),j));
vvar=x(k1,j)+ F*(x(k2,j)-x(k3,j));
if rand<cr
x_new(i,j)=vvar;
else
    x_new(i,j)=x(i,j);
end
    end
end% generation of x_new based on DE
% pop_in=randint(n_p,m,[1,m]);
% for i=1:n_p
% pvar=x_new([pop_in(i,:)],1:m);
% p_new(i,:)=diag(pvar)';
% end
fitx_new=[];
for i=1:n    
fitx_new(i)=bench_func(x_new(i,:),I_fno);
end% fitx_new calculation
% for i=1:n_p  
% fit_new(i)=benchmark_func(p_new(i,:),I_fno);
% end

% fit1=[];pop1=[];
% for i=1:n_p
%     if fit_new(i)>fitb(i)
%      pop1(i,:)=p(i,:);
%      fit1(i)=fitb(i);
%     else
%         pop1(i,:)=p_new(i,:);
%      fit1(i)=fit_new(i);
%     end
% end

kk=fitx>fitx_new;
kkc=fitx<fitx_new;
fitx=fitx.*kkc+fitx_new.*kk;
kk=repmat(kk',1,m);
kkc=repmat(kkc',1,m);
x=x.*kkc+x_new.*kk;



itr
an=min(fitx)
end







