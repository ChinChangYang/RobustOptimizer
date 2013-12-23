% function f=hsdrf(pop11,pop22)
function f=hsdrf(pop)

ipp=pop(end,1);
iqq=pop(end,2);
pop11=pop(1:ipp,:);
pop22=pop(end-iqq-1:end-1,:);
[p11 w]=size(pop11);
[p22 w1]=size(pop22);
A=0;B=0;
for i=1:p11
    edis=squareform(pdist([pop11(i,:);pop22]));   
    ee=min(edis(2:end,1));
    A=A+ee;
end

for i=1:p22
    edis=squareform(pdist([pop22(i,:);pop11]));   
    ee=min(edis(2:end,1));
    B=B+ee;
end    
f=max(A,B);