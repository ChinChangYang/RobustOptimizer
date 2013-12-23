function[]=display_plot(gbest,phi_desired,distance)

dim=length(gbest);
phi=linspace(0,360,5000);
yax(1)=array_factorcir(gbest,(pi/180)*phi(1),phi_desired,distance,dim);
maxi=yax(1);
for i=2:5000%This loop finds out the maximum gain 
    yax(i)=array_factorcir(gbest,(pi/180)*phi(i),phi_desired,distance,dim);
    if maxi<yax(i)
        maxi=yax(i);
    end;
end;
for i=1:5000%This loop normalizes the Y-axis and finds the normalized gain values in decibels 
    yax(i)=yax(i)/maxi;
    yax(i)=20*log10(yax(i));
end;
plot(phi,yax,'g')
xlabel('Azimuth angle(deg)');
ylabel('Gain(db)');