% qc=[0.95,0.955,0.96,1,1.5,2,3.5];
% rc=[0,0.1,0.2,0.3,0.57,0.8,1];
% p=polyfit(rc,qc,3);
% q2=polyval(p,r);


x=linspace(0,1,100);ymin=0.95;x0=0.5;d0=0.8;f=1;
y=(ymin+0.2*(1+(x/x0).^2).^2).*(1+f*exp(-((x-x0)/d0).^2));
plot(x,y)

hold on