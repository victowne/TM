Imax=499;
Lmax=3;
L=1;
rs=310;
drs=rs+L*(Imax+1);
tlast=800;
tstep=50;
ttotal=tlast*tstep;
t=linspace(0,ttotal,tlast+1);
p=load('Profile0.dat');
c=abs(p(rs,1)/(p(rs,2)*p(rs,5)));
% % % % % % % % % % % % % % % --Width-- % % % % % % % % % % % % % % % % %
parpool(4)
parfor i=1000:tlast+1000
                A=load(['xy',num2str(i),'.dat']);
                ww(i-999)=sqrt(sqrt(A(drs,5)^2+A(drs,6)^2)*c)*4;
end


% % % % % % % % Mode Frequence % % % % % % % %
% parpool(4)
parfor i=1000:tlast+1000
A=load(['xy',num2str(i),'.dat']);
a=A(:,5);
b=A(:,6);                                     
a=reshape(a,500,Lmax+1);               
b=reshape(b,500,Lmax+1);                             
a=a+1i*b;
phase(i-999)=angle(a(rs,2));
end
delete(gcp)
for i=1:tlast
    if (phase(i+1)-phase(i)) < 4.5
        omega(i)=(phase(i+1)-phase(i))/tstep;
    else
        omega(i)=(phase(i+1)-2*pi-phase(i))/tstep;
    end
end
omega(tlast+1)=omega(tlast);


subplot(2,1,1)
plot(t,ww,'-.','LineWidth',1.3)
title('fb0.30')
xlabel('t/\tau_a')
ylabel('w/a')
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')
subplot(2,1,2)
plot(t,-omega,'r-','Linewidth',1.)
xlabel('t/\tau_a')
ylabel('\omega')
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')