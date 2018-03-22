rs=310;
Lmax=3;
tlast=1600;
tstep=50;
ttotal=tlast*tstep;
t=linspace(0,ttotal,tlast+1);
% % % % % % % % Mode Frequence % % % % % % % %
% parpool(4)
parfor i=0:tlast
A=load(['xy',num2str(i),'.dat']);
a=A(:,5);
b=A(:,6);                                     
a=reshape(a,500,Lmax+1);               
b=reshape(b,500,Lmax+1);                             
a=a+1i*b;
phase(i+1)=angle(a(rs,2));
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


% plot(t,phase/pi)
% xlabel('t/\tau_a')
% ylabel('\Phi/\pi')
% set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')



% t1=linspace(0,9120,9120);
% y=complex(cos(-0.03688*t1),sin(-0.03688*t1));y=angle(y);hold on;plot(t1,y/pi)






% subplot(2,1,1)
% plot(t,ww,'-.','LineWidth',1.3)
% % title('1e5')
% xlabel('t/\tau_a')
% ylabel('w/a')
% set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')
% subplot(2,1,2)
% plot(t,-omega,'r-','Linewidth',1.)
% xlabel('t/\tau_a')
% ylabel('\omega')
% set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')