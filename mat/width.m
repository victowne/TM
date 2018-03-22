Imax=499;
L=1;
rs=310;
drs=rs+L*(Imax+1);
tlast=1600;
tstep=50;
ttotal=tlast*tstep;
t=linspace(0,ttotal,tlast+1);
p=load('Profile0.dat');
c=abs(p(rs,1)/(p(rs,2)*p(rs,5)));
% % % % % % % % % % % % % % % --Width-- % % % % % % % % % % % % % % % % %
parpool(4)
parfor i=0:tlast
                A=load(['xy',num2str(i),'.dat']);
                ww(i+1)=sqrt(sqrt(A(drs,5)^2+A(drs,6)^2)*c)*4;
end
% plot(t,ww,'-.','LineWidth',1.3)
% xlabel('t/\tau_a')
% ylabel('w/a')
% set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')

% axis([xmin xmax ymin ymax])
% [AX,H1,H2] = plotyy(t,ww,t,phase/pi);
% set(get(AX(1),'Ylabel'),'String','w/a')
% set(get(AX(2),'Ylabel'),'String','\Phi/\pi') 
% delete(gcp)