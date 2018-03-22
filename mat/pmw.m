Imax=499;
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
parfor i=0:tlast
                A=load(['xy',num2str(i),'.dat']);
                ww21(i+1)=sqrt(sqrt(A(drs,5)^2+A(drs,6)^2)*c)*4;
end



Imax=499;
L=2;
rs=246;
drs=rs+L*(Imax+1);
% tlast=150;
% tstep=200;
% ttotal=tlast*tstep;
% t=linspace(0,ttotal,tlast+1);
p=load('Profile0.dat');
c=abs(p(rs,1)/(p(rs,2)*p(rs,5)));
% % % % % % % % % % % % % % % --Width-- % % % % % % % % % % % % % % % % %
parfor i=0:tlast
                A=load(['xy',num2str(i),'.dat']);
                ww32(i+1)=sqrt(sqrt(A(drs,5)^2+A(drs,6)^2)*c)*4;               
end
% delete(gcp)

%  plot(t,ww21,'-.','LineWidth',1.3)
%  hold on
%  plot(t,ww32,'-.','LineWidth',1.3)
%  xlabel('t/\tau_a')
%  ylabel('w/a')
%  set(gca,'FontName','Times New Roman','FontSize',13,'linewidth',1,'Fontweight','bold')