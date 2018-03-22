a=load('./enrdps.dat');
semilogy(a(:,1),a(:,3),'-.');

semilogy(a(:,1),a(:,6),'-.','linewidth',1.3);   
hold on
semilogy(a(:,1),a(:,13),'-.','linewidth',1.3);
xlabel('t/\tau_a')
ylabel('E_m_a_g')
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')
% plot([1,1],[-1e4,1e4])e
% set(gca,'xlim',[0 1]);
% plotyy(x,[y1;y2],x,[y3;y4]) 