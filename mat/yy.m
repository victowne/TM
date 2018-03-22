[AX,H1,H2] = plotyy(t,[ww21;ww32],t,-omega);
set(get(AX(1),'Ylabel'),'String','w/a')
set(get(AX(2),'Ylabel'),'String','\omega')
set(H1(1),'linewidth',1,'linestyle','-.')
set(H1(2),'linewidth',1,'linestyle','-.')
set(H2,'linewidth',1,'linestyle','-')
set(AX(1),'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')
set(AX(2),'FontName','Times New Roman','FontSize',12,'linewidth',1,'Fontweight','bold')
set(gca,'position',[0.1 0.11 0.775 0.815])


% --------------------------rmp-------------------------------- %
% for i=1:10000
%     rmpf(i)=0;
% end
% for i=10001:19800
%     rmpf(i)=1e-7*(t(i)-10000);
% end
% for i=19801:80000
%     rmpf(i)=9.8e-4;
% end
% for i=177:226
%     rmp(i)=0.05-(0.05/10000)*(t(i)-35000);
% end
% for i=227:301
%     rmp(i)=0;
% end
% hold on
% plot(t,rmp,'--','linewidth',1,'color',[0.49 0.18 0.56])