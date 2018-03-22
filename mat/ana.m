% ---------------------------  Width ------------------------------%
s=1;k=1;tlast=1000;tstep=10;ttotal=tlast*tstep;rs=287;
t=linspace(0,ttotal,tlast+1);
p=load('Profile0.dat');c=abs(p(rs,1)/(p(rs,2)*p(rs,5)));
for i=0:tlast
                A=load(['xy',num2str(i),'.dat']);
                ww(s)=sqrt(sqrt(A(rs+500,5)^2+A(rs+500,6)^2)*c)*4;
                s=s+1;
end
plot(t,ww,'-.','color',[0 0.5 0],'LineWidth',1.3)
% t1=linspace(0,200000,200000);
% www=interp1(t,ww,t1,'spline');
% plot(t1,www,'-.','color',[0 0.5 0],'LineWidth',1.3)
xlabel('t/\tau_a')
ylabel('w/a')
set(gca,'FontName','Times New Roman','FontSize',13,'linewidth',1,'Fontweight','bold')
% [AX,H1,H2] = plotyy(t,ww,t,phase/pi);
% set(get(AX(1),'Ylabel'),'String','w/a')
% set(get(AX(2),'Ylabel'),'String','\Phi/\pi') 
figure
%----------------------- Mode Frequence ---------------------------%
for i=0:tlast
A=load(['xy',num2str(i),'.dat']);
r=p(:,1);                                       
r=r';
eq=p(:,4);                                       
a=A(:,5);
b=A(:,6);                                     
a=reshape(a,500,4);               
b=reshape(b,500,4);                             
a=a+1i*b;
% plot(a(201,2),'+');
% hold on
phase(k)=angle(a(rs,2));
k=k+1;
end
plot(t,phase/pi)
xlabel('t/\tau_a')
ylabel('\Phi/\pi')
set(gca,'FontSize',13,'Fontname','Arial')
if k<0
t1=linspace(0,20000,20000);y=complex(cos(-0.01848*t1),sin(-0.01848*t1));y=angle(y);hold on;plot(t1,y/pi,'--')
t1=linspace(5000,6000,1000);y=complex(cos(-0.01848*t1+1.5708),sin(-0.01848*t1+1.5708));y=angle(y);hold on;plot(t1,y/pi,'--')
figure

%---------------------- gif plot --------------------------------%
s=0;kn=3
for  i=0:tlast
A=load(['xy',num2str(i),'.dat']);
r=p(:,1);                                       %径向范围
r=r';
eq=p(:,4);                                       %psi0
a=A(:,5);
b=A(:,6);                                       %选择画图的值  
a=reshape(a,500,kn+1);               
b=reshape(b,500,kn+1);                             %矩阵变形
a=a+1i*b;                     
n=0:kn;                                          %模数
k=3*n;                                          %kn
k=k';
theta=0:0.01:2.005*pi;                              %极向范围
b=2*exp(1i*k*theta);
b(1,:)=b(1,:)/2;
PSI=a*b;                      
PSI=real(PSI);
PSI0=eq'+0.25*r.*r/3;                            %平衡量
PSI0=PSI0';
PSI0=repmat(PSI0,1,630);      
PSI=PSI+PSI0;
r=r(:,1:270);PSI=PSI(1:270,:);             %截取范围
[tt, rr] = meshgrid(theta, r);
[x, y] = pol2cart(tt, rr);
contourf(y,x,PSI,50,'linecolor','none');
%contourf(x,y,PSI,50);
title(sprintf(['t=' num2str(s)]));
saveas(gcf,['t=' num2str(s) '.jpg']);
s=s+tstep;
end

for j=1:9999;  
        im=imread(sprintf('t=%d.jpg',(j-1)*tstep));  
        [I,map]=rgb2ind(im,180);   
        if j==1  
           imwrite(I,map,'psi.gif','gif', 'Loopcount',inf,'DelayTime',0.2) 
       else  
           imwrite(I,map,'psi.gif','gif','WriteMode','append','DelayTime',0.2);  
       end  
    end  
    close all;
end