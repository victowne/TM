tlast=1600;
Lmax=3;        
tstep=50;
m=2;
p=load('./Profile0.dat');
eq=p(:,4); 
r=p(:,1);                                       
r=r';
% =============================GO=========================== %
for  i=1100:1400
A=load(['xy',num2str(i),'.dat']);
a=A(:,5);
b=A(:,6);                                      
a=reshape(a,500,Lmax+1);               
b=reshape(b,500,Lmax+1);                             
a=a+1i*b;                     
k=m*(0:Lmax);                                          
k=k';
theta=0:0.01:2.005*pi;                              
b=2*exp(1i*k*theta);
b(1,:)=b(1,:)/2;
PSI=a*b;                      
PSI=real(PSI);
PSI0=eq'+0.25*r.*r/4;                            
PSI0=PSI0';
PSI0=repmat(PSI0,1,630);      
PSI=PSI+PSI0;
% r=r(:,1:270);PSI=PSI(1:270,:);             %Ωÿ»°∑∂Œß
[tt, rr] = meshgrid(theta, r);
[x, y] = pol2cart(tt, rr);
contourf(y,x,PSI,50,'linecolor','none');
title(sprintf(['t=' num2str((i-1)*tstep)]));
saveas(gcf,['t=' num2str((i-1)*tstep) '.jpg']);
end

for j=1100:9999;  
        im=imread(sprintf('t=%d.jpg',(j-1)*tstep));  
        [I,map]=rgb2ind(im,180);   
        if j==1100  
           imwrite(I,map,'psi.gif','gif', 'Loopcount',inf,'DelayTime',0.2) 
       else  
           imwrite(I,map,'psi.gif','gif','WriteMode','append','DelayTime',0.2);  
       end  
end  
    close all; 