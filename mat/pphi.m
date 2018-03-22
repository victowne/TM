A=load('./xy0050.dat');
p=load('./Profile0.dat');
r=p(:,1);                                       %径向范围
r=r';
eq=p(:,16);                                      %phi0
a=A(:,3);
b=A(:,4);                                       %选择画图的值  
a=reshape(a,500,4);              
b=reshape(b,500,4);                             %矩阵变形
a=a+1i*b;                     
n=0:3;                                          %模数
k=2*n;                                          %kn
k=k';
theta=0:0.01:2.005*pi;                              %极向范围
b=2*exp(1i*k*theta);
b(1,:)=b(1,:)/2;
PHI=a*b;                      
PHI=real(PHI);
PHI0=eq';                                        %平衡量
PHI0=PHI0';
PHI0=repmat(PHI0,1,630);      
% PHI=PHI+PHI0;
%r=r(:,1:300);PHI=PHI(1:300,:);             %截取范围
[tt, rr] = meshgrid(theta, r);
[x, y] = pol2cart(tt, rr);
contourf(y,x,PHI,50,'linecolor','none')
% contourf(x,y,PHI,50)