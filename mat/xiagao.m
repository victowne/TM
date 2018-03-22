% 特征参数
    R0 = 4.0;
    a = 1.0;
    elongation = 1.0;
    triangularity = 0.0;
% 极向角与环向角
    theta = linspace(0, 2*pi, 50);
    phi = linspace(0, 1.8*pi, 50);
% 将数据矩阵化
    [phi2, theta2 ] = meshgrid(phi, theta);
% 绘制3D磁面
bb1=exp(0*i*theta2-0*i*phi2);
bb2=exp(3i*theta2-i*phi2);
bb3=exp(6i*theta2-2i*phi2);
bb4=exp(9i*theta2-3i*phi2);
for s=599:599;
    color=A(s,1)*bb1+A(s,2)*bb2+A(s,3)*bb3+A(s,4)*bb4;
    color=A(s,4)*bb4;
    color=real(color);
    a=s/600;
    z = elongation*a*sin(theta2);
    R = R0 + a*( cos(theta2)-triangularity*sin(theta2).^2 ) ;
    x = R.*cos(phi2);
    y = R.*sin(phi2);
    surf(x,y,z,color); shading interp;
    hold on;
end
axis equal