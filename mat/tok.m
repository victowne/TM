% 绘制3D的托卡马克位形图， (V 0.1 by Jiale Chan for Y. H. Huang)
% Dee Formula
% 特征参数
    rzero = 2.0;
    rmax = 0.75;
    eshape = 2.0;
    xshape = 0.5;
% 极向角与环向角
    theta = linspace(0, 2*pi, 50);
    phi = linspace(0, 1.3*pi, 50);
% 将数据矩阵化
    [phi_g, theta_g ] = meshgrid(phi, theta);
% 参考TOQ网页给的简单公式生成磁面数据：https://fusion.gat.com/THEORY/toq/geometry.html
    z = eshape*rmax*sin(theta_g);
    R = rzero + rmax*( cos(theta_g)-xshape*sin(theta_g).^2 ) ;
    x = R.*cos(phi_g);
    y = R.*sin(phi_g);
% 绘制3D磁面
    surf(x,y,z); shading interp; %或者 surf(x,y,z, 'EdgeColor','none');
%     hold on;
%随便画的一条线
%      R_line = rzero + rmax*( cos(theta)-xshape*sin(theta).^2 ) ;
%      x_line = R_line.*cos(phi*1.5);
%      y_line = R_line.*sin(phi*1.5);
%      z_line = eshape*rmax*sin(theta);
%      plot3(x_line, y_line, z_line, 'k'); 
 
% 设定图的一些其他参数
%     colorbar; %绘制色条
    xlabel('x');   ylabel('y');   zlabel('z');%坐标轴的命名
    axis equal; %将坐标轴比例设为“相等”
    view(2,40); %改变视角