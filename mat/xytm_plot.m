                    % ------------------ %
                    % TM_contour_plot    %   
                    %           twk/2016 %
                    % -------------------%
%---argument---%
irmax = 512;
mode = 5;
lx = 1;
ly = 2;
dy = 0.01;
ny = ly/dy+1;
%--------------%
format shortE
A = load('./xy087.dat');
a = A(:,5);
b = A(:,6);           %----------  chose value(re,im) to plot            
A = reshape(a,irmax+1,mode);               
B = reshape(b,irmax+1,mode); %---------- matrix reform 
A = A+1i*B;                     
k = ((0:mode-1)*pi)';   %----------- ky
x = -lx:2*lx/irmax:lx;  %----------- xgrid
y = 0:dy:ly;   %----------- ygrid                
B = 2*exp(1i*k*y);
B(1,:) = B(1,:)/2;
PSI = real(A*B);  %-------------Fourier Transform                     
PSI0 = repmat((-x+(1.0+0.233509)/2.68298*atan(sinh(2.68298*x)))',1,ny);    %equilibrium
PSI = PSI+PSI0;
contourf(y,x,PSI,100)
xlabel('Y')
ylabel('X')
%view(90,-90)                %reverse axis
%pcolor(y,x,PSI);shading interp