for j=1:500;  
        im=imread(sprintf('t=%d.jpg',(j-1)*10+19000));  
        [I,map]=rgb2ind(im,180);   
        if j==1  
           imwrite(I,map,'psi3.gif','gif', 'Loopcount',inf,'DelayTime',0.9) 
       else  
           imwrite(I,map,'psi3.gif','gif','WriteMode','append','DelayTime',0.2);  
       end  
    end  
    close all;  