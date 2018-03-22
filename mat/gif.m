for j=50:150;  
        im=imread(sprintf('%d.png',(j-1)));  
        [I,map]=rgb2ind(im,180);   
        if j==50  
           imwrite(I,map,'psi2.gif','gif', 'Loopcount',inf,'DelayTime',0.9) 
       else  
           imwrite(I,map,'psi2.gif','gif','WriteMode','append','DelayTime',0.2);  
       end  
    end  
    close all;  