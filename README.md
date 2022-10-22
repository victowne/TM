# TM
=-------------------------------------------=  
 Author: Shuangshuang Lu  
 Maintainer: Weikang Tang (Vic Towne)   
 Mail: wtang@ipp.mpg.de  
=-------------------------------------------=   
This software was developed intended to simulate the neoclassical tearing mode instability in slab geometry based on the full-MHD equation. It was first developed by my girl at the firt two years during pursing her PhD degree. The code is not very complicated and with a good readability, so it is a great resource for the beginers to learn the tearing mode instability and practice MHD simulating. I used to use this code to foster my skill of parallel programing , like Openmp, MPI, Openacc, etc. The GPU version of this code is still far away from completion, because I don't have enough time and energy now. After I understand the gist of Openacc, I just roughly use the simplest constructs 'acc kernels' to accelerate the code. Many variabls should 'copyin' to the GPU device before the start of the main loop. because the frequent communication between the host(CPU) and device(GPU) will decrease the efficiency of parallel computing to a great extent. Perhaps, I will finish this later when I have time or never.  
