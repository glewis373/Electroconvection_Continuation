# Electroconvection_Continuation
Code for the Newton-Krylov Continuation of Amplitude Modulated Waves in Electroconvection

Contributors:  J. Jabbour, G. Lewis, M. Pugh, P.-C. Tsai 

This repository contains the Matlab code for the continuation of amplitude-modulated waves as presented in the manuscript 'Newton-Krylov continuation of amplitude-modulated waves in sheared annular electroconvection' (submitted to Physical Review E).  
To run the code, copy the files to your directory of choice, making sure the subfolder file structure is maintained, execute the 'init' file, and then, to run the code itself, execute the 'Continuation_Caller' file.

'init.m' - executing this adds the required paths for accessing the code in the subdirectories  
'Continuation_Caller.m' - contains the code for the initialization of the continuation and for stepping along the solution branch; this script calls the function 'AmpMod_NewKry' which performs the Newton iterations, etc.
'AmpMod_NewKry.m' - contains the code for the Newton iterations and computation of the preconditioning matrix
The 'SimCode' folder contains the time-stepping code for the simulation of sheared annular electroconvection; the function 'TS_3x_trunc.m' can be called to run a simulation 
The 'mat_files' folder contains the data files used to initialize the computation
The 'Utilities' folder contains code for functions used in the Newton iterations, etc. 
