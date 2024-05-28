# GinzburgLandauLOD
Code used in the paper "Vortex-capturing multiscale spaces for the Ginzburg-Landau equation" 

Results are calculated using the following procedures:

1. Calculate GL minimizers using the LOD method:  
    - run main_GL_LOD.m
    - set/change parameters: kappa, A, beta, l, H_level (coarse mesh size = 2^{-H_level}), h_level (fine mesh size = 2^{-h_level})
    - Results are saved to file groundstate.mat (optional: change filename)
    - plot results by running plot_groundstate.m
    - option: In lines 54, 65 and 79 you can change for parallel execution by using methods getCorrectorMatrixParallel.m and assembleGlobalNonlinearMatrixParallel.m.

2. Calculate GL minimizers with the default P1 FEM:  
    - Execute main_GL_FEM.m
    - set/change parameters: kappa, A and h_level (mesh size = 2^{-h_level})
    - save results to file groundstate.mat (optional: change filename)
    - plot results by running plot_groundstate.m
   
3. Compute FEM best approximation:  
    - First run and save the desired FEM solution (step 2)
    - Run main_bestapproximation.m
