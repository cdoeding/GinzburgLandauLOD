% Ginzburg-Landau equation
% -(nabla - iA)^2 u + kappa^2*(1-|u|^2)*u = 0 in Omega = (-x_a,x_b)^2
% (nabla - iA)u * n = 0 on dOmega
% compute minimizer with LOD method

clearvars
close all

%% model parameter
x_a = 0;
x_b = 1;
area = 1;
boundary = 'Neumann';

kappa = 20;
beta = 0;
A = @(x) sqrt(2)*[sin(pi*x(1)).*cos(pi*x(2)); -sin(pi*x(2)).*cos(pi*x(1))];
u0 = @(x) 0.8 + 0.6i;

%% numerical parameter
% LOD
H_level = 4;
h_level = 6;
ell = 5;

% gradient flow iteration
tau0 = 1;
tau = tau0;
tol = 10^(-12);
i_max = 5000;

%% coarse, fine mesh and patches
[T_H,T_h,P1,P0] = getCoarseFineTriangulation(x_a,x_b,H_level,h_level);
P = P1';

B_H = getBoundaryNodes(T_H.p,x_a,x_b,boundary);
B_h = getBoundaryNodes(T_h.p,x_a,x_b,boundary);
Bd_H = getBoundaryRestriction(B_H);

patches = getPatches(T_H,ell); % patches_ij non-zero iff jth triangle is in patch of ith triangle

N_h = size(T_h.p,1);
N_H = size(T_H.p,1);
NT_H = size(T_H.t,1);

%% global mass and stiffness matricies
disp("start assemble global standard matricies")
S_h = assembleGlobalBilinearForm(A,T_h.t,T_h.p,kappa);
M_h = assembleGlobalMassMatrix(T_h.t,T_h.p);

%% compute LOD basis (corrector)
disp("start computing corrector")
tic;
Q = getCorrectorMatrix(T_H,T_h,patches,A,kappa,beta,S_h,M_h,P1,P0,B_H,B_h); %can change in getCorrectorMatrixParallel for parallel execution
time_Q = toc;

%% compute LOD matricies
S_LOD = Bd_H*(P + Q)*S_h*(P + Q)'*Bd_H';
M_LOD = Bd_H*(P + Q)*M_h*(P + Q)'*Bd_H';

%% iteration
% intialize
disp("start time iteration")
U_LOD = evaluateInitialValue(u0,T_h.t,T_h.p,M_LOD,P+Q); % get intiial value in LOD basis
F_h = assembleGlobalNonlinearMatrix((P + Q)'*U_LOD,T_h.t,T_h.p); %can change in assembleGlobalNonlinearMatrixParallel for parallel execution % assemble matrix for cubic nonlinearity on fine mesh
F_LOD = Bd_H*(P + Q)*F_h*(P + Q)'*Bd_H'; % transform to LOD space
E = (0.5/kappa^2)*real(U_LOD'*S_LOD*U_LOD + (kappa^2/2)*area - kappa^2*U_LOD'*M_LOD*U_LOD + (kappa^2/2)*U_LOD'*F_LOD*U_LOD);

delta = tol;
counter = 0;

% iteration
while delta >= tol && counter < i_max
    E_old = E;
    
    Mat = (1 - tau)*M_LOD + tau*S_LOD + tau*F_LOD;
    
    U_LOD = Mat\(M_LOD*U_LOD);
    F_h = assembleGlobalNonlinearMatrix((P + Q)'*U_LOD,T_h.t,T_h.p); %can change in assembleGlobalNonlinearMatrixParallel for parallel execution
    F_LOD = Bd_H*(P + Q)*F_h*(P + Q)'*Bd_H';
    E = (0.5/kappa^2)*real(U_LOD'*S_LOD*U_LOD + (kappa^2/2)*area - kappa^2*U_LOD'*M_LOD*U_LOD + (kappa^2/2)*U_LOD'*F_LOD*U_LOD);
    delta = abs(E_old - E);
    counter = counter + 1;
    
    disp(counter)
    disp(E)
    disp(delta)
    
    %% reconstruct solution
    u_h = (P + Q)'*U_LOD;
    
end

%% save results
name = 'groundstate';
save(strcat(name,'.mat'),'-v7.3')









