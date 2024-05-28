clearvars
close all

%% model parameter
x_a = 0;
x_b = 1;
area = 1;
boundary = 'Neumann';

A = @(x) sqrt(2)*[sin(pi*x(1)).*cos(pi*x(2)); -sin(pi*x(2)).*cos(pi*x(1))];

%% numerical parameter
% LOD
H_level = 6;
h_level = 10;
ell = 10;
kappa = 20;
beta = 0;
name = 'groundstate'; %FILE NAME HERE

%% load reference solution
u_ref = load(strcat(name,'.mat')).u_h;

%% coarse, fine mesh and patches
[T_H,T_h,P1,P0] = getCoarseFineTriangulation(x_a,x_b,H_level,h_level);
P = P1';
model = createpde();
Geo = geometryFromMesh(model,T_h.p',T_h.t');

B_H = getBoundaryNodes(T_H.p,x_a,x_b,boundary);
B_h = getBoundaryNodes(T_h.p,x_a,x_b,boundary);
Bd_H = getBoundaryRestriction(B_H);

patches = getPatches(T_H,ell); % patches_ij non-zero iff jth triangle is in patch of ith triangle

N_h = size(T_h.p,1);
N_H = size(T_H.p,1);
NT_H = size(T_H.t,1);

%% global mass and stiffness matricies
disp("start assemble global standard matricies")
A_h = assembleGlobalStiffnessMatrix(T_h.t,T_h.p); % carefullly: is the 'normal one'
S_h = assembleGlobalBilinearForm(A,T_h.t,T_h.p,kappa);
M_h = assembleGlobalMassMatrix(T_h.t,T_h.p);

%% compute LOD basis (corrector)
disp("start computing corrector")
tic;
Q = getCorrectorMatrix(T_H,T_h,patches,A,kappa,beta,S_h,M_h,P1,P0,B_H,B_h); %can change in getCorrectorMatrixParallel for parallel execution
time_Q = toc;

%% compute LOD matricies
A_LOD = Bd_H*(P + Q)*A_h*(P + Q)'*Bd_H';
M_LOD = Bd_H*(P + Q)*M_h*(P + Q)'*Bd_H';

%% compute Bestapproximation
Mat = ((1/kappa^2)*A_LOD+M_LOD)';
vek = (P+Q)*((1/kappa^2)*A_h+M_h)'*u_ref;
U_LOD = Mat\vek;

%% reconstruct Bestapproximation
u_best_h = (P + Q)'*U_LOD;

