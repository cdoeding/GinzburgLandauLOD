% Ginzburg-Landau equation
% -(nabla - iA)^2 u + kappa^2*(1-|u|^2)*u = 0 in Omega = (-x_a,x_b)^2
% (nabla - iA)u * n = 0 on dOmega
% compute minimizer with FEM

clearvars
close all

%% model parameter
x_a = 0;
x_b = 1;
area = 1;
boundary = 'Neumann';

kappa = 20;
A = @(x) sqrt(2)*[sin(pi*x(1)).*cos(pi*x(2)); -sin(pi*x(2)).*cos(pi*x(1))];

%% numerical parameter
h_level = 10;

% gradient flow iteration
tau = 1;
tol = 10^(-12);
i_max = 5000;

%% coarse, fine mesh and patches
[~,T_h,~,~] = getCoarseFineTriangulation(x_a,x_b,1,h_level);
T = T_h.t;
Nd = T_h.p;
N_h = size(T_h.p,1);

%% global mass and stiffness matricies
disp("start assemble global standard matricies")
S_h = assembleGlobalBilinearForm(A,T_h.t,T_h.p,kappa);
M_h = assembleGlobalMassMatrix(T_h.t,T_h.p);

u_h = (0.8 + 0.6i)*ones(size(S_h,1),1);

F_h = assembleGlobalNonlinearMatrix(u_h,T_h.t,T_h.p);
E = (0.5/kappa^2)*real(u_h'*S_h*u_h + (kappa^2/2)*area - kappa^2*u_h'*M_h*u_h + (kappa^2/2)*u_h'*F_h*u_h);

delta = tol;
counter = 0;

% iteration
while delta >= tol && counter < i_max
    E_old = E;
    
    Mat = (1 - tau)*M_h + tau*S_h + tau*F_h;
    
    u_h = Mat\(M_h*u_h);
    F_h = assembleGlobalNonlinearMatrix(u_h,T_h.t,T_h.p);
    E = (0.5/kappa^2)*real(u_h'*S_h*u_h + (kappa^2/2)*area - kappa^2*u_h'*M_h*u_h + (kappa^2/2)*u_h'*F_h*u_h);
        
    delta = abs(E_old - E);
    counter = counter + 1;
    
    disp(counter)
    disp(E)
    disp(delta)
    
end

%% save results
name = 'groundstate';
save(strcat(name,'.mat'),'-v7.3')




