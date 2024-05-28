function S = assembleGlobalBilinearForm(A,T,Nd,kappa)
%ASSEMBLEGLOBALBILINEARFORM Summary of this function goes here
%   Detailed explanation goes here

s_i = zeros(9*size(T,1),1);
s_j = zeros(9*size(T,1),1);
s_val = zeros(9*size(T,1),1);
ind = 1;

grad1 = [-1; -1];
grad2 = [1; 0];
grad3 = [0; 1];

phi1 = @(x) -x(1) - x(2) + 1;
phi2 = @(x) x(1);
phi3 = @(x) x(2);

for i = 1:size(T,1)
    tri = T(i,:); %node index of triangle
    z1 = Nd(T(i,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(i,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(i,3),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    b = [z1(1); z1(2)];
    
    BTinv = inv(BT)';
    
    %% assemble
    %grad phi1 grad phi2
    e = detBT*integrate_unit_triangle(@(x) (1i/kappa * BTinv*grad1 + A(BT*x+b)*phi1(x))' * (1i/kappa * BTinv*grad2 + A(BT*x+b)*phi2(x)) ,9);
    
    s_i(ind) = tri(1);
    s_j(ind) = tri(2);
    s_val(ind) = conj(e);
    ind = ind + 1;
    
    s_i(ind) = tri(2);
    s_j(ind) = tri(1);
    s_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi2 grad phi3
    e = detBT*integrate_unit_triangle(@(x) (1i/kappa * BTinv*grad2 + A(BT*x+b)*phi2(x))' * (1i/kappa * BTinv*grad3 + A(BT*x+b)*phi3(x)) ,9);
    
    s_i(ind) = tri(2);
    s_j(ind) = tri(3);
    s_val(ind) = conj(e);
    ind = ind + 1;
    
    s_i(ind) = tri(3);
    s_j(ind) = tri(2);
    s_val(ind) = e;
    ind = ind + 1;
    
    %grad phi3 grad phi1
    e = detBT*integrate_unit_triangle(@(x) (1i/kappa * BTinv*grad3 + A(BT*x+b)*phi3(x))' * (1i/kappa * BTinv*grad1 + A(BT*x+b)*phi1(x)) ,9);
    
    s_i(ind) = tri(3);
    s_j(ind) = tri(1);
    s_val(ind) = conj(e);
    ind = ind + 1;
    
    s_i(ind) = tri(1);
    s_j(ind) = tri(3);
    s_val(ind) = e;
    ind = ind + 1;
    
    %grad phi1 grad phi1
    e = detBT*integrate_unit_triangle(@(x) (1i/kappa * BTinv*grad1 + A(BT*x+b)*phi1(x))' * (1i/kappa * BTinv*grad1 + A(BT*x+b)*phi1(x)) ,9);
    
    s_i(ind) = tri(1);
    s_j(ind) = tri(1);
    s_val(ind) = e;
    ind = ind + 1;
    
    %grad phi2 grad phi2
    e = detBT*integrate_unit_triangle(@(x) (1i/kappa * BTinv*grad2 + A(BT*x+b)*phi2(x))' * (1i/kappa * BTinv*grad2 + A(BT*x+b)*phi2(x)) ,9);
    
    s_i(ind) = tri(2);
    s_j(ind) = tri(2);
    s_val(ind) = e;
    ind = ind + 1;
    
    %grad phi3 grad phi3
    e = detBT*integrate_unit_triangle(@(x) (1i/kappa * BTinv*grad3 + A(BT*x+b)*phi3(x))' * (1i/kappa * BTinv*grad3 + A(BT*x+b)*phi3(x)) ,9);
    
    s_i(ind) = tri(3);
    s_j(ind) = tri(3);
    s_val(ind) = e;
    ind = ind + 1;
end

S = sparse(s_i,s_j,s_val);

end

