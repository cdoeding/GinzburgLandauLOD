function F = assembleGlobalNonlinearMatrix(u,T,Nd)
%ASSEMBLEGLOBALNONLINEARMATRIX Summary of this function goes here
%   Detailed explanation goes here

F_i = zeros(9*size(T,1),1);
F_j = zeros(9*size(T,1),1);
F_val = zeros(9*size(T,1),1);
ind = 1;

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
    
    %% construct f = |u|^2 on triangle
    f1 = @(x) u(tri(1))*phi1(x); % contribution of phi1
    f2 = @(x) u(tri(2))*phi2(x); % contribution of phi2
    f3 = @(x) u(tri(3))*phi3(x); % contribution of phi3
    
    f = @(x) abs(f1(x) + f2(x) + f3(x))^2;
    
    %% assemble
    %grad phi1 grad phi2
    e = detBT*integrate_unit_triangle(@(x) f(x)*phi1(x)*phi2(x),6);
    
    F_i(ind) = tri(1);
    F_j(ind) = tri(2);
    F_val(ind) = e;
    ind = ind + 1;
    
    F_i(ind) = tri(2);
    F_j(ind) = tri(1);
    F_val(ind) = e;
    ind = ind + 1;
    
    %grad phi2 grad phi3
    e = detBT*integrate_unit_triangle(@(x) f(x)*phi2(x)*phi3(x),6);
    
    F_i(ind) = tri(2);
    F_j(ind) = tri(3);
    F_val(ind) = e;
    ind = ind + 1;
    
    F_i(ind) = tri(3);
    F_j(ind) = tri(2);
    F_val(ind) = e;
    ind = ind + 1;
    
    %grad phi3 grad phi1
    e = detBT*integrate_unit_triangle(@(x) f(x)*phi3(x)*phi1(x),6);
    
    F_i(ind) = tri(3);
    F_j(ind) = tri(1);
    F_val(ind) = e;
    ind = ind + 1;
    
    F_i(ind) = tri(1);
    F_j(ind) = tri(3);
    F_val(ind) = e;
    ind = ind + 1;
    
    %grad phi1 grad phi1
    e = detBT*integrate_unit_triangle(@(x) f(x)*phi1(x)*phi1(x),6);
    
    F_i(ind) = tri(1);
    F_j(ind) = tri(1);
    F_val(ind) = e;
    ind = ind + 1;
    
    %grad phi2 grad phi2
    e = detBT*integrate_unit_triangle(@(x) f(x)*phi2(x)*phi2(x),6);
    
    F_i(ind) = tri(2);
    F_j(ind) = tri(2);
    F_val(ind) = e;
    ind = ind + 1;
    
    %grad phi3 grad phi3
    e = detBT*integrate_unit_triangle(@(x) f(x)*phi3(x)*phi3(x),6);
    
    F_i(ind) = tri(3);
    F_j(ind) = tri(3);
    F_val(ind) = e;
    ind = ind + 1;
    
end

F = sparse(F_i,F_j,F_val);

end

