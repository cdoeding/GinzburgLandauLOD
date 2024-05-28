function A = assembleLocalStiffnessMatrix(T_h,P0,l)
%ASSEMBLELOCALSTIFFNESSMATRIX Summary of this function goes here
%   Detailed explanation goes here

active_fine_triangles = find(logical(P0(:,l))); % indicies of fine triangles in current coarse triangle
T = T_h.t(active_fine_triangles,:); % fine trianlges in current coarse triangle
Nd = T_h.p; % fine nodes of fine trianlges
N_h = size(T_h.p,1); % number of fine nodes

a_i = zeros(9*size(T,1),1);
a_j = zeros(9*size(T,1),1);
a = zeros(9*size(T,1),1);
ind = 1;

grad1 = [-1; -1];
grad2 = [1; 0];
grad3 = [0; 1];

for i = 1:size(T,1)
    tri = T(i,:); %node index of triangle
    z1 = Nd(T(i,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(i,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(i,3),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    BTinv = inv(BT);
    D = BTinv*(BTinv');
    
    %% assemble
    %phi1 phi2
    e = (0.5*detBT)*grad1'*D*grad2;

    a_i(ind) = tri(1);
    a_j(ind) = tri(2);
    a(ind) = e;
    ind = ind + 1;
    
    a_i(ind) = tri(2);
    a_j(ind) = tri(1);
    a(ind) = e;
    ind = ind + 1;
    
    % phi2 phi3
    e = (0.5*detBT)*grad2'*D*grad3;

    a_i(ind) = tri(2);
    a_j(ind) = tri(3);
    a(ind) = e;
    ind = ind + 1;
    
    a_i(ind) = tri(3);
    a_j(ind) = tri(2);
    a(ind) = e;
    ind = ind + 1;
    
    % phi3 phi1
    e = (0.5*detBT)*grad3'*D*grad1;

    a_i(ind) = tri(3);
    a_j(ind) = tri(1);
    a(ind) = e;
    ind = ind + 1;
    
    a_i(ind) = tri(1);
    a_j(ind) = tri(3);
    a(ind) = e;
    ind = ind + 1;
    
    % phi1 phi1
    e = (0.5*detBT)*grad1'*D*grad1;

    a_i(ind) = tri(1);
    a_j(ind) = tri(1);
    a(ind) = e;
    ind = ind + 1;
    
    % phi2 phi2
    e = (0.5*detBT)*grad2'*D*grad2;

    a_i(ind) = tri(2);
    a_j(ind) = tri(2);
    a(ind) = e;
    ind = ind + 1;
    
    % phi3 phi3
    e = (0.5*detBT)*grad3'*D*grad3;
    
    a_i(ind) = tri(3);
    a_j(ind) = tri(3);
    a(ind) = e;
    ind = ind + 1;
    
end

A = sparse(a_i,a_j,a,N_h,N_h);

end


