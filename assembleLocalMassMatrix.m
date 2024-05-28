function M = assembleLocalMassMatrix(T_h,P0,l)
%ASSEMBLELOCALBILINEARFORM Summary of this function goes here
%   Detailed explanation goes here

active_fine_triangles = find(logical(P0(:,l))); % indicies of fine triangles in current coarse triangle
T = T_h.t(active_fine_triangles,:); % fine trianlges in current coarse triangle
Nd = T_h.p; % fine nodes of fine trianlges
N_h = size(T_h.p,1); % number of fine nodes

m_i = zeros(9*size(T,1),1);
m_j = zeros(9*size(T,1),1);
m = zeros(9*size(T,1),1);
ind = 1;

phi1 = @(x,y) -x - y + 1;
phi2 = @(x,y) x;
phi3 = @(x,y) y;

for i = 1:size(T,1)
    tri = T(i,:); %node index of triangle
    z1 = Nd(T(i,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(i,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(i,3),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    
    %% assemble
    %phi1 phi2
    e = ((0.5*detBT)/3)*( phi1(0.5,0)*phi2(0.5,0) + phi1(0,0.5)*phi2(0,0.5) + phi1(0.5,0.5)*phi2(0.5,0.5) );
    
    m_i(ind) = tri(1);
    m_j(ind) = tri(2);
    m(ind) = e;
    ind = ind + 1;
    
    m_i(ind) = tri(2);
    m_j(ind) = tri(1);
    m(ind) = e;
    ind = ind + 1;
    
    % phi2 phi3
    e = ((0.5*detBT)/3)*( phi2(0.5,0)*phi3(0.5,0) + phi2(0,0.5)*phi3(0,0.5) + phi2(0.5,0.5)*phi3(0.5,0.5) );
    
    m_i(ind) = tri(2);
    m_j(ind) = tri(3);
    m(ind) = e;
    ind = ind + 1;
    
    m_i(ind) = tri(3);
    m_j(ind) = tri(2);
    m(ind) = e;
    ind = ind + 1;
    
    % phi3 phi1
    e = ((0.5*detBT)/3)*( phi3(0.5,0)*phi1(0.5,0) + phi3(0,0.5)*phi1(0,0.5) + phi3(0.5,0.5)*phi1(0.5,0.5) );
    
    m_i(ind) = tri(3);
    m_j(ind) = tri(1);
    m(ind) = e;
    ind = ind + 1;
    
    m_i(ind) = tri(1);
    m_j(ind) = tri(3);
    m(ind) = e;
    ind = ind + 1;
    
    % phi1 phi1
    e = ((0.5*detBT)/3)*( phi1(0.5,0)*phi1(0.5,0) + phi1(0,0.5)*phi1(0,0.5) + phi1(0.5,0.5)*phi1(0.5,0.5) );
    
    m_i(ind) = tri(1);
    m_j(ind) = tri(1);
    m(ind) = e;
    ind = ind + 1;
    
    % phi2 phi2
    e = ((0.5*detBT)/3)*( phi2(0.5,0)*phi2(0.5,0) + phi2(0,0.5)*phi2(0,0.5) + phi2(0.5,0.5)*phi2(0.5,0.5) );
    
    m_i(ind) = tri(2);
    m_j(ind) = tri(2);
    m(ind) = e;
    ind = ind + 1;
    
    % phi3 phi3
    e = ((0.5*detBT)/3)*( phi3(0.5,0)*phi3(0.5,0) + phi3(0,0.5)*phi3(0,0.5) + phi3(0.5,0.5)*phi3(0.5,0.5) );
    
    m_i(ind) = tri(3);
    m_j(ind) = tri(3);
    m(ind) = e;
    ind = ind + 1;
    
end

M = sparse(m_i,m_j,m,N_h,N_h);

end

