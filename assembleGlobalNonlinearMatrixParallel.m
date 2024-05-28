function F = assembleGlobalNonlinearMatrixParallel(u,T,Nd)
%ASSEMBLEGLOBALNONLINEARMATRIX Summary of this function goes here
%   Detailed explanation goes here

N_h = size(Nd,1);
NT_h = size(T,1);

w = [0.225;
    0.132394152788506;
    0.125939180544827 ];

alpha = [0.059715871789770;
    0.797426985353087 ];

beta  = [0.470142064105115;
    0.101286507323456 ];

z = [1/3, 1/3;
    beta(1), beta(1);
    alpha(1),beta(1);
    beta(1),alpha(1);
    beta(2), beta(2);
    alpha(2),beta(2);
    beta(2),alpha(2)];

spmd
    my_index_start = floor(NT_h*(labindex - 1)/numlabs + 1);
    my_index_end = floor(NT_h*(labindex)/numlabs);
    
    F_i = zeros(9*size(T,1),1);
    F_j = zeros(9*size(T,1),1);
    F_val = zeros(9*size(T,1),1);
    ind = 1;
    
    for i = my_index_start:my_index_end
        
        tri = T(i,:); %node index of triangle
        z1 = Nd(T(i,1),:); %coordinates of 1st triangle node
        z2 = Nd(T(i,2),:); %coordinates of 2nd triangle node
        z3 = Nd(T(i,3),:); %coordinates of 3rd triangle node
        
        %% transformation to refenrence triangle
        BT = [z2(1)-z1(1), z3(1)-z1(1); ...
            z2(2)-z1(2), z3(2)-z1(2)];
        
        detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
        
        %% assemble
        %grad phi1 grad phi2
        e = detBT*0.5*( ...
            w(1)*phi1(z(1,:)')*phi2(z(1,:)')*abs(u(tri(1))*phi1(z(1,:)') + u(tri(2))*phi2(z(1,:)') + u(tri(3))*phi3(z(1,:)'))^2 ...
            + w(2)*phi1(z(2,:)')*phi2(z(2,:)')*abs(u(tri(1))*phi1(z(2,:)') + u(tri(2))*phi2(z(2,:)') + u(tri(3))*phi3(z(2,:)'))^2 ...
            + w(2)*phi1(z(3,:)')*phi2(z(3,:)')*abs(u(tri(1))*phi1(z(3,:)') + u(tri(2))*phi2(z(3,:)') + u(tri(3))*phi3(z(3,:)'))^2 ...
            + w(2)*phi1(z(4,:)')*phi2(z(4,:)')*abs(u(tri(1))*phi1(z(4,:)') + u(tri(2))*phi2(z(4,:)') + u(tri(3))*phi3(z(4,:)'))^2 ...
            + w(3)*phi1(z(5,:)')*phi2(z(5,:)')*abs(u(tri(1))*phi1(z(5,:)') + u(tri(2))*phi2(z(5,:)') + u(tri(3))*phi3(z(5,:)'))^2 ...
            + w(3)*phi1(z(6,:)')*phi2(z(6,:)')*abs(u(tri(1))*phi1(z(6,:)') + u(tri(2))*phi2(z(6,:)') + u(tri(3))*phi3(z(6,:)'))^2 ...
            + w(3)*phi1(z(7,:)')*phi2(z(7,:)')*abs(u(tri(1))*phi1(z(7,:)') + u(tri(2))*phi2(z(7,:)') + u(tri(3))*phi3(z(7,:)'))^2 ...
            );
        
        F_i(ind) = tri(1);
        F_j(ind) = tri(2);
        F_val(ind) = e;
        ind = ind + 1;
        
        F_i(ind) = tri(2);
        F_j(ind) = tri(1);
        F_val(ind) = e;
        ind = ind + 1;
        
        %grad phi2 grad phi3
        e = detBT*0.5*( ...
            w(1)*phi2(z(1,:)')*phi3(z(1,:)')*abs(u(tri(1))*phi1(z(1,:)') + u(tri(2))*phi2(z(1,:)') + u(tri(3))*phi3(z(1,:)'))^2 ...
            + w(2)*phi2(z(2,:)')*phi3(z(2,:)')*abs(u(tri(1))*phi1(z(2,:)') + u(tri(2))*phi2(z(2,:)') + u(tri(3))*phi3(z(2,:)'))^2 ...
            + w(2)*phi2(z(3,:)')*phi3(z(3,:)')*abs(u(tri(1))*phi1(z(3,:)') + u(tri(2))*phi2(z(3,:)') + u(tri(3))*phi3(z(3,:)'))^2 ...
            + w(2)*phi2(z(4,:)')*phi3(z(4,:)')*abs(u(tri(1))*phi1(z(4,:)') + u(tri(2))*phi2(z(4,:)') + u(tri(3))*phi3(z(4,:)'))^2 ...
            + w(3)*phi2(z(5,:)')*phi3(z(5,:)')*abs(u(tri(1))*phi1(z(5,:)') + u(tri(2))*phi2(z(5,:)') + u(tri(3))*phi3(z(5,:)'))^2 ...
            + w(3)*phi2(z(6,:)')*phi3(z(6,:)')*abs(u(tri(1))*phi1(z(6,:)') + u(tri(2))*phi2(z(6,:)') + u(tri(3))*phi3(z(6,:)'))^2 ...
            + w(3)*phi2(z(7,:)')*phi3(z(7,:)')*abs(u(tri(1))*phi1(z(7,:)') + u(tri(2))*phi2(z(7,:)') + u(tri(3))*phi3(z(7,:)'))^2 ...
            );
        
        F_i(ind) = tri(2);
        F_j(ind) = tri(3);
        F_val(ind) = e;
        ind = ind + 1;
        
        F_i(ind) = tri(3);
        F_j(ind) = tri(2);
        F_val(ind) = e;
        ind = ind + 1;
        
        %grad phi3 grad phi1
        e = detBT*0.5*( ...
            w(1)*phi3(z(1,:)')*phi1(z(1,:)')*abs(u(tri(1))*phi1(z(1,:)') + u(tri(2))*phi2(z(1,:)') + u(tri(3))*phi3(z(1,:)'))^2 ...
            + w(2)*phi3(z(2,:)')*phi1(z(2,:)')*abs(u(tri(1))*phi1(z(2,:)') + u(tri(2))*phi2(z(2,:)') + u(tri(3))*phi3(z(2,:)'))^2 ...
            + w(2)*phi3(z(3,:)')*phi1(z(3,:)')*abs(u(tri(1))*phi1(z(3,:)') + u(tri(2))*phi2(z(3,:)') + u(tri(3))*phi3(z(3,:)'))^2 ...
            + w(2)*phi3(z(4,:)')*phi1(z(4,:)')*abs(u(tri(1))*phi1(z(4,:)') + u(tri(2))*phi2(z(4,:)') + u(tri(3))*phi3(z(4,:)'))^2 ...
            + w(3)*phi3(z(5,:)')*phi1(z(5,:)')*abs(u(tri(1))*phi1(z(5,:)') + u(tri(2))*phi2(z(5,:)') + u(tri(3))*phi3(z(5,:)'))^2 ...
            + w(3)*phi3(z(6,:)')*phi1(z(6,:)')*abs(u(tri(1))*phi1(z(6,:)') + u(tri(2))*phi2(z(6,:)') + u(tri(3))*phi3(z(6,:)'))^2 ...
            + w(3)*phi3(z(7,:)')*phi1(z(7,:)')*abs(u(tri(1))*phi1(z(7,:)') + u(tri(2))*phi2(z(7,:)') + u(tri(3))*phi3(z(7,:)'))^2 ...
            );
        
        F_i(ind) = tri(3);
        F_j(ind) = tri(1);
        F_val(ind) = e;
        ind = ind + 1;
        
        F_i(ind) = tri(1);
        F_j(ind) = tri(3);
        F_val(ind) = e;
        ind = ind + 1;
        
        
        %grad phi1 grad phi1
        e = detBT*0.5*( ...
            w(1)*phi1(z(1,:)')*phi1(z(1,:)')*abs(u(tri(1))*phi1(z(1,:)') + u(tri(2))*phi2(z(1,:)') + u(tri(3))*phi3(z(1,:)'))^2 ...
            + w(2)*phi1(z(2,:)')*phi1(z(2,:)')*abs(u(tri(1))*phi1(z(2,:)') + u(tri(2))*phi2(z(2,:)') + u(tri(3))*phi3(z(2,:)'))^2 ...
            + w(2)*phi1(z(3,:)')*phi1(z(3,:)')*abs(u(tri(1))*phi1(z(3,:)') + u(tri(2))*phi2(z(3,:)') + u(tri(3))*phi3(z(3,:)'))^2 ...
            + w(2)*phi1(z(4,:)')*phi1(z(4,:)')*abs(u(tri(1))*phi1(z(4,:)') + u(tri(2))*phi2(z(4,:)') + u(tri(3))*phi3(z(4,:)'))^2 ...
            + w(3)*phi1(z(5,:)')*phi1(z(5,:)')*abs(u(tri(1))*phi1(z(5,:)') + u(tri(2))*phi2(z(5,:)') + u(tri(3))*phi3(z(5,:)'))^2 ...
            + w(3)*phi1(z(6,:)')*phi1(z(6,:)')*abs(u(tri(1))*phi1(z(6,:)') + u(tri(2))*phi2(z(6,:)') + u(tri(3))*phi3(z(6,:)'))^2 ...
            + w(3)*phi1(z(7,:)')*phi1(z(7,:)')*abs(u(tri(1))*phi1(z(7,:)') + u(tri(2))*phi2(z(7,:)') + u(tri(3))*phi3(z(7,:)'))^2 ...
            );
        
        F_i(ind) = tri(1);
        F_j(ind) = tri(1);
        F_val(ind) = e;
        ind = ind + 1;
        
        
        %grad phi2 grad phi2
        e = detBT*0.5*( ...
            w(1)*phi2(z(1,:)')*phi2(z(1,:)')*abs(u(tri(1))*phi1(z(1,:)') + u(tri(2))*phi2(z(1,:)') + u(tri(3))*phi3(z(1,:)'))^2 ...
            + w(2)*phi2(z(2,:)')*phi2(z(2,:)')*abs(u(tri(1))*phi1(z(2,:)') + u(tri(2))*phi2(z(2,:)') + u(tri(3))*phi3(z(2,:)'))^2 ...
            + w(2)*phi2(z(3,:)')*phi2(z(3,:)')*abs(u(tri(1))*phi1(z(3,:)') + u(tri(2))*phi2(z(3,:)') + u(tri(3))*phi3(z(3,:)'))^2 ...
            + w(2)*phi2(z(4,:)')*phi2(z(4,:)')*abs(u(tri(1))*phi1(z(4,:)') + u(tri(2))*phi2(z(4,:)') + u(tri(3))*phi3(z(4,:)'))^2 ...
            + w(3)*phi2(z(5,:)')*phi2(z(5,:)')*abs(u(tri(1))*phi1(z(5,:)') + u(tri(2))*phi2(z(5,:)') + u(tri(3))*phi3(z(5,:)'))^2 ...
            + w(3)*phi2(z(6,:)')*phi2(z(6,:)')*abs(u(tri(1))*phi1(z(6,:)') + u(tri(2))*phi2(z(6,:)') + u(tri(3))*phi3(z(6,:)'))^2 ...
            + w(3)*phi2(z(7,:)')*phi2(z(7,:)')*abs(u(tri(1))*phi1(z(7,:)') + u(tri(2))*phi2(z(7,:)') + u(tri(3))*phi3(z(7,:)'))^2 ...
            );
        
        F_i(ind) = tri(2);
        F_j(ind) = tri(2);
        F_val(ind) = e;
        ind = ind + 1;
        
        
        %grad phi3 grad phi3
        e = detBT*0.5*( ...
            w(1)*phi3(z(1,:)')*phi3(z(1,:)')*abs(u(tri(1))*phi1(z(1,:)') + u(tri(2))*phi2(z(1,:)') + u(tri(3))*phi3(z(1,:)'))^2 ...
            + w(2)*phi3(z(2,:)')*phi3(z(2,:)')*abs(u(tri(1))*phi1(z(2,:)') + u(tri(2))*phi2(z(2,:)') + u(tri(3))*phi3(z(2,:)'))^2 ...
            + w(2)*phi3(z(3,:)')*phi3(z(3,:)')*abs(u(tri(1))*phi1(z(3,:)') + u(tri(2))*phi2(z(3,:)') + u(tri(3))*phi3(z(3,:)'))^2 ...
            + w(2)*phi3(z(4,:)')*phi3(z(4,:)')*abs(u(tri(1))*phi1(z(4,:)') + u(tri(2))*phi2(z(4,:)') + u(tri(3))*phi3(z(4,:)'))^2 ...
            + w(3)*phi3(z(5,:)')*phi3(z(5,:)')*abs(u(tri(1))*phi1(z(5,:)') + u(tri(2))*phi2(z(5,:)') + u(tri(3))*phi3(z(5,:)'))^2 ...
            + w(3)*phi3(z(6,:)')*phi3(z(6,:)')*abs(u(tri(1))*phi1(z(6,:)') + u(tri(2))*phi2(z(6,:)') + u(tri(3))*phi3(z(6,:)'))^2 ...
            + w(3)*phi3(z(7,:)')*phi3(z(7,:)')*abs(u(tri(1))*phi1(z(7,:)') + u(tri(2))*phi2(z(7,:)') + u(tri(3))*phi3(z(7,:)'))^2 ...
            );
        
        F_i(ind) = tri(3);
        F_j(ind) = tri(3);
        F_val(ind) = e;
        ind = ind + 1;
        
    end
    
    del_index = find(F_i==0,1);
    if isempty(del_index)
        del_index = length(F_i)+1;
    end
    
    F_worker = sparse(F_i(1:del_index-1),F_j(1:del_index-1),F_val(1:del_index-1),N_h,N_h); 
end

% collect results from workers
F = sparse(N_h,N_h);
for j = 1:length(my_index_start)
    F = F + F_worker{j};
end

end

function y = phi1(x)
y = -x(1) - x(2) + 1;
end

function y = phi2(x)
y = x(1);
end

function y = phi3(x)
y = x(2);
end


