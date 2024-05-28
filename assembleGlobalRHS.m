function F = assembleGlobalRHS(f,T,Nd)

F = zeros(size(Nd,1),1);

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
    %grad phi1 grad phi1
    e = ((1/6)*detBT)*f(z1);
    index1 = tri(1);
    F(index1) = F(index1) + e;
    
    %grad phi2 grad phi2
    e = ((1/6)*detBT)*f(z2);
    
    index2 = tri(2);
    F(index2) = F(index2) + e;
    
    %grad phi3 grad phi3
    e = ((1/6)*detBT)*f(z3);
    
    index3 = tri(3);
    F(index3) = F(index3) + e;
    
end

end