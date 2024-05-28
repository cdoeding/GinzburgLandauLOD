function U = evaluateInitialValue(u0,T,Nd,M_LOD,Phi)
%EVALUATEINITIALVALUE Summary of this function goes here
%   Detailed explanation goes here

r = assembleGlobalRHS(u0,T,Nd);
U = M_LOD'\(Phi*r);

end

