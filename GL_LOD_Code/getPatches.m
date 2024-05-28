function patch = getPatches(TH,ell)
%GETPATCHES Summary of this function goes here
%   Detailed explanation goes here

d = 2;
NH = size(TH.p,1);
NTH = size(TH.t,1);

Ivt = sparse(TH.t,repmat((1:NTH).',1,d+1),1,NH,NTH);
Itt = spones(Ivt.'*Ivt);
patch = speye(NTH);

for k = 1:ell
    patch = Itt*patch;
end

patch = spones(patch)+speye(NTH);
end

