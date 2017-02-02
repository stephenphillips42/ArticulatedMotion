function [ y ] = sphere_exp_map(x,v,t)
%SPHERE_EXP_MAP Exponential map of the sphere maniofld
%   Takes in x, size N by K, assuming each column is a point on the N-sphere
%   and y, also size N by K. 
%   Does not do checking to ensure input args are on the manifold

v = t*v;
nrm_tv = sqrt(sum(v.^2,1));%norm(tv, 'fro');
y = bsxfun(@times,x,cos(nrm_tv)) + ...
    bsxfun(@times,v,sinc(nrm_tv));
y = mynormc(y);

end


