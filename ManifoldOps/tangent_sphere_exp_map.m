function [y,z] = tangent_sphere_exp_map(x,v,t,L)
% TANGENT_SPHERE_EXP_MAP Numerically computes exponential map of tangent bundle
% of the sphere
%   Takes in x, size 2*N by K, assuming each column is a point on the
%   N-sphere tangent bundle, and v, also size 2*N by K, of the velocities of
%   the exponential map, and t, the expansion of v, and L, the number of
%   iterations 
%   Does not do checking to ensure input args are on the manifold

if nargin < 4
    t = 1;
end
if nargin < 5
    L = 50;
end

% Rename different parts of the exponential map
% assert(mod(size(x,1),2) == 0)
n = size(x,1)/2;
p = x(1:n,:);
u = x(n+1:2*n,:);
w = t*v(n+1:2*n,:);
v = t*v(1:n,:);
% Iterate to get approx. of the exp. map
stepsize = 1/L; % How many iterations
for i=1:L
     pprev = p;
     uprev = u;
     vprev = v;
     wprev = w;

     p = sphere_exp_map(p,v,stepsize);
     u = parallel_transport(pprev,p,uprev+stepsize*wprev);
     v = parallel_transport(pprev,p,vprev-stepsize*curvature_r(uprev,wprev,vprev));
     w = parallel_transport(pprev,p,wprev);
end
y = [p;u];
z = [v;w]; % Don't think this is necessary

end
