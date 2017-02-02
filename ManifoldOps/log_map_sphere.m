function mlog = log_map_sphere(x,y)
% LOG_MAP_SPHERE Log map of the sphere manifold
%   Takes in x, size N by K, assuming each column is a point on the N-sphere
%   and y, also size N by K. 
%   Does not do checking to ensure input args are on the manifold

% Get geodesic distance between x and y
trxy = max(-1,min(dot(x,y),1));
geodist = acos(trxy) ;

% TODO Is there a better way of writing the log map?
mlog = bsxfun(@times, geodist, bsxfun(@rdivide,...
                y-bsxfun(@times,x,trxy),...
                eps+sqrt(1-trxy.^2)));

end
