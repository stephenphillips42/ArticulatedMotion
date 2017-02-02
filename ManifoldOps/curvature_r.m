function [ r ] = curvature_r(u,v,w)
% CURVATURE_R Riemannian curvature for sphere

r = bsxfun(@times,dot(w,u),v)-bsxfun(@times,dot(w,v),u);

end