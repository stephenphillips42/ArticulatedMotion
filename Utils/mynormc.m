function x = mynormc(x)
% MYNORMC Simplified version of normc from MATLAB
    x = bsxfun(@rdivide,x,eps+sqrt(sum(x.^2,1)));
end
