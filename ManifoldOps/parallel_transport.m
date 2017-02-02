function zetat = parallel_transport(x,xt,zeta0)
% PARALLEL_TRANSPORT Parallel transport of the sphere
    x = mynormc(x);
    xt = mynormc(xt);

    % Make sure zeta0 is in the tangent space of x
    zeta0 = zeta0 - bsxfun(@times,x,dot(x,zeta0));

    % Initialize zetat
    zetat = zeta0;
    % Transport zeta0 from the tangent space of x to the tangent space of xt
    xi = log_map_sphere(x,xt);
    nrm_xi = sqrt(sum(xi.^2,1));
    valid_inds = (nrm_xi > 10^-9); % TODO Make threshold not magic number
    if any(valid_inds)
        nrm_xi = nrm_xi(valid_inds);
        u = mynormc(xi(:,valid_inds));
        u_dot_zeta = dot(u,zeta0(:,valid_inds));
        zetat(:,valid_inds) = ...
                -bsxfun(@times,x(:,valid_inds),sin(nrm_xi).*u_dot_zeta) + ...
                 bsxfun(@times,u,cos(nrm_xi).*u_dot_zeta) + ...
                 zeta0(:,valid_inds) - bsxfun(@times,u,u_dot_zeta);   
    end
end
