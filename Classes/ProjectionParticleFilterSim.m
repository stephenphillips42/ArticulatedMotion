classdef ProjectionParticleFilterSim < ParticleFilterSim
%ParentParticleFilterSim Simulation of a particle filter on the
%a generic manifold

    properties
        % Filter Properties
        sigma_h   % Noise level of the projection
        proj_dims % Dimensions to be projected - default: 1:dim_space-1
        div_dim   % Dimension to divide the others by - default: dim_space
        R_h       % Rotation/Linear part before projection
        T_h       % Offset/Affine part before projection
    end
    
    methods
        function sim = ProjectionParticleFilterSim(x0,T,nsamples,dim_space,sigma_h,proj_dims,div_dim,R_h,T_h)
            if nargin < 6
                proj_dims = 1:(dim_space-1);
            end
            if nargin < 7
                div_dim = dim_space;
            end
            if nargin < 8
                R_h = eye(dim_space);
            end
            if nargin < 9
                T_h = [zeros(dim_space-1,1);dim_space];
            end
            sim@ParticleFilterSim(x0,T,nsamples,dim_space,length(proj_dims));
            sim.proj_dims = proj_dims;
            sim.div_dim = div_dim;
            sim.sigma_h = sigma_h;
            sim.R_h = R_h;
            sim.T_h = T_h;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Measurement Model - standard projection model
        function [meas, noise] = h(sim,x,noise)
            if nargin < 3
                noise = randn(sim.dim_meas,size(x,2))*sim.sigma_h;
            end
            meas = sim.Proj(bsxfun(@plus,sim.R_h*x,sim.T_h)) + noise;
        end
        % Measurement Likelihood
        function l = h_likelihood(sim,x,z)
            d = bsxfun(@minus,sim.h(x,0),z);
            l = prod(normpdf(d,0,sim.sigma_h),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%
        function v = Proj(sim,x)
            v = bsxfun(@rdivide,x(sim.proj_dims,:),x(sim.div_dim,:));
        end
        
    end
    
end

