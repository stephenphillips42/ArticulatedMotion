classdef ProjectionParticleFilterSim < ParticleFilterSim
%ProjectionParticleFilterSim Simulation of a particle filter on the
%a generic manifold, specifically with a projection model of measurement

    properties
        % Filter Properties
        sigma_h   % Noise level of the projection
        proj_dims % Dimensions to be projected - default: 1:dim_space-1
        div_dim   % Dimension to divide the others by - default: dim_space
        R_h       % Rotation/Linear part before projection
        T_h       % Offset/Affine part before projection
    end
    
    methods
        function sim = ProjectionParticleFilterSim(varargin)
            definedOrDefault = @(name,default) ...
                                definedOrDefault_long(name,default,varargin);
            sim@ParticleFilterSim(varargin{:});
            sim.proj_dims = definedOrDefault('proj_dims',1:(sim.dim_space-1));
            sim.dim_meas = length(sim.proj_dims);
            sim.div_dim = definedOrDefault('div_dim', sim.dim_space);
            sim.sigma_h = definedOrDefault('sigma_h', 0.1);
            sim.R_h = definedOrDefault('R_h', eye(sim.dim_space));
            sim.T_h = definedOrDefault('T_h',[ zeros(length(sim.proj_dims),1);
                                               length(sim.dim_space)+1       ]);
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
            l = sum(lognormpdf(d,0,sim.sigma_h),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%
        function v = Proj(sim,x)
            v = bsxfun(@rdivide,x(sim.proj_dims,:),x(sim.div_dim,:));
        end
        
    end
    
end

