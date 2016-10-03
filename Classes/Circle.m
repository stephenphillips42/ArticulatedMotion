classdef Circle < ProjectionParticleFilterSim
%ParentParticleFilterSim Simulation of a particle filter on the
%a generic manifold

    properties
        % Filter Properties
        v        % Average velocity of motion
        sigma_f  % Variance of the noise for movement
        % Helping parameters
        J = [0, -1; 1, 0];
        circ_viz = [ cos(linspace(0,2*pi,200));
                     sin(linspace(0,2*pi,200)) ];
    end
    
    methods
        function sim = Circle(x0,T)
            % Parameters for parent classes
            nparticles = 300; % Number of particles
            dim = 2; % Dimension of the space - R^2
            sigma_h = sqrt(0.1*(pi/180));
            sim@ProjectionParticleFilterSim(x0,T,nparticles,dim,sigma_h);
            % Motion model parameters
            sim.v = (1*(pi/180));
            sim.sigma_f = (1*(pi/180)); 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Motion Model
        function [v, noise] = f(sim,x,noise)
            if nargin < 3
                noise = randn(size(x))*sim.sigma_f;
                noise = noise - bsxfun(@times,x,dot(x,noise));
            end
            vel = sim.v*sim.J*x + noise;
            v = sim.CircleExpMap(x,vel,1);
        end

        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function v = estimate(~,samples,w)
            v = normc(sum(bsxfun(@times,samples,w),2));
        end
        
        function [samples,w] = create_samples(sim)
            samples = normc(randn(sim.dim_space,sim.n_samples));
            w = sim.l1normalize(ones(1,length(samples)));
        end
        
        function [e,e_full] = compute_error(~,x_true,x_est)
            e = acos(dot(normc(x_true),normc(x_est)));
            e_full = e;
        end

        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,w,est)
            % Plot world
            l = [0.95 1.05];
            figure(1);
            plot(sim.circ_viz(1,:),sim.circ_viz(2,:),'g-')
            hold on
            % Plot ground truth position
            plot(x_gt(1)*l,x_gt(2)*l,'k-')
            % Plot measurement
            start_point = -sim.R_h.'*sim.T_h;
            end_point = sim.R_h.'*([meas;1]*10 - sim.T_h);
            plot([start_point(1),end_point(1)], ...
                 [start_point(2),end_point(2)],'r-')
             % Get and plot estimate
            plot(est(1)*l,est(2)*l,'c-')
            scatter(samples(1,:),samples(2,:),w*(12*2^11),'b.')
            hold off
            axis equal
            axis([-1.2 1.2 -1.2 1.2])
        end


        %%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = CircleExpMap(~,x,v,t)
            v = t*v;
            nrm_tv = sqrt(sum(v.^2,1));%norm(tv, 'fro');
            y = bsxfun(@times,x,cos(nrm_tv)) + ...
                bsxfun(@times,v,sinc(nrm_tv));
            ynrms = sqrt(sum(y.^2,1));
            y = bsxfun(@rdivide, y, ynrms);
        end
        
    end
    
end

