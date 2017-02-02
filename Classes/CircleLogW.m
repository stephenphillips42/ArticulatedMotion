classdef CircleLogW < ProjectionParticleFilterSim
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
        function sim = CircleLogW(varargin)
            % Parameters for parent classes
            varargin = [
                {'dim_space',2,'dim_meas',1}, varargin
                ];
            sim@ProjectionParticleFilterSim(varargin{:});
            definedOrDefault = @(name,default) ...
                                definedOrDefault_long(name,default,varargin);
            % Motion model parameters
            sim.v = (1*(pi/180));
            sim.sigma_h = definedOrDefault('sigma_h',sqrt(0.1*(pi/180)));
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
        % Measurement Likelihood
        function l = h_likelihood(sim,x,z)
            d = bsxfun(@minus,sim.h(x,0),z);
            l = sum(sim.lognormpdf(d,0,sim.sigma_h),1);
        end

        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function v = estimate(~,samples,w)
            v = normc(sum(bsxfun(@times,samples,w),2));
        end
        
        function [samples,w] = create_samples(sim)
            samples = normc(randn(sim.dim_space,sim.n_samples));
            w = sim.uniform();
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
            scatter(samples(1,:),samples(2,:),exp(w)*(12*2^11),'b.')
            hold off
            axis equal
            axis([-1.2 1.2 -1.2 1.2])
        end


        %%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%
        function w = normalize(~,w)
            % Working in log space
            w = w - log(sum(exp(w)));
        end
        function w = uniform(sim)
            w = ones(1,sim.n_samples)*(-log(sim.n_samples));
        end
        function Neff = compute_Neff(~,w)
            Neff = 1/sum(exp(2*w));
        end

        function [samples,w] = resample(sim,samples,w)
            % Basically our standard algorithm but in log 
            hist_edges = min([-inf sim.log_cum_sum(w)],0);
            % get the upper edge exact
            hist_edges(end) = 0;
            U1 = rand/sim.n_samples;
            log_intervals = log(U1:(1/sim.n_samples):1);
            % this works like the inverse of the empirical distribution and returns
            % the interval where the sample is to be found
            [~, idx] = histc(log_intervals, hist_edges);
            samples = samples(:,idx);
            w = sim.uniform();
        end

        function y = CircleExpMap(~,x,v,t)
            v = t*v;
            nrm_tv = sqrt(sum(v.^2,1));%norm(tv, 'fro');
            y = bsxfun(@times,x,cos(nrm_tv)) + ...
                bsxfun(@times,v,sinc(nrm_tv));
            ynrms = sqrt(sum(y.^2,1));
            y = bsxfun(@rdivide, y, ynrms);
        end

        function p = lognormpdf(~,x,mu,sigma)
            p = -log(sqrt(2*pi)*sigma)+(-(x - mu).^2/(2*sigma.^2));
        end
        
        function z = log_sum(~,x,y)
            m1 = max(x,y);
            m2 = min(x,y);
            if m1 > m2 + 200
                z = m1;
            else
                z = m1 + log1p(exp(m2-m1));
            end
        end
        
        function W = log_cum_sum(sim,w)
            W = zeros(size(w));
            W(1) = w(1);
            for i = 2:length(w)
                W(i) = sim.log_sum(W(i-1),w(i));
            end
        end

    end
    
end

