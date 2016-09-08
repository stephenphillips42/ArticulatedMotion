classdef CircleRiemannianParticleFilterSim
%ManOptParticleFilterSim Simulation of a particle filter on the circle
%using generic manifold methods
%   This will provide the template for later particle filters - right now
%   it has the following motion model:
%     Space: x in { z in R^2 | norm(z) == 1 }
%     Motion model: f(x) = exp(x,v + N(0,sigma_f))
%     Measurement model: h(x) = Proj(R_h*x + T_h) + N(0,sigma_h)
%   where N denotes the normal distribution, Proj(x) = x(1)/x(2). 

    properties
        % Filter Properties
        v         % tangent speed of the particle on the sphere
        sigma_f   % inverse variance of the noise for movement
        sigma_h   % standard devation of noise for the measurement
        R_h       % Orientation of the camera
        T_h       % Translation of camera
        n_samples % Number of particles in the filter
        x0        % Initial position of the true particle
        % Simulation Properties
        T         % Number of time steps simulated
        dt        % Time between steps (for visualization)
        % Filter variables
        x_gt      % Ground truth of the particle
        meas      % Measurements given to the particle
        noise_f   % Movement noise
        noise_h   % Measurement noise
        % Visualization variables
        circ_viz 
    end
    
    methods
        function sim = CircleRiemannianParticleFilterSim(theta0,T)
            if nargin < 1
                theta0 = pi/2;
            end
            if nargin < 2
                T = 100;
            end
            % Filter properties
            sim.v = (1*(pi/180));
            sim.sigma_f = (1*(pi/180)); 
            sim.sigma_h = sqrt(0.1*(pi/180));
            sim.R_h = eye(2);%sim.R(pi/3);
            sim.T_h = [0;3];
            sim.n_samples = 300;
            sim.x0 = [ cos(theta0); sin(theta0) ];
            
            % Simulation properties
            sim.T = T;
            sim.dt = 0.01;
            
            % Filter variables (all precomputed before simulation)
            sim.x_gt = zeros(2,sim.T);
            sim.meas = zeros(1,sim.T);
            sim.noise_f = normrnd(0, sim.sigma_f, 1, sim.T-1);
            sim.noise_h = normrnd(0, sim.sigma_h, 1, sim.T);
            sim.x_gt(:,1) = sim.x0;
            sim.meas(1) = sim.h(sim.x_gt(:,1),sim.noise_h(1));
            for i = 2:sim.T
                sim.x_gt(:,i) = sim.f(sim.x_gt(:,i-1),sim.noise_f(i-1));
                sim.meas(i) = sim.h(sim.x_gt(:,i),sim.noise_h(i));
            end
            
            % Visualization Variables
            th_viz = linspace(0,2*pi,200);
            sim.circ_viz = [ cos(th_viz); sin(th_viz) ];
        end
        
        
        function [results] = simulate(sim,plotting)
            if nargin < 2
                plotting = true;
            end
            [samples,w] = sim.create_samples();
            Neff = zeros(1,sim.T);
            est = zeros(2,sim.T);
            % Visualization
            for i = 1:sim.T
                % Get metrics
                est(:,i) = sim.estimate(samples,w);
                Neff(i) = 1/sum(w.^2);
                % Plotting code
                if plotting
                    sim.plot_simulation(sim.x_gt(:,i),sim.meas(i),samples,w,est,Neff,i)
                end
                pause(sim.dt)
                % Propogate
                [samples,w] = sim.propogate_samples(samples,w,sim.meas(i));
                % Resampling
                if Neff(i) < sim.n_samples/2
                    [samples,w] = resample(sim,samples,w);
                end
            end
            results.Neff = Neff;
            results.est = est;
            results.simrun.x_gt = sim.x_gt;
        end
        
        %%% Model of the system
        % Motion Model
        function v = f(sim,x,noise)
            if nargin < 3
                noise = (eye(2)-x*x.')*randn(2,1)*sim.sigma_f;
            end
            vel = sim.v*[ 0 -1; 1 0 ]*x + noise;
            v = sim.circle_exponential(x,vel,1);
        end
        % Measurement Model
        function v = h(sim,x,noise)
            if nargin < 3
                noise = normrnd(0,sim.sigma_h);
            end
            v = sim.Proj(sim.R_h*x+sim.T_h) + noise;
        end
        % Measurement Likelihood
        function v = h_likelihood(sim,x,z)
            v = normpdf(sim.h(x,0)-z,0,sim.sigma_h);
        end
        
        %%% Particle filter functions
        function [samples,w] = create_samples(sim)
            sample_thetas = vmrand(atan2(sim.x0(2),sim.x0(1)), 0.01, 1, sim.n_samples);
            samples = [ cos(sample_thetas); sin(sample_thetas) ];
            w = sim.normalize(ones(1,length(samples)));
        end
        
        function [samples,w] = propogate_samples(sim,samples,w,meas)
            for j = 1:sim.n_samples
                samples(:,j) = sim.f(samples(:,j));
                w(j) = sim.h_likelihood(samples(:,j),meas);
            end
            w = sim.normalize(w);
        end
        
        function [samples,w] = resample(sim,samples,w)
            % Taken from Diego Andrés Alvarez Marín's Particle Filter code
            % this is performing latin hypercube sampling on w
            edges = min([0 cumsum(w)],1); % protect against accumulated round-off
            edges(end) = 1;                 % get the upper edge exact
            U1 = rand/sim.n_samples;
            % this works like the inverse of the empirical distribution and returns
            % the interval where the sample is to be found
            [~, idx] = histc(U1:(1/sim.n_samples):1, edges);
            samples = samples(:,idx);
            w = sim.normalize(ones(1,length(samples)));
        end

        %%% Visualization functions
        function plot_simulation(sim,x_gt,meas,samples,w,est,Neff,i)
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
            plot(est(1,i)*l,est(2,i)*l,'c-')
            scatter(samples(1,:),samples(2,:),w*(12*2^11),'b.')
            hold off
            axis equal
            axis([-1.2 1.2 -1.2 1.2])
            %legend({'Circle','True','Measurement','Estimate','Samples'},...
            %        'Location','northoutside')
            % figure(2)
            % plot(Neff)
        end
        
        %%% Helper functions
        function v = Proj(~,x)
            v = x(1)/x(2);
        end
        
        function outR = R(~,theta)
            outR = [ cos(theta), -sin(theta); sin(theta), cos(theta) ];
        end
        
        function v = normalize(~,x)
            v = x/norm(x,1);
        end
        
        function v = estimate(~,samples,w)
            v = normc(sum(samples.*(ones(2,1)*w),2));
        end
        
        function y = circle_exponential(~, x, v, t)
            if nargin < 4
                t = 1;
            end

            tv = t*v;

            nrm_tv = norm(tv, 'fro');

            if nrm_tv > 4.5e-8
                y = normc(x*cos(nrm_tv) + tv*(sin(nrm_tv)/nrm_tv));
            else
                % If the step is too small to accurately evaluate sin(x)/x,
                % then sin(x)/x is almost indistinguishable from 1.
                y = x + tv;
                y = y / norm(y, 'fro');
            end
        end
        
    end
    
end

