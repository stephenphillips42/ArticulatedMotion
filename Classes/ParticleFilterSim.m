classdef ParticleFilterSim < handle
%Simulation of a particle filter on the a generic manifold

    properties
        %%% Filter Properties
        
        dim_space % Dimenstion of the embedding space of the manifold
        dim_meas  % Dimension of the measurement model
        n_samples % Number of particles in the filter
        x0        % Initial position of the true particle
        
        %%% Simulation Properties
        
        T         % Number of time steps simulated
        dt        % Time between steps (for visualization)
        
        %%% Filter variables
        
        simrun    % Stores information about current run
    end
    
    methods
        function sim = ParticleFilterSim(x0,T,nsamples,dim_space,dim_meas)
            % Builds the basic necessities for running the simulation

            % Filter properties
            sim.n_samples = nsamples;
            sim.x0 = x0;
            sim.dim_space = dim_space;
            sim.dim_meas = dim_meas;
            
            % Simulation properties
            sim.T = T;
            sim.dt = 0.1;

            % Get template simrun
            sim.simrun.x_gt = 0;    % Ground truth position
            sim.simrun.meas = 0;    % Measurements of position
            sim.simrun.noise_f = 0; % Movement noise
            sim.simrun.noise_h = 0; % Measurment noise
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Run Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%
        function [results] = simulate(sim,plotting,simrun)
            % Runs the simulation of the particle filter.
            % Based on parameters set in the constructor. 
            % Input:
            %   plotting - (optional) level of plotting desired, from 0-3
            %   simrun   - (optional) precomputed simulation run
            % Output:
            %   results - numerical results of the simulation
            if nargin < 2
                plotting = 1;
            end
            if nargin < 3
                sim.simrun = sim.create_simrun();
            else
                sim.simrun = simrun;
            end

            [samples,w] = sim.create_samples();
            Neff = zeros(1,sim.T);
            est = zeros(size(sim.simrun.x_gt,1),sim.T);
            % Visualization
            for i = 1:sim.T
                % Get metrics
                w = sim.l1normalize(sim.h_likelihood(samples,sim.simrun.meas(:,i)));
                Neff(i) = 1/sum(w.^2);
                % Plot reweighted particles
                if plotting > 1
                    sim.plot_simulation(...
                        sim.simrun.x_gt(:,i),...
                        sim.simrun.meas(:,i),...
                        samples,w,...
                        est(:,k));
                end
                % Resampling
                if sim.should_resample(Neff(i),i)
                    [samples,w] = sim.resample(samples,w);
                end
                % Estimate
                est(:,i) = sim.estimate(samples,w);
                % Plot final resampled (or not) particles
                if plotting > 0
                    sim.plot_simulation(...
                        sim.simrun.x_gt(:,i),...
                        sim.simrun.meas(:,i),...
                        samples,w,...
                        est(:,i));
                    % Pause for viewing
                    pause(sim.dt)
                end
                % Propogate
                samples = sim.f(samples);
                k = max(1,i-1);
                % Plot propgated particles
                if plotting > 2
                    sim.plot_simulation(...
                        sim.simrun.x_gt(:,i),...
                        sim.simrun.meas(:,i),...
                        samples,w,...
                        est(:,k));
                end
            end
            % Save results
            results.Neff = Neff;
            results.est = est;
            results.simrun = sim.simrun;
        end
        
        function [simrun] = create_simrun(sim)
            % Precompute motion before simulation
            simrun.x_gt = zeros(size(sim.x0,1),sim.T);
            simrun.x_gt(:,1) = sim.x0;
            % Get appropriate sized storage for noise and measurements
            [~,ntst] = sim.f(simrun.x_gt(:,1));
            simrun.noise_f = zeros(size(ntst,1),sim.T);
            [mtst,ntst] = sim.h(simrun.x_gt(:,1));
            simrun.meas = zeros(size(mtst,1),sim.T);
            simrun.noise_h = zeros(size(ntst,1),sim.T);
            % Compute the motion, measurements, measurements, and noise
            [simrun.meas(:,1), simrun.noise_h(:,1)] = ...
                sim.h(simrun.x_gt(:,1));
            for i = 2:sim.T
                [simrun.x_gt(:,i), simrun.noise_f(:,i)] = ...
                    sim.f(simrun.x_gt(:,i-1));
                [simrun.meas(:,i), simrun.noise_h(:,i)] = ...
                    sim.h(simrun.x_gt(:,i));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Motion Model
        function [v, noise] = f(sim,x,noise)
            if nargin < 3
                noise = zeros(sim.dim_space,size(x,2));
            end
            v = x;
        end
        % Measurement Model
        function [meas, noise] = h(sim,x,noise)
            if nargin < 3
                noise = zeros(sim.dim_meas,size(x,2))*sim.sigma_h;
            end
            meas = x(1:sim.dim_meas,:);
        end
        % Measurement Likelihood
        function l = h_likelihood(sim,x,z)
            l = 1;
        end
        
        % Error function
        function e = compute_error(~,x1,x2)
            e = norm(x1-x2);
        end
        
        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function [samples,w] = create_samples(sim)
            samples = (randn(sim.dim_space,sim.n_samples));
            w = sim.l1normalize(ones(1,length(samples)));
        end
        
        function [b] = should_resample(sim,Neff_cur,i)
            b =  Neff_cur < sim.n_samples*(1/2);
        end
        
        function [samples,w] = resample(sim,samples,w)
            % Taken from Diego Andrés Alvarez Marín's Particle Filter code
            % this is performing latin hypercube sampling on w
            % protect against accumulated round-off
            edges = min([0 cumsum(w)],1);
            % get the upper edge exact
            edges(end) = 1;
            U1 = rand/sim.n_samples;
            % this works like the inverse of the empirical distribution and returns
            % the interval where the sample is to be found
            [~, idx] = histc(U1:(1/sim.n_samples):1, edges);
            samples = samples(:,idx);
            w = sim.l1normalize(ones(1,length(samples)));
        end
        
        function v = estimate(~,samples,w)
            v = sum(bsxfun(@times,samples,w),2);
        end
        
        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,w,est)
            % To implement
        end

        %%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%
        function v = l1normalize(~,x)
            v = x/norm(x,1);
        end

        
    end
    
end

