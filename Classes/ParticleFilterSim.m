classdef ParticleFilterSim < handle
%Simulation of a particle filter on the a generic manifold

    properties
        %%% Filter Properties
        
        dim_space % Dimenstion of the embedding space of the manifold
        dim_meas  % Dimension of the measurement model
        n_samples % Number of particles in the filter
        x0        % Initial position of the true particle
        sigma_s   % Resample variance (could be multidimensional)
        
        %%% Simulation Properties
        
        T         % Number of time steps simulated
        dt        % Time between steps (for visualization)
        
        %%% Filter variables
        
        simrun    % Stores information about current run
        pows      % Powers for annealing
    end
    
    methods
        function sim = ParticleFilterSim(varargin)
            % Builds the basic necessities for running the simulation
            definedOrDefault = @(name,default) ...
                         definedOrDefault_long(name,default,varargin);
            % Filter properties
            sim.n_samples = definedOrDefault('nsamples',300);
            sim.x0 = definedOrDefault('x0',[]);
            sim.dim_space = definedOrDefault('dim_space',3);
            sim.dim_meas = definedOrDefault('dim_meas',sim.dim_space-1);
            sim.sigma_s = definedOrDefault('sigma_s',0.0);

            % Simulation properties
            sim.T = definedOrDefault('nsteps',10);
            sim.dt = definedOrDefault('dt',0.1);

            % Get template simrun
            sim.simrun.x_gt = 0;    % Ground truth position
            sim.simrun.meas = 0;    % Measurements of position
            sim.simrun.noise_f = 0; % Movement noise
            sim.simrun.noise_h = 0; % Measurment noise
            
            % For annealing
            sim.pows = definedOrDefault('pows',[0.1,0.3,0.5,0.7,1.0]);
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
                fprintf('Step %d\n',i)
                k = max(1,i-1);
                % Get metrics
                w = sim.normalize(sim.h_likelihood(samples,sim.simrun.meas(:,i)));
                Neff(i) = sim.compute_Neff(w);
                % Plot reweighted particles
                if plotting > 1
                    sim.plot_simulation(...
                        sim.simrun.x_gt(:,i),...
                        sim.simrun.meas(:,i),...
                        samples,w,...
                        est(:,k));
                    pause(sim.dt);
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
                % Plot propgated particles
                if plotting > 2 && i < sim.T
                    sim.plot_simulation(...
                        sim.simrun.x_gt(:,i),...
                        sim.simrun.meas(:,i),...
                        samples,w,...
                        est(:,k));
                    pause(sim.dt);
                end
            end
            % Save results
            results.Neff = Neff;
            results.est = est;
            results.simrun = sim.simrun;
        end

        function [results] = simulate_annealing(sim,plotting,simrun)
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
                fprintf('Step %d\n',i)
                l = max(1,i-1);
                % Get metrics
                w = sim.normalize(...
                        sim.h_likelihood(samples,sim.simrun.meas(:,i)));
                Neff(i) = sim.compute_Neff(w);
                % Resampling
                if sim.should_resample(Neff(i),i)
                    if plotting > 1
                        sim.plot_simulation(...
                            sim.simrun.x_gt(:,i),...
                            sim.simrun.meas(:,i),...
                            samples,w,...
                            est(:,i));
                        % Pause for viewing
                        pause(sim.dt)
                    end
                    w_prev = w;
                    for k = 1:length(sim.pows)
                        if plotting > 2
                            sim.plot_simulation(...
                                sim.simrun.x_gt(:,i),...
                                sim.simrun.meas(:,i),...
                                samples,w,...
                                est(:,i));
                            % Pause for viewing
                            pause(sim.dt)
                        end
                        w_cur = sim.weight_power(w_prev,sim.pows(k));
                        samples = sim.add_noise(sim.resample_weighted(samples,w_cur));
                        w_prev = sim.normalize(sim.h_likelihood(samples,sim.simrun.meas(:,i)));
                    end
                    w = sim.uniform();
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
                % Plot propgated particles
                if plotting > 2 && i < sim.T
                    sim.plot_simulation(...
                        sim.simrun.x_gt(:,i),...
                        sim.simrun.meas(:,i),...
                        samples,w,...
                        est(:,l));
                    pause(sim.dt);
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
        
        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function [samples,w] = create_samples(sim)
            samples = (randn(sim.dim_space,sim.n_samples));
            w = sim.uniform();
        end
        
        function [b] = should_resample(sim,Neff_cur,i)
            b =  Neff_cur < sim.n_samples*(1/2);
        end
        
        function [samples,w] = resample(sim,samples,w)
            samples = sim.add_noise(sim.resample_weighted(samples,w));
            w = sim.uniform();
        end
        
        function [ samples ] = resample_weighted(sim,samples,w)
            % Taken from Diego Andrés Alvarez Marín's Particle Filter code
            % Basically our standard algorithm but in log 

            % this is performing latin hypercube sampling on w
            % protect against accumulated round-off
            hist_edges = min([-inf sim.log_cum_sum(w)],0);
            % get the upper edge exact
            hist_edges(end) = 0;
            U1 = rand/sim.n_samples;
            log_intervals = log(U1:(1/sim.n_samples):1);
            % this works like the inverse of the empirical distribution and returns
            % the interval where the sample is to be found
            [~, idx] = histc(log_intervals, hist_edges);
            samples = samples(:,idx);
        end

        function samples = add_noise(sim,samples)
            samples = samples + sim.sigma_s*randn(size(samples));
        end

        function v = estimate(~,samples,w)
            v = sum(bsxfun(@times,samples,w),2);
        end
        
        function [e, e_full] = compute_error(~,x_true,x_est)
            % e is a single number for each sample
            % e_full is more detailed, giving a vector of errors. Sometimes
            % these two are the same
            e = sqrt(sum((x_true-x_est).^2,1));
            e_full = abs(x_true-x_est);
        end
        
        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,w,est)
            % To implement
        end

        %%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% 
        % Thank you stack overflow
        function z = log_sum(~,x,y)
        % Computes log(exp(x) + exp(y)) in a numerical stable fashion
            m1 = max(x,y);
            m2 = min(x,y);
            if m1 > m2 + 200
                z = m1;
            else
                z = m1 + log1p(exp(m2-m1));
            end
        end
        function W = log_cum_sum(sim,w)
        % Cumulative sum of log weights w using numerically stable methods
            W = zeros(size(w));
            W(1) = w(1);
            for i = 2:length(w)
                W(i) = sim.log_sum(W(i-1),w(i));
            end
        end

        function w = normalize(sim,w)
            % v = x/norm(x,1);
            % Working in log space
            W = sim.log_cum_sum(w);
            w = w - W(end);
        end
        function w = uniform(sim)
            w = ones(1,sim.n_samples)*(-log(sim.n_samples));
            % w = ones(1,sim.n_samples)/sim.n_samples;
        end
        function Neff = compute_Neff(~,w)
            % Neff = 1/sum(w.^2);
            Neff = 1/sum(exp(2*w));
        end
        function w_k = weight_power(sim,w,p)
            w_k = sim.normalize(w*p);
        end
    end
end

