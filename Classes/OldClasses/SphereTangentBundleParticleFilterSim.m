classdef SphereTangentBundleParticleFilterSim
%SphereRiemannianParticleFilterSim Simulation of a particle filter on the
%sphere using generic manifold methods
%   This will provide the template for later particle filters - right now
%   it has the following motion model:
%     Space: x in { z in R^3 | norm(z) == 1 }
%     Motion model: f(x) = exp(x,N(0,sigma_f))
%     Measurement model: h(x) = Proj(R_h*x + T_h) + N(0,sigma_h)
%   where N denotes the normal distribution (isotropic), Proj(x) = x(1:2)/x(3). 

    properties
        % Filter Properties
        v         % tangent speed of the particle on the sphere
        sigma_f   % inverse variance of the noise for movement
        omega     % angular velocity of the movement
        omegahat  % skew-symmetric version of omega for conveninece
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
        X
        Y
        Z
    end
    
    methods
        function sim = SphereRiemannianParticleFilterSim(x0,T)
            if nargin < 1
                x0 = normc([1;1;1]);
            end
            if nargin < 2
                T = 100;
            end
            % Filter properties
            sim.sigma_f = sqrt(1*(pi/180)); % Motion
            sim.omega = [0;0;1];
            sim.omegahat = sim.hat(sim.omega);
            sim.sigma_h = sqrt(0.02*(pi/180)); % Measurement
            sim.R_h = eye(3);
            sim.T_h = [0;0;3];
            sim.n_samples = 10000;
            sim.x0 = x0;
            
            % Simulation properties
            sim.T = T;
            sim.dt = 0.1;
            
            % Filter variables (all precomputed before simulation)
            sim.x_gt = zeros(3,sim.T);
            sim.meas = zeros(2,sim.T);
            sim.noise_f = zeros(3, sim.T-1);
            sim.noise_h = zeros(2, sim.T);
            sim.x_gt(:,1) = normc(sim.x0);
            [sim.meas(:,1), sim.noise_h(:,1)] = sim.h(sim.x_gt(:,1));
            for i = 2:sim.T
                [sim.x_gt(:,i), sim.noise_f(:,i)] = sim.f(sim.x_gt(:,i-1));
                [sim.meas(:,i), sim.noise_h(:,i)] = sim.h(sim.x_gt(:,i));
            end
            
            % Visualization Variables
            [sim.X, sim.Y, sim.Z] = sphere(20);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Run Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%
        function [results] = simulate(sim,plotting)
            if nargin < 2
                plotting = true;
            end
            [samples,w] = sim.create_samples();
            Neff = zeros(1,sim.T);
            est = zeros(size(sim.x_gt,1),sim.T);
            % Visualization
            for i = 1:sim.T
                % Get metrics
                Neff(i) = 1/sum(w.^2);
                % Propogate
                samples = sim.f(samples);
                w = sim.l1normalize(sim.h_likelihood(samples,sim.meas(:,i)));
                % Resampling
                if sim.should_resample(Neff(i),i)
                    [samples,w] = sim.resample(samples,w);
                end
                % Estimate
                est(:,i) = sim.estimate(samples,w);
                % Plotting code
                if plotting
                    sim.plot_simulation_2D(sim.x_gt(:,i),sim.meas(:,i),samples,w,est(:,i));
                end
                % Pause for viewing
                pause(sim.dt)
            end
            results.Neff = Neff;
            results.est = est;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Motion Model
        function [v, noise] = f(sim,x,noise)
            if nargin < 3
                % noise = zeros(3,1);
                noise = randn(3,size(x,2))*sim.sigma_f;
            end
            mv = sim.dt*sim.omegahat*x + (eye(3)-x*x.')*noise;
            v = sim.ExpMap(x,mv,1);
        end
        % Measurement Model
        function [meas, noise] = h(sim,x,noise)
            if nargin < 3
                noise = randn(2,size(x,2))*sim.sigma_h;
            end
            meas = sim.Proj(bsxfun(@plus,sim.R_h*x,sim.T_h)) + noise;
        end
        % Measurement Likelihood
        function l = h_likelihood(sim,x,z)
            d = bsxfun(@minus,sim.h(x,0),z);
            l = normpdf(d(1,:),0,sim.sigma_h).*normpdf(d(2,:),0,sim.sigma_h);
        end
        
        %%% Particle filter functions
        function [samples,w] = create_samples(sim)
            samples = normc(randn(3,sim.n_samples));
            w = sim.l1normalize(ones(1,length(samples)));
        end
        
        function [b] = should_resample(sim,Neff_cur,i)
            b =  Neff_cur < sim.n_samples*(1/2); % || mod(i,40) == 0;
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

        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,w,est)
            sim.plot_simulation_2D(x_gt,meas,samples,w,est)
        end

        function plot_simulation_3D(sim,x_gt,meas,samples,w,est)
            % Plot world
            colormap summer;
            C = ones(size(sim.X));
            C(:) = sim.h_likelihood([sim.X(:) sim.Y(:) sim.Z(:)].',meas);
            surf(sim.X,sim.Y,sim.Z,C);
            axis equal
            hold on
            % Plot samples
            scatter3(samples(1,:),samples(2,:),samples(3,:),eps+w*(2^14),'b.')
            % Plot measurement
            start_point = -sim.R_h.'*sim.T_h;
            end_point = sim.R_h.'*([meas;1]*10 - sim.T_h);
            plot3([start_point(1),end_point(1)], ...
                 [start_point(2),end_point(2)], ...
                 [start_point(3),end_point(3)],'r-')
            % Plot estimate
            scatter3(est(1),est(2),est(3),500,'c.')
            % Plot ground truth
            scatter3(x_gt(1),x_gt(2),x_gt(3),500,'r.')
            axis equal
            axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])
            hold off
        end
        
        function plot_simulation_2D(sim,x_gt,meas,samples,w,est)
            % Plot likelihood
            [a,b] = meshgrid(-pi:0.01:pi,0:0.01:pi);
            z = zeros(size(a));
            z(:) = sim.h_likelihood(sim.cartesian_coords([a(:) b(:)].'),meas);
            figure(1)
            contourf(a,b,z,5)
            hold on
            % Get spherical coordinates
            samples_s = sim.spherical_coords(samples);
            % Plot samples
            scatter(samples_s(1,:),samples_s(2,:),eps+w*(2^14),'b.')
            hold on
            % Plot measurements
            %sim.plot_meas_2D(meas)
            % Plot estimate
            est_s = sim.spherical_coords(est);
            scatter(est_s(1),est_s(2),500,'c.')
            % Plot ground truth
            x_gt_s = sim.spherical_coords(x_gt);
            scatter(x_gt_s(1),x_gt_s(2),600,'r.')
            axis equal
            axis([(-pi), (pi), (0), (pi) ])
            hold off
        end

        function plot_meas_2D(sim,meas)
            s = -sim.R_h.'*sim.T_h; % Start point
            e = sim.R_h.'*([meas;1]*10 - sim.T_h); % End point e
            t1 = (-s.'*e + sqrt((s.'*e)^2 - (e.'*e)*(s.'*s-1)))/(e.'*e);
            if isreal(t1)
                t2 = (-s.'*e - sqrt((s.'*e)^2 - (e.'*e)*(s.'*s-1)))/(e.'*e);
                m1 = sim.spherical_coords(e*t1 + s); % Measurement 1
                m2 = sim.spherical_coords(e*t2 + s); % Measurement 2
                scatter([m1(1),m2(1)],[m1(2),m2(2)],600,'rx')
            else
                t_steps = 0:0.01:1;
                line = e*t_steps + s*ones(size(t_steps));
                line_norms = sqrt(sum(line.^2,1));
                meas_s = sim.spherical_coords(bsxfun(@rdivide,line,line_norms));
                % plot(meas_s(1,:),meas_s(2,:))
                surface([meas_s(1,:); meas_s(1,:)],...
                        [meas_s(2,:); meas_s(2,:)],...
                        [line_norms; line_norms],...
                        'facecol','no',...
                        'edgecol','interp',...
                        'linew',2);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%
        function v = Proj(~,x)
            v = bsxfun(@rdivide,x(1:end-1,:),x(end,:));
        end
        
        function v = l1normalize(~,x)
            v = x/norm(x,1);
        end
        
        function v = tbproj(~,X)
            v = zeros(size(X));
            v(1:3,:) = normc(X(1:3,:));
            v(4:6,:) = X(4:6) - bsxfun(@times,X(1:3),dot(X(1:3),X(4:6)));
        end
        
        function xhat = hat(~,x)
            xhat = [     0, -x(3),  x(2); 
                      x(3),     0, -x(1);
                     -x(2),  x(1),     0 ];
        end
        
        function v = estimate(~,samples,w)
            % if nargin < 4
            %     v = sim.spherical_mean(samples,w);
            % else
            %     v = sim.spherical_mean(samples,w,est_prev);
            % end
            v = normc(sum(samples.*(ones(3,1)*w),2));
        end
        
        function [C] = spherical_coords(~,X)
            C = [ atan2(X(2,:),X(1,:));
                  acos(X(3,:)) ];
        end
        
        function [C] = cartesian_coords(~,X)
            C = [ cos(X(1,:)).*sin(X(2,:));
                  sin(X(1,:)).*sin(X(2,:));
                  cos(X(2,:)) ];
        end
        
        function q = spherical_mean(sim,X,w,q0)
            if nargin < 3
                w = ones(1,size(X,2))/size(X,2);
            end
            W = (ones(size(X,1),1)*w);
            if nargin < 4
                q0 = normc(sum(X.*W,2));
            end
            step_size = inf;
            q = q0;
            Xl = zeros(size(X));
            step = 1;
            while step_size > 10^(-2)
                dist = 0;
                for i = 1:size(X,2)
                    Xl(:,i) = sim.sphere_log(q,X(:,i));
                    dist = dist + w(i)*real(acos(Xl(:,i).'*q))^2;
                end
                u = (eye(3)-q*q.')*sum(bsxfun(@minus,Xl,q).*W,2);
                q = sim.ExpMap(q,u);
                step_size = norm(u);
                if mod(step,100) == 0
                    disp(step)
                    disp(step_size)
                end
                step = step + 1;
            end
        end
        
        function y = ExpMap(~, x, v, t)
            if nargin < 4
                t = 1;
            end

            tv = t*v;
            nrm_tv = sqrt(sum(tv.^2,1));%norm(tv, 'fro');
            y = bsxfun(@times,x,cos(nrm_tv)) + ...
                bsxfun(@times,tv,sinc(nrm_tv));
            ynrms = sqrt(sum(y.^2,1));
            y = bsxfun(@rdivide, y, ynrms);
        end

        function y = sphere_exponential(~, x, v, t)
            if nargin < 4
                t = 1;
            end

            tv = t*v;
            nrm_tv = sqrt(sum(tv.^2,1));%norm(tv, 'fro');
            y = bsxfun(@times,x,cos(nrm_tv)) + ...
                bsxfun(@times,tv,sinc(nrm_tv));
            ynrms = sqrt(sum(y.^2,1));
            y = bsxfun(@rdivide, y, ynrms);
        end

        function v = sphere_log(~, x1, x2)
            v = (eye(3)-x1*x1.')*(x2 - x1);
            di = real(acos(x1.'*x2));
            % If the two points are "far apart", correct the norm.
            if di > 1e-6
                nv = norm(v, 'fro');
                v = v * (di / nv);
            end
        end
        
    end
    
end

