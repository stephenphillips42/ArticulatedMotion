classdef Sphere < ProjectionParticleFilterSim
%SphereParticleFilterSim Simulation of a particle filter on the
%a sphere using the projection equation for the measurement

    properties
        %%% Filter Properties
        
        v        % Average velocity of motion
        sigma_f  % Variance of the noise for movement
        
        %%% Plotting parameters
        
        plot_type % Plot 2D or 3D?
        X % X coordinates of the sphere surface
        Y % Y coordinates of the sphere surface
        Z % Z coordinates of the sphere surface
        angle_offset % In 3D view, where we are in spin around the sphere
        angle_speed % In 3D view, how fast to spin around the sphere
        
        %%% Helping parameters
        
        omega = [0; 0; 1];
        omegahat
    end
    
    methods
        function sim = Sphere(x0,T,plot_type)
            if nargin < 3
                plot_type = true;
            end
            % Parameters for parent classes
            nparticles = 1000; % Number of particles
            dim = 3; % Dimension of the space - R^3
            sigma_h = sqrt(0.02*(pi/180));
            sim@ProjectionParticleFilterSim(x0,T,nparticles,dim,sigma_h);
            % Motion model parameters
            sim.omegahat = sim.hat(sim.omega);
            sim.sigma_f = sqrt(1*(pi/180)); 
            sim.plot_type = plot_type;
            
            % Visualization Variables
            [sim.X, sim.Y, sim.Z] = sphere(20);
            sim.angle_offset = 0;
            sim.angle_speed = 5;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Motion Model
        function [v, noise] = f(sim,x,noise)
            if nargin < 3
                noise = randn(size(x))*sim.sigma_f;
                noise = noise - bsxfun(@times,x,dot(x,noise));
            end
            vel = sim.dt*sim.omegahat*x + noise;
            v = sim.sphere_exp_map(x,vel,1);
        end

        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function [samples,w] = create_samples(sim)
            samples = normc(randn(sim.dim_space,sim.n_samples));
            w = sim.uniform();
        end
        
        function v = estimate(~,samples,w)
            v = normc(sum(bsxfun(@times,samples,w),2));
        end

        function [e, e_full] = compute_error(~,x_true,x_est)
            e = acos(dot(normc(x_true),normc(x_est)));
            e_full = e;
        end
        
        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,w,est)
            if sim.plot_type
                sim.plot_simulation_3D(x_gt,meas,samples,w,est)
            else
                sim.plot_simulation_2D(x_gt,meas,samples,w,est)
            end
        end

        function plot_simulation_3D(sim,x_gt,meas,samples,w,est)
            % Plot world
            colormap summer;
            C = ones(size(sim.X));
            C(:) = sim.h_likelihood([sim.X(:) sim.Y(:) sim.Z(:)].',meas);
            surf(sim.X,sim.Y,sim.Z,C);
            % Change view of world
            [az,el] = view;
            sim.angle_offset = sim.angle_offset + sim.angle_speed;
            view(az + sim.angle_offset, el)
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
            sim.plot_meas_2D(meas)
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
                scatter([m1(1),m2(1)],[m1(2),m2(2)],400,'rx')
            end
        end


        %%%%%%%%%%%%%%%% Differential Geometry functions %%%%%%%%%%%%%%%%%%%
        function y = sphere_exp_map(~,x,v,t)
            v = t*v;
            nrm_tv = sqrt(sum(v.^2,1));%norm(tv, 'fro');
            y = bsxfun(@times,x,cos(nrm_tv)) + ...
                bsxfun(@times,v,sinc(nrm_tv));
            ynrms = sqrt(sum(y.^2,1));
            y = bsxfun(@rdivide, y, ynrms);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%
        function xhat = hat(~,x)
            xhat = [     0, -x(3),  x(2); 
                      x(3),     0, -x(1);
                     -x(2),  x(1),     0 ];
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

        
    end
    
end

