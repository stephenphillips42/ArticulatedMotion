classdef TangentSphereLength < ParticleFilterSim
%TangentSphereLength Simulation of a particle filter on the a sphere using the
%projection equation for the measurement, with the velocity attached in
%the state, as well as the scale of the sphere

    properties
        %%% Filter Properties
        v        % Average velocity of motion
        sigma_f  % Variance of the noise for movement
        sigma_h  % Variance of the measurements
        R_h      % Rotation of measurement frame
        T_h      % Translation of the measurement frame
        
        %%% Plotting parameters
        
        plot_type % Plot 2D or 3D?
        X % X coordinates of the sphere surface
        Y % Y coordinates of the sphere surface
        Z % Z coordinates of the sphere surface
        angle_offset % In 3D view, where we are in spin around the sphere
        angle_speed % In 3D view, how fast to spin around the sphere

        %%% Helping parameters
    end
    
    methods
        function sim = TangentSphereLength(x0,T,plot_type)
            if nargin < 3
                plot_type = true;
            end
            % Parameters for parent classes
            nparticles = 2000; % Number of particles
            dim = 7; % Dimension of the space - R^3 x R^3
            meas_dim = 2; % Only 2D perspective projection
            sim@ParticleFilterSim(x0,T,nparticles,dim,meas_dim);
            % Measurement model parameters
            sim.sigma_h = sqrt(0.05*(pi/180));
            sim.R_h = eye(3);
            sim.T_h = [0;0;4];
            % Motion model parameters
            sim.sigma_f = sqrt(0.05*(pi/180)); 
            % Visualization Variables
            sim.plot_type = plot_type;
            [sim.X, sim.Y, sim.Z] = sphere(20);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Motion Model
        function [v, noise] = f(sim,x,noise)
            if nargin < 3
                noise = randn(size(x(4:7,:)))*sim.sigma_f;
            end
            v = [sim.tangent_sphere_exp_map(x(1:6,:),[x(4:6,:);noise(1:3,:)],1); ...
                 x(7,:) + noise(4,:)];
        end
        % Measurement Model - standard projection model
        function [meas, noise] = h(sim,x,noise)
            if nargin < 3
                noise = randn(sim.dim_meas,size(x,2))*sim.sigma_h;
            end
            X_scale = bsxfun(@times,x(1:3,:),x(7,:));
            meas = sim.Proj(bsxfun(@plus,sim.R_h*X_scale,sim.T_h)) + noise;
        end
        % Measurement Likelihood
        function l = h_likelihood(sim,x,z)
            d = bsxfun(@minus,sim.h(x,0),z);
            l = prod(normpdf(d,0,sim.sigma_h),1);
        end
        
        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function [samples,w] = create_samples(sim)
            s_pos = normc(randn(3,sim.n_samples));
            s_vel = randn(3,sim.n_samples);
            s_vel = s_vel - bsxfun(@times,dot(s_vel,s_pos),s_pos);
            s_scale = 0.5+1*rand(1,sim.n_samples);
            samples = [s_pos; s_vel; s_scale];
            w = sim.l1normalize(ones(1,length(samples)));
        end

        function v = estimate(~,samples,w)
            pos = normc(sum(bsxfun(@times,samples(1:3,:),w),2));
            % This is very much not correct but I'll worry about it later
            vel = (sum(bsxfun(@times,samples(4:6,:),w),2));
            scale = (sum(bsxfun(@times,samples(7,:),w),2));
            v = [pos; vel; scale];
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
%             colormap summer;
%             C = ones(size(sim.X));
%             C(:) = sim.h_likelihood([sim.X(:) sim.Y(:) sim.Z(:) zeros(length(sim.Z(:)),4)].', meas);
%             surf(sim.X,sim.Y,sim.Z,C);
            % Change view of world
            [az,el] = view;
            % [az_gt,el_gt,~] = cart2sph(x_gt(1),x_gt(2),x_gt(3));
            % disp([az,el,az_gt,el_gt,az_gt*(180/pi)])
            view(az, el)
            axis equal
            % Plot samples
            scatter3(samples(1,:).*samples(7,:),...
                     samples(2,:).*samples(7,:),...
                     samples(3,:).*samples(7,:),...
                     eps+w*(2^14),'b.')
            hold on
            quiver3(samples(1,:).*samples(7,:),...
                    samples(2,:).*samples(7,:),...
                    samples(3,:).*samples(7,:),...
                    samples(4,:).*samples(7,:),...
                    samples(5,:).*samples(7,:),...
                    samples(6,:).*samples(7,:),0,'b')
            % Plot measurement
            start_point = -sim.R_h.'*sim.T_h;
            end_point = sim.R_h.'*([meas;1]*10 - sim.T_h);
            plot3([start_point(1),end_point(1)], ...
                 [start_point(2),end_point(2)], ...
                 [start_point(3),end_point(3)],'r-')
            % Plot estimate
            scatter3(est(1)*est(7),est(2)*est(7),est(3)*est(7),500,'c.')
            quiver3(est(1)*est(7),est(2)*est(7),est(3)*est(7),...
                    est(4)*est(7),est(5)*est(7),est(6)*est(7),0,'c')
            % Plot ground truth
            scatter3(x_gt(1)*x_gt(7),x_gt(2)*x_gt(7),x_gt(3)*x_gt(7),500,'r.')
            quiver3(x_gt(1)*x_gt(7),x_gt(2)*x_gt(7),x_gt(3)*x_gt(7),...
                    x_gt(4)*x_gt(7),x_gt(5)*x_gt(7),x_gt(6)*x_gt(7),0,'r')
            axis equal
            axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])
            hold off
        end
        
        function plot_simulation_2D(sim,x_gt,meas,samples,w,est)
            % Plot likelihood
            [a,b] = meshgrid(-pi:0.01:pi,0:0.01:pi);
            z = zeros(size(a));
            z(:) = sim.h_likelihood([sim.cartesian_coords([a(:) b(:)].'); zeros(4,length(a(:)))],meas);
            figure(1)
            colormap summer
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
            s = -sim.R_h(1:3,1:3).'*sim.T_h(1:3); % Start point
            e = sim.R_h(1:3,1:3).'*([meas;1]*10 - sim.T_h(1:3)); % End point e
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
            % Exponential map of the sphere
            v = t*v;
            nrm_tv = sqrt(sum(v.^2,1));%norm(tv, 'fro');
            y = bsxfun(@times,x,cos(nrm_tv)) + ...
                bsxfun(@times,v,sinc(nrm_tv));
            y = normc(y);
        end
        
        function r = curvature_r(~,u,v,w)
            % Riemannian curvature for sphere
            % r = (w'*u)*v-(w'*v)*u;
            r = bsxfun(@times,dot(w,u),v)-bsxfun(@times,dot(w,v),u);
        end

        function zetat = parallel_transport(sim,x,xt,zeta0)
            % Parallel transport of the sphere
            x = normc(x);
            xt = normc(xt);

            % Make sure zeta0 is in the tangent space of x
            zeta0 = zeta0 - bsxfun(@times,x,dot(x,zeta0));

            % Transport zeta0 from the tangent space of x to the tangent space of xt
            xi =  sim.log_map_sphere(x,xt);
            if ~all(isreal(xi))
                disp('ohno parallel_transport')
            end
            nrm_xi = sqrt(sum(xi.^2,1));
            u = normc(xi);
            u_dot_zeta = dot(u,zeta0);
            zetat = -bsxfun(@times,x,sin(nrm_xi).*u_dot_zeta) + ...
                     bsxfun(@times,u,cos(nrm_xi).*u_dot_zeta) + ...
                     zeta0 - bsxfun(@times,u,u_dot_zeta);   
        end
        
        function mlog = log_map_sphere(~,x,y)
            % Log map of the sphere

            % Get geodesic distance between x and y
            trxy = dot(x,y);
            geodist = acos(trxy) ;
            if ~all(isreal(geodist))
                disp('ohno log_map_sphere')
                disp(geodist)
            end

            % TODO Is there a better way of writing the log map?
            mlog = bsxfun(@times, geodist, bsxfun(@rdivide,...
                            y-bsxfun(@times,x,trxy),...
                            sqrt(1-trxy.^2)));
            % if geodist< eps
            %     mlog = zeros(size(x));
            % else
            %     mlog = ((y-x*trxy)/sqrt(1-trxy^2)) *geodist;
            % end
        end
        
        % function [p,u,v,w] = ExpMapTS(p,u,v,w,L)
        
        function [y,z] = tangent_sphere_exp_map(sim,x,v,t,L)
            % Numerically computes exponential map of tangent bundle of
            % the sphere
            if nargin < 4
                t = 1;
            end
            if nargin<5
                L = 50;
            end
            % Rename different parts of the exponential map
            p = x(1:3,:);
            u = x(4:6,:);
            w = t*v(4:6,:);
            v = t*v(1:3,:);
            
            % How many iterations
            stepsize = 1/L;

            for i=1:L
                 pprev = p;
                 uprev = u;
                 vprev = v;
                 wprev = w;

                 p = sim.sphere_exp_map(p,v,stepsize);
                 u = sim.parallel_transport(pprev,p,uprev+stepsize*wprev);
                 v = sim.parallel_transport(pprev,p,vprev-stepsize*sim.curvature_r(uprev,wprev,vprev));
                 w = sim.parallel_transport(pprev,p,wprev);
            end
            y = [p;u];
            z = [v;w]; % Don't think this is necessary
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
        
        function v = Proj(~,x)
            v = bsxfun(@rdivide,x(1:2,:),x(3,:));
        end

        
    end
    
end

