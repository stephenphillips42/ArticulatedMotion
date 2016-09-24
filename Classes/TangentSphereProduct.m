classdef TangentSphereProduct < ParticleFilterSim
%TangentSphereProduct Simulation of a particle filter on the a sphere using the
%projection equation for the measurement, with the velocity attached in
%the state

    properties
        %%% Filter Properties
        v        % Average velocity of motion
        sigma_f  % Variance of the noise for movement
        sigma_h  % Variance of the measurements
        R_h      % Rotation of measurement frame
        T_h      % Translation of the measurement frame
        lengths  % Lengths of the arms

        %%% Helping parameters
    end
    
    methods
        function sim = TangentSphereProduct(x0,T)
            % Parameters for parent classes
            nparticles = 2000; % Number of particles
            dim = 12; % Dimension of the space - R^3 x R^3
            meas_dim = 4; 
            sim@ParticleFilterSim(x0,T,nparticles,dim,meas_dim);
            % Measurement model parameters
            sim.sigma_h = sqrt(0.01*(pi/180));
            sim.R_h = eye(3);
            sim.T_h = [0;0;4];
            % Motion model parameters
            sim.sigma_f = sqrt(0.01*(pi/180)); 
            % Arm lengths
            sim.lengths = [2;1];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Motion Model
        function [v, noise] = f(sim,x,noise)
            if nargin < 3
                noise = randn(6, size(x,2))*sim.sigma_f;
            end
            P1 = x(1:6,:);
            V1 = [x(4:6,:);noise(1:3,:)];
            P2 = x(7:12,:);
            V2 = [x(10:12,:);noise(4:6,:)];
            v = [ sim.tangent_sphere_exp_map(P1,V1,1); ...
                  sim.tangent_sphere_exp_map(P2,V2,1) ];
        end
        % Measurement Model - standard projection model
        function [meas, noise] = h(sim,x,noise)
            if nargin < 3
                noise = randn(sim.dim_meas,size(x,2))*sim.sigma_h;
            end
            P1 = x(1:3,:)*sim.lengths(1);
            P2 = x(7:9,:)*sim.lengths(2) + x(1:3,:)*sim.lengths(1);
            X1 = bsxfun(@plus,sim.R_h*P1,sim.T_h);
            X2 = bsxfun(@plus,sim.R_h*P2,sim.T_h);
            meas = [sim.Proj(X1); sim.Proj(X2)] + noise;
        end
        % Measurement Likelihood
        function l = h_likelihood(sim,x,z)
            d = bsxfun(@minus,sim.h(x,0),z);
            l = prod(normpdf(d,0,sim.sigma_h),1);
        end

        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function v = estimate(~,samples,w)
            % This is very VERY much not correct but I'll worry about it later
            pos1 = normc(sum(bsxfun(@times,samples(1:3,:),w),2));
            pos2 = normc(sum(bsxfun(@times,samples(7:9,:),w),2));
            vel1 = (sum(bsxfun(@times,samples(4:6,:),w),2));
            vel2 = (sum(bsxfun(@times,samples(10:12,:),w),2));
            v = [pos1; vel1; pos2; vel2];
        end
        

        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,w,est)
            % Change view of world
            scatter3(0,0,0);
            % [az,el] = view;
            % view(az, el)
            hold on
            % Plot samples
            sim.plot_joint(samples,'b',w)
            % Plot measurement
            for i = 1:length(sim.lengths)
                inds = (1:2) + 2*(i-1);
                start_point = -sim.R_h.'*sim.T_h;
                end_point = sim.R_h.'*([meas(inds);1]*10 - sim.T_h);
                plot3([start_point(1),end_point(1)], ...
                     [start_point(2),end_point(2)], ...
                     [start_point(3),end_point(3)],'m-')
            end
            % Plot estimate
            % sim.plot_joint(est,'c');
            % Plot ground truth
            sim.plot_joint(x_gt,'r');
            hold off
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis equal
            plot_scale = sum(sim.lengths)*1.2;
            axis([-1 1 -1 1 -1 1]*plot_scale)
        end
        
        function plot_joint(sim,X,c,w)
            if nargin < 4
                w = [];
            end
            PosPrev = zeros(3,size(X,2));
            for i = 1:length(sim.lengths)
                inds = (1:3)+6*(i-1);
                Pos = X(inds,:)*sim.lengths(i)+PosPrev;
                Vel = X(inds+3,:)*sim.lengths(i);
                if ~isempty(w)
                    scatter3(Pos(1,:),Pos(2,:),Pos(3,:),eps+w*(2^14),[c '.'])
                else
                    scatter3(Pos(1,:),Pos(2,:),Pos(3,:),[c '.'])
                end
                quiver3(PosPrev(1,:),PosPrev(2,:),PosPrev(3,:),...
                    Pos(1,:)-PosPrev(1,:),...
                    Pos(2,:)-PosPrev(2,:),...
                    Pos(3,:)-PosPrev(3,:),0,c)
                quiver3(Pos(1,:),Pos(2,:),Pos(3,:),...
                        Vel(1,:),Vel(2,:),Vel(3,:),0,c)
                PosPrev = Pos;
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

