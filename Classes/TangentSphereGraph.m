classdef TangentSphereGraph < ParticleFilterSim
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
        root_pos % Position of the root node (assumed known)

        %%% Graph Properties
        n_edges      % Number of joints to consider
        lengths      % Joint lengths of connectivity graph
        edges        % Edges of the graph - should be a joint graph (out-tree)
                     % Assumes that root vertex is labeled 1
        edges_constr % Additional edges - the non-tree edges (NOT IMPLEMENTED)
        
        %%% Helping parameters
    end
    
    % TODO: Make this so that it uses the PointClasses
    
    methods
        function sim = TangentSphereGraph(x0,T,L,E,Ed,root_pos)
            if nargin < 6
                root_pos = [0;0;0];
            end
            % Parameters for parent classes
            nparticles = 10000; % Number of particles
            dim = 12; % Dimension of the space - R^3 x R^3
            meas_dim = length(L)*2;
            sim@ParticleFilterSim(x0,T,nparticles,dim,meas_dim);
            % Measurement model parameters
            sim.sigma_h = 0.005;
            sim.R_h = eye(3);
            sim.T_h = [0;0;4];
            % Motion model parameters
            sim.sigma_f = sqrt(0.005*(pi/180));
            % Graph parameters
            sim.n_edges = size(E,1);
            sim.lengths = L;
            sim.edges = E;
            sim.edges_constr = Ed;
            sim.root_pos = root_pos;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Motion Model
        function [v, noise] = f(sim,x,noise)
            if nargin < 3
                noise = randn(3*sim.n_edges, size(x,2))*sim.sigma_f;
            end
            v = zeros(6*sim.n_edges,size(x,2));
            for i = 1:sim.n_edges
                inds = (1:6) + 6*(i-1);
                P = x(inds,:);
                V = [ x(inds(4:6),:); noise((1:3)+3*(i-1),:) ];
                v(inds,:) = sim.tangent_sphere_exp_map(P,V,1);
            end
        end
        % Measurement Model
        % Assumes that vertex 1 is root and its position is known
        function [meas, noise] = h(sim,x,noise)
            if nargin < 3
                noise = randn(sim.dim_meas,size(x,2))*sim.sigma_h;
            end
            % Get joint positions
            [P,~] = sim.build_positions_from_state(x);
            % Make projections and noise
            projs = reshape(sim.Proj(P),[],size(P,3));
            meas = projs(3:end,:) + noise; % Get rid of root node
        end
        % Measurement Likelihood
        function l = h_likelihood(sim,x,z)
            d = bsxfun(@minus,sim.h(x,0),z);
            % l = prod(normpdf(d,0,sim.sigma_h+0.0001),1)+eps;
            % Use log probabilities due to numerical instabiity
            l = sum(sim.lognormpdf(d,0,sim.sigma_h),1);
        end

        
        function [P, J] = build_positions_from_state(sim,X)
            % Get 3D positions of the vertecies and joints from state X
            % Input: X - (3*n_edges x N) size matrix of the states
            % Output: P - (3 x n_edges+1 x N) size matrix of the positions
            %             for each of the states
            %         J - (6 x n_edges x N) size matrix of the joints roots
            %             and vectors
            xind = @(i) (1:3) + 6*(i-2); % Position indeces for graph
            P = zeros(3,sim.n_edges+1,size(X,2));
            J = zeros(6,sim.n_edges,size(X,2));
            P(1,:,:) = sim.root_pos(1);
            P(2,:,:) = sim.root_pos(2);
            P(3,:,:) = sim.root_pos(3);
            for k = 1:sim.n_edges
                i = sim.edges(k,1);
                j = sim.edges(k,2);
                % Compute position (stored for later)
                J(1:3,j-1,:) = P(:,i,:);
                J(4:6,j-1,:) = sim.lengths(j-1)*X(xind(j),:);
                P(:,j,:) = P(:,i,:) + J(4:6,j-1,:);
            end
        end

        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function [samples,w] = create_samples(sim)
            nlengths = length(sim.lengths);
            samples = zeros(6*nlengths,sim.n_samples);
            my_x0 = sim.simrun.x_gt(:,1);
            for i = 1:nlengths
                inds = (1:6) + 6*(i-1);
                pos_noise = 0.2*randn(3,sim.n_samples);
                % vel_noise = 0.1*randn(3,sim.n_samples);
                s_pos = sim.normc(bsxfun(@plus,my_x0(inds(1:3)), pos_noise));
                % s_vel = bsxfun(@plus, my_x0(inds(4:6)), vel_noise);
                % s_vel = s_vel - bsxfun(@times,dot(s_vel,s_pos),s_pos);
                samples(inds,:) = [s_pos; 0*randn(3,sim.n_samples)];
            end
            w = sim.uniform();
        end

        function v = estimate(sim,samples,w)
            % So I don't know how to measure 'average' in the case
            v = zeros(6*sim.n_edges,1);
            for i = 1:sim.n_edges
                v((1:3)+6*(i-1)) = ...
                    sim.normc(sum(bsxfun(@times,samples((1:3)+6*(i-1),:),w),2));
                v((4:6)+6*(i-1)) = ...
                    (sum(bsxfun(@times,samples((4:6)+6*(i-1),:),w),2));
                v((4:6)+6*(i-1)) = v((4:6)+6*(i-1)) ...
                    - dot(v((1:3)+6*(i-1)),v((4:6)+6*(i-1)))*v((1:3)+6*(i-1));
            end
        end
        
        function [e, e_full] = compute_error(sim,x_true,x_est)
            % Position error (degrees)
            e1 = zeros(sim.n_edges,size(x_true,2));
            for i = 1:sim.n_edges
                e1(i,:) = acos(dot(normc(x_true((1:3)+6*(i-1),:)),...
                                   normc(x_est((1:3)+6*(i-1),:))));
            end
            % Velocity percent error
            e2 = zeros(sim.n_edges,size(x_true,2));
            for i = 1:sim.n_edges
                true_norms = sqrt(sum(x_true((4:6)+6*(i-1),:).^2,1));
                diff_norms = sqrt(sum((x_true((4:6)+6*(i-1),:) - ...
                                       x_est((4:6)+6*(i-1),:)).^2,1));
                e2(i,:) = diff_norms ./ true_norms;
            end
            e = sum(e1,1);
            e_full = [e1;e2];
        end

        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,w,est)
            % Change view of world
            Zroot = sim.root_pos(3);
            scatter3(0,0,0);
            % [az,el] = view;
            % view(az, el)
            hold on
            % Plot samples
            sim.plot_joint_graph(samples,'b',w,0)
            % Plot measurement
            meas = reshape(meas,2,sim.n_edges);
            msize = size(meas(1:2:end,:)); % One dimension of measurement size
            plot3([ zeros(msize); ones(msize)*(Zroot+2)     ],...
                  [ zeros(msize); meas(1:2:end,:)*(Zroot+2) ],...
                  [ zeros(msize); meas(2:2:end,:)*(Zroot+2) ],'--')
            % Plot estimate
            % sim.plot_joint_graph(est,'c');
            sim.plot_joint_graph(sim.x0,'c');
            % Plot ground truth
            sim.plot_joint_graph(x_gt,'r');
            hold off
            axis equal
            axis([[-0.1, Zroot+1],[-1,1]*1.3,[-1,1]*3 ])
            xlabel('z')
            ylabel('x')
            zlabel('y')
        end
        
        function plot_joint_graph(sim,X,c,w,withquiver)
            if nargin < 4
                w = [];
            end
            if nargin < 5
                withquiver = true;
            end
            % Drawing shapes
            [P,J] = sim.build_positions_from_state(X);
            if withquiver
                quiver3(reshape(J(3,:,:),1,[]),...
                        reshape(J(1,:,:),1,[]),...
                        reshape(J(2,:,:),1,[]),...
                        reshape(J(6,:,:),1,[]),...
                        reshape(J(4,:,:),1,[]),...
                        reshape(J(5,:,:),1,[]),0,c);
            end
            % if ~isempty(w)
            %     scatter3(reshape(P(3,:,:),1,[]), ...
            %              reshape(P(1,:,:),1,[]), ...
            %              reshape(P(2,:,:),1,[]),...
            %              200*reshape(ones(size(P,2),1)*exp(w),1,[]),[c '.']);
            % else
                scatter3(reshape(P(3,:,:),1,[]), ...
                         reshape(P(1,:,:),1,[]), ...
                         reshape(P(2,:,:),1,[]), 100, [c '.']);
            % end
        end

        %%%%%%%%%%%%%%%% Differential Geometry functions %%%%%%%%%%%%%%%%%%%
        function y = sphere_exp_map(sim,x,v,t)
            % Exponential map of the sphere
            v = t*v;
            nrm_tv = sqrt(sum(v.^2,1));%norm(tv, 'fro');
            y = bsxfun(@times,x,cos(nrm_tv)) + ...
                bsxfun(@times,v,sinc(nrm_tv));
            y = sim.normc(y);
        end
        
        function r = curvature_r(~,u,v,w)
            % Riemannian curvature for sphere
            % r = (w'*u)*v-(w'*v)*u;
            r = bsxfun(@times,dot(w,u),v)-bsxfun(@times,dot(w,v),u);
        end

        function zetat = parallel_transport(sim,x,xt,zeta0)
            % Parallel transport of the sphere
            x = sim.normc(x);
            xt = sim.normc(xt);

            % Make sure zeta0 is in the tangent space of x
            zeta0 = zeta0 - bsxfun(@times,x,dot(x,zeta0));

            % Initialize zetat
            zetat = zeta0;
            % Transport zeta0 from the tangent space of x to the tangent space of xt
            xi =  sim.log_map_sphere(x,xt);
            nrm_xi = sqrt(sum(xi.^2,1));
            valid_inds = (nrm_xi > 10^-9); % TODO Make threshold not magic number
            if any(valid_inds)
                nrm_xi = nrm_xi(valid_inds);
                u = sim.normc(xi(:,valid_inds));
                u_dot_zeta = dot(u,zeta0(:,valid_inds));
                zetat(:,valid_inds) = ...
                        -bsxfun(@times,x,sin(nrm_xi).*u_dot_zeta) + ...
                         bsxfun(@times,u,cos(nrm_xi).*u_dot_zeta) + ...
                         zeta0 - bsxfun(@times,u,u_dot_zeta);   
            end
        end
        
        function mlog = log_map_sphere(~,x,y)
            % Log map of the sphere

            % Get geodesic distance between x and y
            trxy = max(-1,min(dot(x,y),1));
            geodist = acos(trxy) ;

            % TODO Is there a better way of writing the log map?
            mlog = bsxfun(@times, geodist, bsxfun(@rdivide,...
                            y-bsxfun(@times,x,trxy),...
                            eps+sqrt(1-trxy.^2)));
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
        function w = normalize(sim,w)
            % Working in log space
            W = sim.log_cum_sum(w);
            w = w - W(end);
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

        function v = Proj(~,x)
            if length(size(x)) == 2
                v = bsxfun(@rdivide,x(1:2,:),x(3,:));
            elseif length(size(x)) == 3
                v = bsxfun(@rdivide,x(1:2,:,:),x(3,:,:));
            else
                error('Input x on projections wrong size');
            end
        end
        
        function x = normc(~,x)
            x = bsxfun(@rdivide,x,eps+sqrt(sum(x.^2,1)));
        end

        % Dealing with log formats (TODO make this better somehow)
        function p = lognormpdf(~,x,mu,sigma)
            p = -min(log(sqrt(2*pi)*sigma)+(x - mu).^2/(2*sigma.^2),realmax);
        end
        
        % Thank you stack overflow
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

