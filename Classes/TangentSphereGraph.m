classdef TangentSphereGraph < ParticleFilterSim
%TangentSphereProduct Simulation of a particle filter on the a sphere using the
%projection equation for the measurement, with the velocity attached in
%the state

    properties
        %%% Filter Properties
        sigma_f  % Variance of the noise for movement
        sigma_h  % Variance of the measurements
        sigma_ip % Variance of position sampling for initalization
        sigma_iv % Variance of velocity sampling for initalization
        sigma_rp % Variance of resampling (position)
        sigma_rv % Variance of resampling (velocity)
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
        sphere_dim % Dimension of the underlying space of the tangent bundle
    end
    
    % TODO: Make this so that it uses the PointClasses
    
    methods
        function sim = TangentSphereGraph(varargin)
            % Parameters for parent classes
            sim@ParticleFilterSim(varargin{:});
            definedOrDefault = @(name,default) ...
                         definedOrDefault_long(name,default,varargin);
            % Graph parameters
            sim.edges = definedOrDefault('edges',[]);
            sim.n_edges = size(sim.edges,1);
            sim.lengths = definedOrDefault('lengths',ones(1,sim.n_edges));
            sim.edges_constr = definedOrDefault('edges_constr',[1,1]);
            % Other basic parameters
            sim.sphere_dim = definedOrDefault('sphere_dim',3);
            sim.dim_space = 2*sim.sphere_dim*sim.n_edges;
            sim.root_pos = definedOrDefault('root_pos',zeros(sim.sphere_dim,1));
            assert(size(sim.x0,1) == sim.dim_space)
            sim.n_samples = definedOrDefault('nsamples',5000);
            % Measurement model parameters
            sim.dim_meas = definedOrDefault('dim_meas',(sim.sphere_dim-1)*sim.n_edges);
            sim.dt = definedOrDefault('dt',0.4);
            sim.sigma_h = definedOrDefault('sigma_h',0.001);
            sim.R_h = definedOrDefault('R_h',eye(3));
            sim.T_h = definedOrDefault('T_h',[0;0;4]);
            % Motion model parameters
            sim.sigma_f = definedOrDefault('sigma_f',sqrt(0.001*(pi/180)));
            % Other noise parameters
            sim.sigma_ip = definedOrDefault('sigma_ip',0.1);
            sim.sigma_iv = definedOrDefault('sigma_iv',0.0);
            sim.sigma_rp = definedOrDefault('sigma_rp',0.016);
            sim.sigma_rv = definedOrDefault('sigma_rv',0.001);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Model of the system %%%%%%%%%%%%%%%%%%%%%%%
        % Motion Model
        function [v, noise] = f(sim,x,noise)
            if nargin < 3
                noise = randn(3*sim.n_edges, size(x,2))*sim.sigma_f;
            end
            v = zeros(6*sim.n_edges,size(x,2));
            k = sim.sphere_dim;
            for i = 1:sim.n_edges
                inds = (1:2*k) + 2*k*(i-1);
                P = x(inds,:);
                V = [ x(inds(k+1:2*k),:); noise((1:k)+k*(i-1),:) ];
                v(inds,:) = tangent_sphere_exp_map(P,V,1);
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
            meas = projs(sim.sphere_dim:end,:) + noise; % Get rid of root node
        end
        % Measurement Likelihood
        function l = h_likelihood(sim,x,z)
            d = bsxfun(@minus,sim.h(x,0),z);
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
        
        function [samples,w] = create_samples(sim)
            nlengths = length(sim.lengths);
            samples = zeros(6*nlengths,sim.n_samples);
            my_x0 = sim.simrun.x_gt(:,1);
            for i = 1:nlengths
                inds = (1:6) + 6*(i-1);
                pos_noise = sim.sigma_ip*randn(3,sim.n_samples);
                vel_noise = sim.sigma_iv*randn(3,sim.n_samples);
                s_pos = sim.normc(bsxfun(@plus,my_x0(inds(1:3)), pos_noise));
                s_vel = bsxfun(@plus, my_x0(inds(4:6)), vel_noise);
                s_vel = s_vel - bsxfun(@times,dot(s_vel,s_pos),s_pos);
                samples(inds,:) = [s_pos; s_vel];
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
        function plot_simulation(sim,x_gt,meas,samples,~,~)
            subplot(2,1,1)
            % Change view of world
            Zroot = sim.root_pos(3);
            scatter3(0,0,0);
            hold on
            % Plot samples
            sim.plot_joint_graph(samples,'b',0)
            % Plot measurement
            meas = reshape(meas,2,sim.n_edges);
            msize = size(meas(1:2:end,:)); % One dimension of measurement size
            plot3([ zeros(msize); ones(msize)*(Zroot+2)     ],...
                  [ zeros(msize); meas(1:2:end,:)*(Zroot+2) ],...
                  [ zeros(msize); meas(2:2:end,:)*(Zroot+2) ],'--')
            % Plot estimate
            sim.plot_joint_graph(sim.x0,'c');
            % Plot ground truth
            sim.plot_joint_graph(x_gt,'r');
            hold off
            axis equal
            axis([[-0.1, Zroot+1],[-1,1]*1.3,[-1,1]*3 ])
            xlabel('z')
            ylabel('x')
            zlabel('y')
            subplot(2,1,2)
            
            % quiver(reshape(J(1,:,:),1,[])./reshape(J(3,:,:),1,[]),...
            %        reshape(J(2,:,:),1,[])./reshape(J(3,:,:),1,[]),...
            %        reshape(J(4,:,:),1,[])./reshape(J(6,:,:),1,[]),...
            %        reshape(J(5,:,:),1,[])./reshape(J(6,:,:),1,[]),0,'r');
            [P,~] = sim.build_positions_from_state(samples);
            scatter(reshape(P(1,:,:),1,[])./reshape(P(3,:,:),1,[]),...
                    reshape(P(2,:,:),1,[])./reshape(P(3,:,:),1,[]),10,'b.');
            hold on;
            [P,~] = sim.build_positions_from_state(x_gt);
            scatter(reshape(P(1,:,:),1,[])./reshape(P(3,:,:),1,[]),...
                    reshape(P(2,:,:),1,[])./reshape(P(3,:,:),1,[]),10,'ro');
            hold off;
            axis equal;
        end
        
        function plot_joint_graph(sim,X,c,withquiver)
            if nargin < 4
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
            scatter3(reshape(P(3,:,:),1,[]), ...
                     reshape(P(1,:,:),1,[]), ...
                     reshape(P(2,:,:),1,[]), 100, [c '.']);
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
        % Resample based on log probabilities w
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
            % Add noise...
            for i = 1:sim.n_edges
                inds = (1:6) + 6*(i-1);
                pos_noise = sim.sigma_rp*randn(3,sim.n_samples);
                vel_noise = sim.sigma_rv*randn(3,sim.n_samples);
                s_pos = sim.normc(samples(inds(1:3),:) + pos_noise);
                s_vel = samples(inds(4:6),:) + vel_noise;
                s_vel = s_vel - bsxfun(@times,dot(s_vel,s_pos),s_pos);
                samples(inds,:) = [s_pos; s_vel];
            end
            
            w = sim.uniform();
        end

        function v = Proj(sim,x)
            k = sim.sphere_dim;
            if length(size(x)) == 2
                v = bsxfun(@rdivide,x(1:(k-1),:),x(k,:));
            elseif length(size(x)) == 3
                v = bsxfun(@rdivide,x(1:(k-1),:,:),x(k,:,:));
            else
                error('Input x on projections wrong size');
            end
        end
        
        function x = normc(~,x)
        % Simplified version of normc from MATLAB
            x = bsxfun(@rdivide,x,eps+sqrt(sum(x.^2,1)));
        end

        function p = lognormpdf(~,x,mu,sigma)
        % Dealing with log probability, use this when the normal pdf comes up
            p = -min(log(sqrt(2*pi)*sigma)+(x - mu).^2/(2*sigma.^2),realmax);
        end
        
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
    end
end

