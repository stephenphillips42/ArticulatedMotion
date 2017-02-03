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
            sim.R_h = definedOrDefault('R_h',eye(sim.sphere_dim));
            sim.T_h = definedOrDefault('T_h',[zeros(sim.sphere_dim-1,1);
                                              sim.sphere_dim+1]);
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
            p = sim.sphere_dim;
            if nargin < 3
                noise = randn(p*sim.n_edges, size(x,2))*sim.sigma_f;
            end
            v = zeros(2*p*sim.n_edges,size(x,2));
            for i = 1:sim.n_edges
                pinds = (1:p) + 2*p*(i-1);
                vinds = (1:p) + p + 2*p*(i-1);
                P = x([pinds,vinds],:);
                V = [ x(vinds,:); noise((1:p)+p*(i-1),:) ];
                v([pinds,vinds],:) = tangent_sphere_exp_map(P,V,1);
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
            l = sum(lognormpdf(d,0,sim.sigma_h),1);
        end

        
        function [P, J] = build_positions_from_state(sim,X)
            % Get 3D positions of the vertecies and joints from state X
            % Input: X - (sphere_dim*n_edges x N) size matrix of the states
            % Output: P - (sphere_dim x n_edges+1 x N) size matrix of the positions
            %             for each of the states
            %         J - (2*sphere_dim x n_edges x N) size matrix of the joints roots
            %             and vectors
            p = sim.sphere_dim;
            xind = @(i) (1:p) + 2*p*(i-2); % Position indeces for graph
            P = zeros(p,sim.n_edges+1,size(X,2));
            J = zeros(2*p,sim.n_edges,size(X,2));
            for i = 1:p
                P(i,:,:) = sim.root_pos(i) + sim.T_h(i);
            end
            for k = 1:sim.n_edges
                i = sim.edges(k,1);
                j = sim.edges(k,2);
                % Compute position (stored for later)
                J(1:p,j-1,:) = P(:,i,:);
                J((1:p)+p,j-1,:) = sim.lengths(j-1)*sim.R_h*X(xind(j),:);
                P(:,j,:) = P(:,i,:) + J((1:p)+p,j-1,:);
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
            p = sim.sphere_dim;
            samples = zeros(2*p*sim.n_edges,sim.n_samples);
            my_x0 = sim.simrun.x_gt(:,1);
            for i = 1:sim.n_edges
                pinds = (1:p) + 2*p*(i-1);
                vinds = (1:p) + p + 2*p*(i-1);
                pos_noise = sim.sigma_ip*randn(p,sim.n_samples);
                vel_noise = sim.sigma_iv*randn(p,sim.n_samples);
                s_pos = mynormc(bsxfun(@plus,my_x0(pinds), pos_noise));
                s_vel = bsxfun(@plus, my_x0(vinds), vel_noise);
                s_vel = s_vel - bsxfun(@times,dot(s_vel,s_pos),s_pos);
                samples([pinds,vinds],:) = [s_pos; s_vel];
            end
            w = sim.uniform();
        end

        function v = estimate(sim,samples,w)
            % So I don't know how to measure 'average' in the case
            p = sim.sphere_dim;
            v = zeros(2*p*sim.n_edges,1);
            for i = 1:sim.n_edges
                pinds = (1:p) + 2*p*(i-1);
                vinds = (1:p) + p + 2*p*(i-1);
                v(pinds) = mynormc(sum(...
                                bsxfun(@times,samples(1:p,:),w),2));
                v(vinds) = sum(bsxfun(@times,samples((1:p)+p,:),w),2);
                v(vinds) = v(vinds) - dot(v(pinds),v(vinds))*v(pinds);
            end
        end
        
        function [e, e_full] = compute_error(sim,x_true,x_est)
            % Position error (degrees)
            p = sim.sphere_dim;
            e1 = zeros(sim.n_edges,size(x_true,2));
            for i = 1:sim.n_edges
                xind = (1:p) + 2*p*(i-1);
                e1(i,:) = acos(dot(mynormc(x_true(xind,:)),...
                                   mynormc(x_est(xind,:))));
            end
            % Velocity percent error
            e2 = zeros(sim.n_edges,size(x_true,2));
            for i = 1:sim.n_edges
                vind = (1:p) + p + 2*p*(i-1);
                true_norms = sqrt(sum(x_true(vind,:).^2,1));
                diff_norms = sqrt(sum((x_true(vind,:) - ...
                                       x_est(vind,:)).^2,1));
                e2(i,:) = diff_norms ./ true_norms;
            end
            e = sum(e1,1);
            e_full = [e1;e2];
        end

        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,~,~)
            subplot(2,1,1)
            % Change view of world
            Zroot = sim.root_pos(3) + sim.T_h(3);
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
        function [samples,w] = resample(sim,samples,w)
        % Resample based on log probabilities w
            [samples, ~] = sim.resample_weighted(samples,w);
            % Add noise...
            p = sim.sphere_dim;
            for i = 1:sim.n_edges
                pinds = (1:p) + 2*p*(i-1);
                vinds = (1:p) + p + 2*p*(i-1);
                pos_noise = sim.sigma_rp*randn(p,sim.n_samples);
                vel_noise = sim.sigma_rv*randn(p,sim.n_samples);
                s_pos = mynormc(samples(pinds,:) + pos_noise);
                s_vel = samples(vinds,:) + vel_noise;
                s_vel = s_vel - bsxfun(@times,dot(s_vel,s_pos),s_pos);
                samples([pinds,vinds],:) = [s_pos; s_vel];
            end
            
            w = sim.uniform();
        end

        function v = Proj(sim,x)
            p = sim.sphere_dim;
            if length(size(x)) == 2
                v = bsxfun(@rdivide,x(1:(p-1),:),x(p,:));
            elseif length(size(x)) == 3
                v = bsxfun(@rdivide,x(1:(p-1),:,:),x(p,:,:));
            else
                error('Input x on projections wrong size');
            end
        end
        
        function w_k = weight_power(sim,w,p)
            w_k = sim.normalize(sim.normalize(w)*p);
        end
    end
end

