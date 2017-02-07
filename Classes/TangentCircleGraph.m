classdef TangentCircleGraph < TangentSphereGraph
%TangentCircleGraph Simulation of a particle filter on the product of
%circle tangent bundles

    properties
    end
    
    % TODO: Make this so that it uses the PointClasses
    
    methods
        function sim = TangentCircleGraph(varargin)
            % Parameters for parent classes
            varargin = [
                {'sphere_dim',2}, varargin
                ];
            sim@TangentSphereGraph(varargin{:});
            definedOrDefault = @(name,default) ...
                         definedOrDefault_long(name,default,varargin);
            % Measurement model parameters
            sim.sigma_h = definedOrDefault('sigma_h',0.001);
            % Motion model parameters
            sim.sigma_f = definedOrDefault('sigma_f',sqrt(0.1*(pi/180)));
            % Other noise parameters
            sim.sigma_ip = definedOrDefault('sigma_ip',0.1);
            sim.sigma_iv = definedOrDefault('sigma_iv',0.0);
            sim.sigma_rp = definedOrDefault('sigma_rp',0.016);
            sim.sigma_rv = definedOrDefault('sigma_rv',0.001);
        end
        
        %%%%%%%%%%%%%%%%%%% Particle filter functions %%%%%%%%%%%%%%%%%%%%
        function [samples,w] = create_samples(sim)
            p = sim.sphere_dim;
            samples = zeros(2*p*sim.n_edges,sim.n_samples);
            % my_x0 = sim.simrun.x_gt(:,1);
            meas0 = sim.simrun.meas(:,1);
            % TODO: Make this not a loop
            pinds = @(k) (1:p) + 2*p*(k-2);
            minds = @(k) (1:(p-1)) + (p-1)*(k-2);
            for n = 1:sim.n_samples
                s = zeros(2*p*sim.n_edges,1);
                positions = zeros(p,(sim.n_edges + 1));
                positions(:,1) = sim.R_h*sim.root_pos + sim.T_h;
                for k = 1:sim.n_edges
                    e1 = sim.edges(k,1);
                    e2 = sim.edges(k,2);
                    r = positions(:,e1);
                    v = [meas0(minds(e2));1];
                    disc = (r'*v)^2 - (v'*v)*(r'*r - sim.lengths(e2-1)^2);
                    t = ((r'*v) + (disc > 0)*sign(randn)*sqrt(disc))/(v'*v);
                    s(pinds(e2)) = normc(t*v - r);
                    positions(:,e2) = r + sim.lengths(e2-1)*s(pinds(e2));
                end
                samples(:,n) = s;
            end
            samples = sim.add_noise(samples);
            w = sim.uniform();
        end
        
        % TODO: I like this plotting... but I don't know what to do with it
        function [samples,w] = create_samples_plot(sim)
            p = sim.sphere_dim;
            samples = zeros(2*p*sim.n_edges,sim.n_samples);
            % my_x0 = sim.simrun.x_gt(:,1);
            meas0 = sim.simrun.meas(:,1);
            % TODO: Make this not a loop
            pinds = @(k) (1:p) + 2*p*(k-2);
            minds = @(k) (1:(p-1)) + (p-1)*(k-2);
            % vinds = @(k) (1:p) + p + 2*p*(k-1);
            plt = true;
            state = sim.build_positions_from_state(sim.x0);
            for n = 1:sim.n_samples
                s = zeros(2*p*sim.n_edges,1);
                positions = zeros(p,(sim.n_edges + 1));
                positions(:,1) = sim.R_h*sim.root_pos + sim.T_h;
                if plt
                    figure(1);
                    subplot(2,1,1)
                    scatter(0,0);
                    hold on;
                    scatter(state(1,:),state(2,:));
                    axis equal;
                    subplot(2,1,2)
                    scatter(0,0);
                end
                for k = 1:sim.n_edges
                    e1 = sim.edges(k,1);
                    e2 = sim.edges(k,2);
                    r = positions(:,e1);
                    v = [meas0(minds(e2));1];
                    disc = (r'*v)^2 - (v'*v)*(r'*r - sim.lengths(e2-1)^2);
                    t = ((r'*v) + (disc > 0)*sign(randn)*sqrt(disc))/(v'*v);
                    if plt
                        subplot(2,1,1)
                        hold on
                        t1 = ((r'*v) + sqrt(disc))/(v'*v);
                        t2 = ((r'*v) - sqrt(disc))/(v'*v);
                        t3 = (r'*v)/(v'*v);
                        scatter([t1,t2]*v(1),[t1,t2]*v(2),'rx');
                        scatter(t3*v(1),t3*v(2),'mx');
                        plot([positions(1,e1),t*v(1)],[positions(2,e1),t*v(2)],'k*-')
                        plot([0,5*meas0(e2-1)],[0,5])
                        plot(sim.lengths(e2-1)*cos(0:0.01:2*pi) + r(1),sim.lengths(e2-1)*sin(0:0.01:2*pi) + r(2))
                        hold off
                    end
                    s(pinds(e2)) = normc(t*v - r);
                    positions(:,e2) = r + sim.lengths(e2-1)*s(pinds(e2));
                    if plt
                        subplot(2,1,2)
                        z = normc(t*v - r);
                        hold on
                        plot([0, z(1)], [0, z(2)])
                        axis equal
                        xlim([-1 1])
                        ylim([-1 1])
                        hold off
                    end
                end
                samples(:,n) = s;
            end
            samples = sim.add_noise(samples);
            w = sim.uniform();
        end
        
        %%%%%%%%%%%%%%%%%%%% Visualization functions %%%%%%%%%%%%%%%%%%%%%%
        function plot_simulation(sim,x_gt,meas,samples,~,est)
            % Change view of world
            Zroot = sim.root_pos(2) + sim.T_h(2);
            scatter(0,0);
            hold on
            % Plot samples
            sim.plot_joint_graph(samples,'b',0)
            % Plot measurement
            meas = reshape(meas,1,sim.n_edges);
            msize = size(meas); % One dimension of measurement size
            plot([ zeros(msize);        meas*(Zroot+2) ],...
                 [ zeros(msize); ones(msize)*(Zroot+2) ],'--')
            % Plot estimate
            sim.plot_joint_graph(est,'c');
            % Plot ground truth
            sim.plot_joint_graph(x_gt,'r');
            % Plot measurement stuff
            plot([-10,10],[1,1])
            [P,~] = sim.build_positions_from_state(samples);
            sample_meas = reshape(P(1,:,:)./P(2,:,:),1,[]);
            scatter(sample_meas,ones(size(sample_meas)),'gx')
            hold off
            axis equal
            axis([[-1,1]*1.3,[-0.1, Zroot+1] ])
            xlabel('x')
            ylabel('y')
        end
        
        function plot_joint_graph(sim,X,c,withquiver)
            if nargin < 4
                withquiver = true;
            end
            % Drawing shapes
            [P,J] = sim.build_positions_from_state(X);
            if withquiver
                quiver(reshape(J(1,:,:),1,[]),...
                       reshape(J(2,:,:),1,[]),...
                       reshape(J(3,:,:),1,[]),...
                       reshape(J(4,:,:),1,[]),0,c);
            end
            scatter(reshape(P(1,:,:),1,[]), ...
                    reshape(P(2,:,:),1,[]), 100, [c '.']);
        end
    end
end

