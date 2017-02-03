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

