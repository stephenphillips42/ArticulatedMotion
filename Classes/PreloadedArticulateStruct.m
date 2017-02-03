classdef PreloadedArticulateStruct < StructureSim
    %PreloadedArticulateStruct - Loading an articulated structure from file
    %   Loads an articulated structure from file and build the appropriate
    %   data augmentation structures from file. Expects a directed acylcic
    %   structure of the graph on the joints.
    
    properties
        % Position properties
        X        % Time organized aggregate data (n_steps x 1 cell)
        root_pos % Position of root node
        % Graph properties
        edges   % Edges of the graph, N x 2 matrix
        n_edges % Number of edges
        lengths % Lengths of each articulated joint, N x 1 matrix
        % Camera properties
        R_c % Rotation of camera
        T_c % Translation of camera
        % Constants
        n_steps = 50  % Number of timesteps saves
        deg_step = 0  % Angle rotated per step
        scale = 6     % Scaling of articulated structure
    end
    
    methods
        function s = PreloadedArticulateStruct(matname,jointsname,rootpos,R_c,T_c)
            s@StructureSim(1);
            if nargin < 3
                rootpos = [0;0;0];
            end
            if nargin < 5
                R_c = eye(3);
                T_c = [0;0;4];
            end
            s.root_pos = rootpos;
            s.R_c = R_c;
            s.T_c = T_c;
            % Load data
            jointsld = load(jointsname);
            structld = load(matname);
            s.buildStructure(jointsld,structld);
            s.time_step = 0;
            s.structure = s.X{1};
        end
        function [] = buildStructure(s,jointsdata,structdata)
            % Create the structure of the articulated object from the
            % loaded data jointsdata and structdata
            
            % Build edges
            skel = jointsdata.skel;
            mocapIndices = jointsdata.mocapIndices;
            s.edges = s.build_edges(skel);
            s.n_edges = size(s.edges,1);
            % Build structure
            Sl = structdata.S(:,mocapIndices)/s.scale;
            % Construct the verticies
            sind = @(i) (1:3) + 3*(i-1); % Position indeces for graph
            for t=1:s.n_steps+1
                % Mean center the data
                means = mean(Sl(sind(t),:),2);
                Sl(sind(t),:) = Sl(sind(t),:) - means*ones(1,size(Sl,2));
                % Rotate properly to get extra movement
                Ry = expm(hat([0,1,0])*(t*s.deg_step/180*pi));
                Sl(sind(t),:) = Ry*Sl(sind(t),:);
            end
            pinds = @(i) (1:3) + 6*(i-1); % Position indecies for the state
            vinds = @(i) (4:6) + 6*(i-1); % Velocity indecies for the state
            s.lengths = zeros(s.n_edges,1);
            for t = 1:s.n_steps
                s.X{t} = zeros(6*s.n_edges,1);
                Pm = Sl(sind(t),:);
                Pm1 = Sl(sind(min(s.n_steps,t+1)),:);
                for k = 1:s.n_edges
                    i = s.edges(k,1);
                    j = s.edges(k,2);
                    s.X{t}(pinds(j-1)) = normc(Pm(:,j) - Pm(:,i));
                    x1 = normc(Pm1(:,j) - Pm1(:,i));
                    % Get aggregate lengths
                    s.lengths(k) = s.lengths(k) + norm(Pm(:,i)-Pm(:,j));
                    % Use sphere log-map to get velocity of this time step
                    % Log map of the sphere 
                    trxy = s.X{t}(pinds(j-1)).'*x1;
                    geodist = acos(trxy);
                    s.X{t}(vinds(j-1)) = ...
                        geodist*(x1-s.X{t}(pinds(j-1))*trxy)/sqrt(1-trxy.^2);
                end
            end
            s.lengths = s.lengths/s.n_steps;
        end
        % Model functions
        function motion(s)
            s.time_step = s.time_step + 1; % Take a time step
            s.structure = s.data.x{min(s.time_step,length(s.data.x))};
        end
        function m = meas(s,noise)
            if nargin < 3
                noise = randn(s.dim_meas,size(x,2))*s.sigma_h;
            end
            [P, ~] = s.build_positions();
            projs = reshape(s.proj(P),[],size(P,3));
            m = projs(3:end,:) + noise; % Get rid of root node
        end
        % Positions from state
        function [P, J] = build_positions(s)
            % Get 3D positions of the vertecies and joints from state X
            % Input: X - (3*n_edges x 1) size matrix of the states
            % Output: P - (3 x n_edges+1) size matrix of the positions
            %             for each of the states
            %         J - (6 x n_edges) size matrix of the joints roots
            %             and vectors
            xind = @(i) (1:3) + 6*(i-2); % Position indeces for graph
            P = zeros(3,s.n_edges+1);
            J = zeros(6,s.n_edges);
            P(1,:) = s.root_pos(1);
            P(2,:) = s.root_pos(2);
            P(3,:) = s.root_pos(3);
            for k = 1:s.n_edges
                i = s.edges(k,1);
                j = s.edges(k,2);
                % Compute position (stored for later)
                J(1:3,j-1) = P(:,i);
                J(4:6,j-1) = s.lengths(j-1)*s.structure(xind(j));
                P(:,j) = P(:,i) + J(4:6,j-1);
            end
            P = bsxfun(@plus,s.R_c*P,s.T_c);
        end
        % Visualization functions
        function plot_state(s)
            % Drawing shapes
            [~,J] = s.build_positions();
            quiver3(reshape(J(3,:,:),1,[]),...
                    reshape(J(1,:,:),1,[]),...
                    reshape(J(2,:,:),1,[]),...
                    reshape(J(6,:,:),1,[]),...
                    reshape(J(4,:,:),1,[]),...
                    reshape(J(5,:,:),1,[]),0,'b');
            axis equal;
        end

        % Helper functions
        function E = build_edges(~,skel)
            E = [];
            for i=1:numel(skel.tree)
                for j=1:numel(skel.tree(i).children)
                   E = [E;i skel.tree(i).children(j) ] ;
                end
            end
            [~,idx] = sort(E(:,2));
            E = E(idx,:);
        end
    end
end

