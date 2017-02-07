function ParticleSimTest(name, nsteps, plot_error, plot, anneal)
% Test various particle filters
if nargin < 3
    plot_error = false;
end
if nargin < 4
    plot = 1;
end
if nargin < 5
    anneal = false;
end

if strcmp(name, 'Circle')
    sim = Circle('x0',normc([-1;1]),'nsteps',nsteps);
elseif strcmp(name, 'Sphere')
    sim = Sphere(normc([-1;-1;1]),nsteps);
elseif strcmp(name, 'TangentSphere')
    sim = TangentSphere('x0',[normc([-1;-1;1]);0;0.5;0.5],'nsteps',nsteps,'plot_type',true);
elseif strcmp(name, 'TangentSphereLength')
    sim = TangentSphereLength('x0',[normc([-1;-1;1]);0;0.5;0.5;1.1],'nsteps',nsteps,'plot_type',true);
elseif strcmp(name, 'TangentSphere2')
    sim = TangentSphere('x0',[normc([1;1;0]);0;0;0.5],'nsteps',nsteps,'plot_type',true);
    sim.sigma_f = 0.00000000001;
    sim.n_samples = 2000;
elseif strcmp(name, 'TangentSphereProduct')
    sim = TangentSphereProduct(...
            'x0',[normc([1;1;0]);0;0;0.5;normc([1;0;0]);0;0;0.5],...
            'nsteps',nsteps);
    sim.sigma_f = 0.001;
    sim.n_samples = 6000;
elseif strcmp(name, 'TangentSphereGraph')
    % Get results
    tic
    matname = 'Data/15_02.mat';
    jointsname = 'Data/joints15.mat';
    data = LoadFormattedData(matname, jointsname, nsteps);
    toc
    th = pi/4;
    R_h = [ cos(th), 0, sin(th);
                  0, 1,       0;
           -sin(th), 0, cos(th) ];
    sim = TangentSphereGraph(...
                'x0',data.x{1},...
                'nsteps',nsteps,...
                'lengths',data.L,...
                'edges',data.E,...
                'root_pos',data.root_pos,...
                'T_h',[0;0;0],...
                'R_h',R_h);
else
    error('Particle filter %s not implemented yet',name)
end
% results = sim.simulate(0);
tic
if anneal
    sim.simulate_annealing(plot);
else
    sim.simulate(plot);
end
toc

% Plot error
if plot_error
    error('Plot errors not implemented yet')
    % TODO: Make ParticleFilterSim have a plot_error function not just compute
    % error
    % figure;
    % [err, err_full] = sim.compute_error(results.simrun.x_gt, results.est);
    % plotyy(1:nsteps,rad2deg(err_full(1,:)),1:nsteps,results.Neff)
end

end