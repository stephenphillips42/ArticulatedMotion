function ParticleSimTest(name, nsteps, plot_error, plot)
% Test various particle filters
if nargin < 3
    plot_error = false;
end
if nargin < 4
    plot = 1;
end

if name == 'Circle'
    sim = Circle(normc([-1;1]),nsteps);
elseif name == 'Sphere'
    sim = Sphere(normc([-1;-1;1]),nsteps);
elseif name == 'TangentSphere'
    sim = TangentSphere([normc([-1;-1;1]);0;0.5;0.5],nsteps,true);
elseif name == 'TangentSphereLength'
    sim = TangentSphereLength([normc([-1;-1;1]);0;0.5;0.5;1.1],nsteps,true);
elseif name == 'TangentSphere'
    sim = TangentSphere([normc([1;1;0]);0;0;0.5],nsteps,true);
    sim.sigma_f = 0.00000000001;
    sim.n_samples = 2000;
elseif name == 'TangentSphereProduct'
    sim = TangentSphereProduct(...
            [normc([1;1;0]);0;0;0.5;normc([1;0;0]);0;0;0.5],...
            nsteps);
    sim.sigma_f = 0.001;
    sim.n_samples = 6000;
else
    error('Particle filter %s not implemented yet',name)
end
% results = sim.simulate(0);
sim.simulate(plot);

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