% profile on
rng(1561)
tic
% TODO: Make this a test script by doing if statements
% sim = CircleRiemannianParticleFilterSim(pi/2+pi/4,100);
% sim = SphereRiemannianParticleFilterSim(normc([-1;-1;1]),30);
% sim = Circle(normc([-1;1]),100);
% sim = Sphere(normc([-1;-1;1]),20);
% sim = TangentSphere([normc([-1;-1;1]);0;0.5;0.5],30,true);
% sim = TangentSphereLength([normc([-1;-1;1]);0;0.5;0.5;1.1],20,true);
% sim = TangentSphere([normc([1;1;0]);0;0;0.5],30,true);
% sim.sigma_f = 0.00000000001;
% sim.n_samples = 2000;
sim = TangentSphereProduct(...
        [normc([1;0;0]);0;0;0.5;normc([1;0;0]);0;0;0.5],...
        20);
sim.sigma_f = 0.00000000001;
sim.n_samples = 6000;
results = sim.simulate(1);
toc
% profile off
% profile viewer

figure;
err = acosd(dot(results.simrun.x_gt(1:3,:), results.est(1:3,:)));
plot(err)
hold on
plot(results.Neff/max(results.Neff)*max(err))
hold off


