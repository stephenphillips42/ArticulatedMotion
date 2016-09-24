% profile on
rng(1561)
tic
% sim = CircleRiemannianParticleFilterSim(pi/2+pi/4,100);
% sim = SphereRiemannianParticleFilterSim(normc([-1;-1;1]),30);
% sim = Circle(normc([-1;1]),100);
% sim = Sphere(normc([-1;-1;1]),20);
% sim = TangentSphere([normc([-1;-1;1]);0;0.5;0.5],30,true);
sim = TangentSphereLength([normc([-1;-1;1]);0;0.5;0.5;1.1],20,true);
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


