% profile on
rng(1561)
tic
% sim = CircleRiemannianParticleFilterSim(pi/2+pi/4,100);
% sim = CircleParticleFilterSim(normc([-1;1]),100);
% sim = SphereRiemannianParticleFilterSim(normc([-1;-1;1]),30);
sim = SphereParticleFilterSim(normc([-1;-1;1]),30);
results = sim.simulate(true);
toc
% profile off
% profile viewer

figure;
err = acosd(diag(results.simrun.x_gt.'*results.est));
plot(err)
hold on
plot(results.Neff/max(results.Neff)*max(err))
hold off


