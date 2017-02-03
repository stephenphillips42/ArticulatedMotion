sim = TangentCircleGraph(...
            'x0',[0;0;0;0],...
            'nsteps',212,...
            'lengths',[1,2],...
            'edges',[1,2]);

sim2 = Circle('x0',normc([-1;1]),'nsteps',nsteps);

ntests = 100;
wsize = [1 8000];
tic
for i = 1:ntests
    w2 = sim.weight_power(randn(wsize)*0.6 - 2,rand);
end
time_log = toc;

tic
for i = 1:ntests
    w3 = sim.normalize((randn(wsize)*0.6 - 2)*rand);
end
time_log2 = toc;

tic
for i = 1:ntests
    w1 = sim2.weight_power(rand(wsize),rand);
end
time_sum = toc;

err = zeros(1,ntests);
w0 = zeros(ntests,wsize(2));
w1 = zeros(ntests,wsize(2));
w2 = zeros(ntests,wsize(2));
w3 = zeros(ntests,wsize(2));
for i = 1:ntests
    ww = randn(wsize)*0.6 - 150;
    p = rand;
    w0(i,:) = sim2.weight_power(exp(ww),p);
    w1(i,:) = (sim2.normalize(sim2.normalize(exp(ww)).^p));
    w2(i,:) = sim.weight_power(ww,p);
    w3(i,:) = sim.normalize(ww*p);
    err(i) = norm(w2(i,:)-w3(i,:));
end

fprintf('Log time: %.03f; Log time 2: %.03f; Sum time: %.03f; Ratio: %e\n',time_log,time_log2,time_sum,time_log/time_sum)

