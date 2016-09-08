%%% Particle filter on the plane (super easy)
% Space: C = { x | x in R^1 }
% Noise Distribution: Gaussian (denoted N(mu,sigma^2))
% Motion model: f(x) = x + v + N(0,sigma_f^2)
% Measurement model: h(x) = x + N(0,sigma_h^2)

% Tools
normalize = @(x) x/norm(x,1);

% Parameters
v = 1; % meters per time unit
sigma_f = sqrt(0.5); % degree standard devation for movement
sigma_h = sqrt(2); % degree standard devation for measurements
sigma_s = sigma_h; % degree standard devation for sampling
n = 100;
T = 100;
x0 = 0;

% Create samples
samples = normrnd(0,10,n,1);
w = normalize(ones(size(samples)));

% Create motion and measurements
x_gt = zeros(T,1);
x_gt(1) = x0;
motion_noise = normrnd(0, sigma_f, T-1, 1);
for i = 2:T
    x_gt(i) = x_gt(i-1) + v + motion_noise(i-1);
end
meas_noise = normrnd(0, sigma_h, T, 1);
meas = x_gt + meas_noise;

% Visualization of filter
x_viz = linspace(-2,100,200);
Neff = zeros(T,1);
for i = 1:T
    % Drawing code
    figure(1);
    subplot(1,2,1)
    plot(x_viz,zeros(size(x_viz)),'g-')
    hold on
    scatter(x_gt(i),0,'k*')
    scatter(meas(i),0,'rx')
    scatter(samples,w,'b.')
    hold off
%     axis equal
    pause(0.05)
    subplot(1,2,2)
    plot(Neff)
    % Particle Filter code
    samples = samples + normrnd(v,sigma_f,size(samples));
    w = normalize(normpdf(samples,meas(i),sigma_h));
    % Resampling
    Neff(i) = 1/sum(w.^2);
    if Neff(i) < n/3
        % Taken from Diego Andrés Alvarez Marín's Particle Filter code
        % this is performing latin hypercube sampling on wk
        edges = min([0 cumsum(w)'],1); % protect against accumulated round-off
        edges(end) = 1;                 % get the upper edge exact
        U1 = rand/length(w);
        % this works like the inverse of the empirical distribution and returns
        % the interval where the sample is to be found
        [~, idx] = histc(U1:(1/n):1, edges);
        samples = samples(idx) + normrnd(0,sigma_s,size(samples));
        w = normalize(ones(size(samples)));
    end
end
