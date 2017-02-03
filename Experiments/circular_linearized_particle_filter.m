%%% Trying to do particle filter on circle

% Space: C = { x | norm(x) == 1 }
% Noise Distribution: von Mises (denoted vM(mu,kappa))
% Motion model: f(x) = R(vM(dtheta,kappa_f))x
% Measurement model: h(x) = R(vM(0,kappa_h))x

% Tools
R = @(theta) [ cos(theta), -sin(theta); sin(theta), cos(theta) ];
angle2x = @(th) [ cos(th); sin(th) ];
x2angle = @(x) atan2(x(2,:),x(1,:));
normalize = @(x) x/norm(x,1);
proj = @(x, d) d - x*(x.'*d);
dist = @(x, y) real(acos(x.'*y));
log_map = @(x, y) normc(proj(x,y-x))*dist(x,y);
% vmlikelihood_angle = @(x,mu,kappa) exp(kappa*cos(x-mu))./(2*pi*besseli(0,kappa));
% vmlikelihood_x = @(x,mu,kappa) exp(kappa*(x.'*mu))./(2*pi*besseli(0,kappa));

% Parameters
dtheta = (1*(pi/180)); % degrees per time unit
kappa_f = 1/(0.1*(pi/180)); % degree variance for movement
kappa_h = 1/(0.2*(pi/180)); % degree variance for measurements
kappa_s = 1/(0.1*(pi/180)); % degree variance for sampling
n = 200;
T = 100;
theta0 = 0;
x0 = angle2x(theta0);

% Build model:
f = @(x) R(dtheta + vmrand(0,kappa_f))*x;
h = @(x) x;

% Create motion and measurements
x_gt = zeros(2,T);
meas = zeros(2,T);
motion_noise = vmrand(0, kappa_f, 1, T-1);
meas_noise = vmrand(0, kappa_h, 1, T);

x_gt(:,1) = x0;
meas(:,1) = R(meas_noise(1))*x_gt(:,1);
for i = 2:T
    x_gt(:,i) = R(dtheta + motion_noise(i-1))*x_gt(:,i-1);
    meas(:,i) = R(meas_noise(i))*x_gt(:,i);
end

% Create samples
sample_thetas = vmrand(theta0, 0.01, 1, n);
samples = [ cos(sample_thetas); sin(sample_thetas) ];
w = normalize(ones(1,length(samples)));
% Visualization
th_viz = linspace(0,2*pi,200);
circ_x_viz = cos(th_viz);
circ_y_viz = sin(th_viz);
Neff = zeros(1,T);
est = zeros(2,T);
l = [0.95 1.05];

for i = 1:T
    figure(1);
    plot(circ_x_viz,circ_y_viz,'g-')
    hold on
    scatter(x_gt(1,i),x_gt(2,i),'k*')
    scatter(meas(1,i),meas(2,i),'rx')
    scatter(samples(1,:).*(1+w),samples(2,:).*(1+w),'b.')
    hold off
    axis equal
    axis([-1.2 1.2 -1.2 1.2])
    pause(0.01)
    % Particle Filter code
    y_meas = meas(:,i);
    for j = 1:n
        samples(:,j) = f(samples(:,j));
        % Linearization method
        x_sample = samples(:,j);
        z = log_map(x_sample, y_meas);
        w(j) = normpdf(norm(z),0,0.1);
    end
    w = normalize(w);
    % samples = samples + normrnd(v,sigma_f,size(samples));
    % w = normalize(normpdf(samples,meas(i),sigma_h));
    % Resampling
    Neff(i) = 1/sum(w.^2);
    if Neff(i) < n/3
        % Taken from Diego Andrés Alvarez Marín's Particle Filter code
        % this is performing latin hypercube sampling on wk
        edges = min([0 cumsum(w)],1); % protect against accumulated round-off
        edges(end) = 1;                 % get the upper edge exact
        U1 = rand/length(w);
        % this works like the inverse of the empirical distribution and returns
        % the interval where the sample is to be found
        [~, idx] = histc(U1:(1/n):1, edges);
        samples = samples(:,idx);
%         for j = 1:n
%             samples(:,j) = R(vmrand(0,kappa_s))*samples(:,j);
%         end
        w = normalize(ones(1,length(samples)));
    end

end


