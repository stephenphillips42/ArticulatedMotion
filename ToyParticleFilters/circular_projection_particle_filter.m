%%% Trying to do particle filter on circle

% Space: C = { x in R^3 | norm(x) == 1 }
% Noise Distribution: Gaussian on the tangent space
% Motion model: f(x) = R(vM(dtheta,kappa_f))x
% Measurement model: h(x) = R(vM(0,kappa_h))x

% Tools
R = @(theta) [ cos(theta), -sin(theta); sin(theta), cos(theta) ];
angle2x = @(th) [ cos(th); sin(th) ];
x2angle = @(x) atan2(x(2,:),x(1,:));
normalize = @(x) x/norm(x,1);
proj = @(x,y) @(x, d) d - x*(x.'*d);
projlikelihood = @(x,z,offset,sigma) normpdf(x(1)/(x(2)+offset)-z,0,sigma);

circmean = @(X) normc(sum(X,2));

% Parameters
dtheta = (1*(pi/180)); % degrees per time unit
kappa_f = 1/(0.1*(pi/180)); % degree variance for movement
sigma_h = sqrt(0.1*(pi/180)); % degree variance for measurements
projoffset = 3;
n = 300;
T = 100;
theta0 = pi/2;
x0 = angle2x(theta0);

% Create motion and measurements
x_gt = zeros(2,T);
meas = zeros(1,T);
motion_noise = vmrand(0, kappa_f, 1, T-1);
meas_noise = normrnd(0, sigma_h, 1, T);

x_gt(:,1) = x0;
meas(1) = x_gt(1,1)/(x_gt(2,1) + projoffset);
for i = 2:T
    x_gt(:,i) = R(dtheta + motion_noise(i-1))*x_gt(:,i-1);
    meas(i) = x_gt(1,i)/(x_gt(2,i) + projoffset) + normrnd(0,sigma_h);
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
    subplot(1,2,1)
    plot(circ_x_viz,circ_y_viz,'g-')
    hold on
    plot(x_gt(1,i)*l,x_gt(2,i)*l,'k-')
    plot([0,meas(i)*10],[-projoffset,(-projoffset+10)],'r-')
    est(:,i) = circmean(samples);
    plot(est(1,i)*l,est(2,i)*l,'c-')
    scatter(samples(1,:).*(1+w),samples(2,:).*(1+w),'b.')
    hold off
    axis equal
    axis([-1.2 1.2 -1.2 1.2])
    subplot(1,2,2)
    plot(w)
    pause(0.01)
    % Particle Filter code
    for j = 1:n
        samples(:,j) = R(dtheta + vmrand(0,kappa_f))*samples(:,j);
        w(j) = projlikelihood(samples(:,j),meas(i),projoffset,sigma_h);
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
        w = normalize(ones(1,length(samples)));
    end
end


