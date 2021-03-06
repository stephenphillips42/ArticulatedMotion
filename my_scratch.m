% Scratch code
%rng(1561)
rng(5747)
% Begin timing
% tic
nsteps = 4;

% Create graph
% E = [ 1, 2;
%       2, 3;
%       1, 4;
%       4, 5;
%       1, 6 ];
% L = [ 2, 1, 1, 1, 1 ]*0.25;
E = [ 1, 2;
      2, 3;
      1, 4 ];
L = [ 2, 1, 1 ]*0.25;

% Create initial condition
root_pos = [0;0];
thetas = (pi/4)*[1;7;5;3;7];
vels = (pi/64)*[2;-3;3;-3;-1];
% thetas = (pi/4)*[1;7;7];
% vels = (pi/64)*[2;-3;-1];
x0 = zeros(4*size(E,1),1);
for i = 1:size(E,1)
    pinds = (1:2) + 4*(i-1);
    vinds = (3:4) + 4*(i-1);
    x0(pinds) = [cos(thetas(i)); sin(thetas(i))];
    x0(vinds) = vels(i)*[cos(thetas(i)+pi/2); sin(thetas(i)+pi/2)];
end

% profile on
sim = TangentCircleGraph(...
            'x0',x0,...
            'nsteps',nsteps,...
            'lengths',L,...
            'edges',E,...
            'root_pos',root_pos,...
            'nsamples',10000,...
            'pows',[1e-5,1e-4,1e-3,1e-2,1e-1,1e0]);


% sim.simulate(3);
sim.simulate_annealing(3);

% profile off
% profile viewer

% End timing
% toc
