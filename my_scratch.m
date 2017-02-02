% Scratch code
%rng(1561)
rng(5747)
% Begin timing
tic
nsteps = 10;

% Get results
if ~exist('data','var')
    matname = 'Data/15_02.mat';
    jointsname = 'Data/joints15.mat';
    data = LoadFormattedData(matname, jointsname);
end

% profile on
sim = TangentSphereGraph(data.x{1},nsteps,data.L,data.E,[],data.root_pos);
results = sim.simulate(3);
% sim = CircleLogW(normc([-1;1]),500);
% results = sim.simulate(0);
% profile off
% profile viewer

% End timing
toc
