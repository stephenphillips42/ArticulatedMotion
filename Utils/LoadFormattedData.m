function [ x, meas ] = LoadFormattedData(matname, jointsname)
%LOADFORMATTEDDATA Puts ICRA data in format that we want

Zdes = 3;
AnglePerFrame = 0.1;
[PosTrue,Meas,E] = IcraCreateData(matname,jointsname,AnglePerFrame,Zdes);
% Get info for initial configuration of the human
edges = E;
num_edges = size(E,1);
P1 = squeeze(PosTrue(:,1,:));
P2 = squeeze(PosTrue(:,2,:));
x = cell(nsteps,1);
meas = cell(nsteps,1);
x0 = zeros(6*num_edges,1);
root_pos = PosTrue(:,1,1);
% Extract lengths
Ls = zeros(8000,num_edges);
for j = 1:8000;
    Xj = squeeze(PosTrue(:,j,:));
    for i = 1:length(E);
        Ls(j,i) = norm(Xj(:,E(i,1))-Xj(:,E(i,2)));
    end
end
stdL = std(Ls);
L = mean(Ls);
% Extract initial position and velocity (roughly)
pinds = @(i) (1:3) + 6*(i-1); % Position indecies for x0
vinds = @(i) (4:6) + 6*(i-1); % Velocity indecies for x0
for t = 1:nsteps
    x{t} = zeros(6*num_edges,1);
    Pm = squeeze(PosTrue(:,t,:));
    Pm1 = squeeze(PosTrue(:,min(nsteps,t+1),:));
    for k = 1:num_edges
        i = E(k,1);
        j = E(k,2);
        x{t}(pinds(j-1)) = normc(Pm(:,j) - Pm(:,i));
        x1 = normc(Pm1(:,j) - Pm1(:,i));
        % Log map of the sphere 
        % TODO: functionize this later
        trxy = x{t}(pinds(j-1)).'*x1;
        geodist = acos(trxy);
        x{t}(vinds(j-1)) = geodist*(x1-x{t}(pinds(j-1))*trxy)/sqrt(1-trxy.^2);
    end
    meas{t} = squeeze(Meas(:,t,:));
end

end

