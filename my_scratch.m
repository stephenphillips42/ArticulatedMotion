% Scratch code
rng(1561)
tic
nsteps = 20;
% TODO: Make this a test script by doing if statements
% sim = Circle(normc([-1;1]),nsteps);
% sim = Sphere(normc([-1;-1;1]),nsteps);
% sim = TangentSphere([normc([-1;-1;1]);0;0.5;0.5],nsteps,true);
% sim = TangentSphereLength([normc([-1;-1;1]);0;0.5;0.5;1.1],nsteps,true);
% sim = TangentSphere([normc([1;1;0]);0;0;0.5],nsteps,true);
% sim.sigma_f = 0.00000000001;
% sim.n_samples = 2000;
% sim = TangentSphereProduct(...
%         [normc([1;1;0]);0;0;0.5;normc([1;0;0]);0;0;0.5],...
%         nsteps);
% sim.sigma_f = 0.001;
% sim.n_samples = 6000;
% results = sim.simulate(0);
% % Plot error
% figure;
% [err, err_full] = sim.compute_error(results.simrun.x_gt, results.est);
% plotyy(1:nsteps,rad2deg(err_full(1,:)),1:nsteps,results.Neff)

% 
% % Draw 3D figure using kinematic chain :D
% T = 200;
% nsamples = 20;
% runtimes = zeros(T-nsamples,1);
% root_pos = PosTrue(:,t,1);
% for t = 1:T-nsamples
%     Xt = cell2mat(x(t:t+nsamples-1).');
%     tic
%     scatter3(0,0,0,'k.')
%     hold on
%     P = zeros(3,num_edges+1,nsamples);
%     J = zeros(6,num_edges,nsamples);
%     P(1,1,:) = root_pos(1);
%     P(2,1,:) = root_pos(2);
%     P(3,1,:) = root_pos(3);
%     for k = 1:size(E)
%         i = E(k,1);
%         j = E(k,2);
%         % Compute position and joint
%         J(1:3,j-1,:) = P(:,i,:);
%         J(4:6,j-1,:) = L(j-1)*Xt((1:3)+6*(j-2),:);
%         P(:,j,:) = P(:,i,:) + J(4:6,j-1,:);
%     end
%     % Drawing shapes
%     quiver3(reshape(J(3,:,:),1,[]),...
%             reshape(J(1,:,:),1,[]),...
%             reshape(J(2,:,:),1,[]),...
%             reshape(J(6,:,:),1,[]),...
%             reshape(J(4,:,:),1,[]),...
%             reshape(J(5,:,:),1,[]),0)
%     scatter3(P(3,:),P(1,:),P(2,:),100,'.')
% 
%     % Drawing measurements
%     plot3([zeros(size(meas{t}(1,:)));ones(size(meas{t}(1,:)))*(Zdes+2)],...
%             [zeros(size(meas{t}(1,:)));meas{t}(1,:)*(Zdes+2)],...
%             [zeros(size(meas{t}(1,:)));meas{t}(2,:)*(Zdes+2)])
%     % scatter3(ones(size(meas{t}(1,:))),meas{t}(1,:),meas{t}(2,:),'k.')
%     hold off
%     axis equal
%     axis([[-0.1, Zdes+1],[-1,1]*1.3,[-1,1]*3 ])
%     xlabel('z')
%     ylabel('x')
%     zlabel('y')
%     runtimes(t) = toc;
%     drawnow()
%     pause(0.01)
% end

% T = 200;
% nsamples = 60;
% runtimes = zeros(T-nsamples,1);
% root_pos = PosTrue(:,t,1);
% for t = 1:T-nsamples
%     Xt = cell2mat(x(t:t+nsamples-1).');
%     tic
%     hold on
%     P = zeros(3,nsamples,num_edges+1);
%     % J = zeros(6,nsamples,num_edges);
%     P(1,:,1) = root_pos(1);
%     P(2,:,1) = root_pos(2);
%     P(3,:,1) = root_pos(3);
%     for k = 1:size(E)
%         i = E(k,1);
%         j = E(k,2);
%         % Compute position (stored for later)
%         P(:,:,j) = P(:,:,i) + L(j-1)*Xt((1:3)+6*(j-2),:);
%         %J(1:3,:,j-1) = P(:,:,i);
%         %J(4:6,:,j-1) = L(j-1)*Xt((1:3)+6*(j-2),:);
%     end
%     runtimes(t) = toc;
% end
% 
% 
% figure
% subplot(1,2,1)
% hist(runtimes)
% subplot(1,2,2)
% hist(1./max(runtimes,2*eps))
% disp(mean(runtimes))
% disp(mean(1./max(runtimes,2*eps)))

%% Get results
if ~exist('x','var')
    matname = 'Data/15_02.mat';
    jointsname = 'Data/joints15.mat';
    [x, meas] = LoadFormattedData(matname, jointsname);
end

% profile on
sim = TangentSphereGraph(x{1},10,L,E,[],root_pos);
results = sim.simulate(3);
% sim = CircleLogW(normc([-1;1]),500);
% results = sim.simulate(0);
% profile off
% profile viewer
