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
sim = TangentSphereProduct(...
        [normc([1;1;0]);0;0;0.5;normc([1;0;0]);0;0;0.5],...
        nsteps);
sim.sigma_f = 0.001;
sim.n_samples = 6000;
results = sim.simulate(0);
%% Plot error
figure;
[err, err_full] = sim.compute_error(results.simrun.x_gt, results.est);
% yyaxis left
% plot(err)
% yyaxis right
% plot(results.Neff/max(results.Neff)*max(err))
% plotyy(1:nsteps,err_full(1:2,:),1:nsteps,results.Neff)
plotyy(1:nsteps,rad2deg(err_full(1:2,:)),1:nsteps,results.Neff)

% Create data TODO: Make this a function wrapping IcraCreateData
% tic
% Zdes = 3;
% AnglePerFrame = 0.1;
% [PosTrue,Meas,E] = IcraCreateData('Data/15_02.mat','Data/joints15.mat',AnglePerFrame,Zdes);
% % Get info for initial configuration of the human
% edges = E;
% num_edges = size(E,1);
% P1 = squeeze(PosTrue(:,1,:));
% P2 = squeeze(PosTrue(:,2,:));
% x = cell(nsteps,1);
% meas = cell(nsteps,1);
% x0 = zeros(6*num_edges,1);
% root_pos = PosTrue(:,t,1);
% % Extract lengths
% Ls = zeros(8000,num_edges);
% for j = 1:8000;
%     Xj = squeeze(PosTrue(:,j,:));
%     for i = 1:length(E);
%         Ls(j,i) = norm(Xj(:,E(i,1))-Xj(:,E(i,2)));
%     end
% end
% stdL = std(Ls);
% L = mean(Ls);
% % Extract initial position and velocity (roughly)
% pinds = @(i) (1:3) + 6*(i-1); % Position indecies for x0
% vinds = @(i) (4:6) + 6*(i-1); % Velocity indecies for x0
% for t = 1:nsteps
%     x{t} = zeros(6*num_edges,1);
%     Pm = squeeze(PosTrue(:,t,:));
%     Pm1 = squeeze(PosTrue(:,min(nsteps,t+1),:));
%     for k = 1:num_edges
%         i = E(k,1);
%         j = E(k,2);
%         x{t}(pinds(j-1)) = normc(Pm(:,j) - Pm(:,i));
%         x1 = normc(Pm1(:,j) - Pm1(:,i));
%         % Log map of the sphere 
%         % TODO: functionize this later
%         trxy = x{t}(pinds(j-1)).'*x1;
%         geodist = acos(trxy);
%         x{t}(vinds(j-1)) = geodist*(x1-x{t}(pinds(j-1))*trxy)/sqrt(1-trxy.^2);
%     end
%     meas{t} = squeeze(Meas(:,t,:));
% end
% toc

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

% Get results
% profile on
% sim = TangentSphereGraph(x{1},nsteps,L,E,[],root_pos);
% results = sim.simulate(0);
% profile off
% profile viewer
% toc
% 
% % TODO Make the error printing... better, more modular
% figure;
% err = acosd(dot(results.simrun.x_gt(1:ndims,:), results.est(1:ndims,:)));
% plot(err)
% hold on
% plot(results.Neff/max(results.Neff)*max(err))
% hold off


