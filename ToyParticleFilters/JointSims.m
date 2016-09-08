% Simulations

Rx = @(th) [ 1, 0, 0; 0, cos(th), -sin(th); 0, sin(th), cos(th)];
Ry = @(th) [ cos(th), 0, sin(th); 0, 1, 0; -sin(th), 0, cos(th)];
Rz = @(th) [ cos(th), -sin(th), 0; sin(th), cos(th), 0; 0, 0, 1];

% Create graph 
% TODO: Make a class??
% Format: The graph g is a cell array
% g{i} gives the adjacency list of node i, which is a k x 2 array
% The first column show what edge it is connected to. The second column shows
% the length of the edge. This unfortunately means there is not consistency
% Enforced between the two
g = {
    [2, 3],
    [1, 3; 3, 1],
    [2, 1]
};

% Instance, must follow the form set by the graph g
% TODO: Make a class as well??
sph = @(u,v) [ sin(u)*cos(v), sin(u)*sin(v), cos(u) ];
points_ = @(u1,v1,u2,v2) [                               zeros(1,3);
                                               g{1}(1,2)*sph(u1,v1);
                        g{1}(1,2)*sph(u1,v1) + g{2}(2,2)*sph(u2,v2) ];

nspacings = 50;
u1 = linspace(0,pi,nspacings);
v1 = linspace(0,pi,nspacings);
u2 = linspace(0,pi/2,nspacings);
v2 = linspace(0,pi/2,nspacings);

% R = Rx(-pi/4)*Ry(-pi/4)*Rx(-pi/2);
R = Rx(-pi/4)*Ry(-pi/4)*Rx(pi/2)*Rz(pi);
T = [0,0,-4];

testpoints = [ 0  0  0;
               0  0  1;
               0  0 -1;
               0  0  0;
               0  1  0;
               0 -1  0;
               0  0  0;
               1  0  0;
              -1  0  0 ];
colors = [  0, 0, 0;
            0, 1, 1;
            1, 1, 1;
            0, 0, 0;
            0, 1, 0;
            1, 1, 1;
            0, 0, 0;
            1, 0, 0;
            1, 1, 1 ];
testcampoints = testpoints*R.' + ones(length(testpoints),1)*T;
testcampoints = bsxfun(@rdivide, testcampoints, testcampoints(:,3));
% hold on
% plot(testcampoints(:,1),testcampoints(:,2),'b-','LineWidth',1);
% scatter(testcampoints(:,1),testcampoints(:,2),100,colors,'.');
% hold off
% axis equal;
% axis([-3,3,-3,3]);
% xlabel('x')
% ylabel('y')


points = [];
campoints = [];
for i = 1:50
    curpoints = points_(u1(i),v1(i),u2(i),v2(i));
    points = [points; curpoints];
    figure(1)
    subplot(1,2,1)
    scatter3(testpoints(:,1),testpoints(:,2),testpoints(:,3),150,colors,'.');
    hold on;
    plot3(testpoints(:,1),testpoints(:,2),testpoints(:,3),'b-');
    plot3(curpoints(:,1),curpoints(:,2),curpoints(:,3),'b-','LineWidth',1);
    scatter3(points(:,1),points(:,2),points(:,3),150,'b.');
    hold off;
    axis equal;
    axis([-4 4 -4 4 -4 4]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    % hold on
    subplot(1,2,2);
    curcampoints = curpoints*R.' + ones(length(curpoints),1)*T;
    curcampoints = bsxfun(@rdivide, curcampoints, curcampoints(:,3));
    campoints = [campoints; curcampoints];
    scatter(testcampoints(:,1),testcampoints(:,2),100,colors,'.');
    hold on;
    plot(testcampoints(:,1),testcampoints(:,2),'b-','LineWidth',1);
    plot(curcampoints(:,1),curcampoints(:,2),'b-','LineWidth',1);
    scatter(campoints(:,1),campoints(:,2),100,'b.');
    hold off;
    axis equal;
    axis([-3,3,-3,3]);
    xlabel('x')
    ylabel('y')

    pause(0.08)
end




