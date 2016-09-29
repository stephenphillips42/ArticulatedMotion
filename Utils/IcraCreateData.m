function [X,x,E] = IcraCreateData(matname,jointsname,AnglePerFrame,Zdes)
% Written by Spiros Leonardos
if nargin <3
    Zdes = 0; % desired depth
end

if nargin < 2
AnglePerFrame = 1;
end

% Load data
ld = load(jointsname);
skel = ld.skel;
mocapIndices = ld.mocapIndices;
ld = load(matname);
S = ld.S;

Sl = S(:,mocapIndices);

%ATTENTION
scale =6;
Sl = Sl/scale;

E = [];
for i=1:numel(skel.tree)
    for j=1:numel(skel.tree(i).children)
       E = [E;i skel.tree(i).children(j) ] ;
    end
end

[~,idx] = sort(E(:,2));
E = E(idx,:);

% Change depth to make it approximately constant
F= size(Sl,1)/3;
P = size(Sl,2);

Angle = (1:F)*(AnglePerFrame/180*pi);

for i=1:F
    
    meanX = mean(Sl(3*i-2,:));
    meanY = mean(Sl(3*i-1,:));
    meanZ = mean(Sl(3*i,:));
    
    
    Sl(3*i-2,:) =  Sl(3*i-2,:) - meanX;
    Sl(3*i-1,:) =  Sl(3*i-1,:) - meanY;
    Sl(3*i,:)   =  Sl(3*i,:)   - meanZ;
    
    % rotate properly here
    
    rotaxis = [0 1 0]';
    Ry = expm(hat(rotaxis)*Angle(i));
    Sl(3*i-2:3*i,:) = Ry*Sl(3*i-2:3*i,:);
    
    
    % translate to the desired depth
      Sl(3*i,:) =  Sl(3*i,:) + Zdes;
end

X = zeros([3 F P]);
x = zeros([2 F P]);
for p=1:P, 
    X(:,:,p) = reshape(Sl(:,p),[3 F]); 
    x(:,:,p) =  [ X(1,:,p)./X(3,:,p);...
                  X(2,:,p)./X(3,:,p)];
end

end
