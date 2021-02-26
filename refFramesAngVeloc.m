function [wGlobal,wGlobalMag]=refFramesAngVeloc(frame1,frame2,frameRate,degOrRad)

% Inputs: two coordinate frames
% Frame 1 and Frame 2 are two 3x3 each representing a coordinate frame in global coordinates, where:
%     The rows represent each individual axis
%     The columns are the X, Y, and Z coordinates of the endpoint of each axis.
%     The origin of both coordinate frames is assumed to be at [0 0 0].
% frameRate: 1/time between the two frames.
% degOrRad: Specifies preference of outputs to be in degrees (0) or radians (1).

% Outputs:
% wGlobal: The 1x3 angular velocity vector in global coordinates
% wGlobalMag: The magnitude of the wGlobal vector (e.g. amount of rotation/second).

% Function:
% Computes angular velocity from frame 1 to frame 2.

% Linear distance from frame 1 to frame 2 between each axes' endpoints.
d1=[frame2(1,1)-frame1(1,1) frame2(1,2)-frame1(1,2) frame2(1,3)-frame1(1,3)];
d2=[frame2(2,1)-frame1(2,1) frame2(2,2)-frame1(2,2) frame2(2,3)-frame1(2,3)];
d3=[frame2(3,1)-frame1(3,1) frame2(3,2)-frame1(3,2) frame2(3,3)-frame1(3,3)];

% Checks if rotation is occurring purely about one axis. If so, uses the other two to define the axis of rotation.
if isequal(d3,zeros(size(d3))) % Case of rotation about purely one of the axes.
    axOfRot=(cross(d1,d2)/norm(cross(d1,d2)));
elseif isequal(d1,zeros(size(d1)))
    axOfRot=(cross(d3,d2)/norm(cross(d3,d2)));
elseif isequal(d2,zeros(size(d2)))
    axOfRot=(cross(d1,d3)/norm(cross(d1,d3)));
end

% Isolate axes by frame. 1st row is current frame, 2nd is prev frame.
ax1=[frame2(1,:)'; frame1(1,:)']; 
ax2=[frame2(2,:)'; frame1(2,:)'];
ax3=[frame2(3,:)'; frame1(3,:)'];

% Project the (1) current, and (2) previous axes onto the plane perpendicular to the axis of rotation.
ax1RotPlaneProj=[ax1(1,:)-dot(ax1(1,:),axOfRot)*axOfRot; ax1(2,:)-dot(ax1(2,:),axOfRot)*axOfRot]; % Axis 1
ax2RotPlaneProj=[ax2(1,:)-dot(ax2(1,:),axOfRot)*axOfRot; ax2(2,:)-dot(ax2(2,:),axOfRot)*axOfRot]; % Axis 2
ax3RotPlaneProj=[ax3(1,:)-dot(ax3(1,:),axOfRot)*axOfRot; ax3(2,:)-dot(ax3(2,:),axOfRot)*axOfRot]; % Axis 3

% Find rotation (theta) of each axes in the plane of rotation. These should be the same, except in the case of uniaxial rotation!
th1 = atan2(norm(cross(ax1RotPlaneProj(1,:),ax1RotPlaneProj(2,:))),dot(ax1RotPlaneProj(1,:),ax1RotPlaneProj(2,:)));
th2 = atan2(norm(cross(ax2RotPlaneProj(1,:),ax2RotPlaneProj(2,:))),dot(ax2RotPlaneProj(1,:),ax2RotPlaneProj(2,:))); 
th3 = atan2(norm(cross(ax3RotPlaneProj(1,:),ax3RotPlaneProj(2,:))),dot(ax3RotPlaneProj(1,:),ax3RotPlaneProj(2,:))); 

% Check the theta values
if ~(isequal(d3,zeros(size(d3))) || isequal(d1,zeros(size(d1))) || isequal(d2,zeros(size(d2)))) % Ensure not rotation about purely one axis.
    assert((abs(th1-th2)<0.0001 && abs(th3-th2)<0.0001 && abs(th3-th1)<0.0001) || any(isnan([th1 th2 th3]))); % Rotation axis not aligned with global.
    thUse=th3;
elseif isequal(d1,zeros(size(d1))) % Rotation axis aligned with X.
    assert((abs(th2-th3)<0.0001));
    vUse=v2;
    thUse=th3;
elseif isequal(d2,zeros(size(d2))) % Rotation axis aligned with Y.
    assert((abs(th1-th3)<0.0001));
    vUse=v1;
    thUse=th3;
elseif isequal(d3,zeros(size(d3)))% Rotation axis aligned with Z.
    assert((abs(th2-th1)<0.0001));
    vUse=v1;
    thUse=th2;
end

% Determine whether rotation is positive or negative (follows the right hand rule!)
wDir=cross(ax1(2,:),ax1(1,:))/norm(cross(ax1(2,:),ax1(1,:)));

% Determine the magnitude of the rotation.
wGlobalMag=thUse*frameRate;
wGlobal=wGlobalMag*wDir; % Multiply magnitude (deg/s)*direction of cross product. IN GLOBAL.

