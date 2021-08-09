function [wGlobal,wGlobalMag]=refFramesAngVeloc(frame1,frame2,frameRate,degOrRad)

% Mitchell Tillman and Jun Liu, Feb. 25 2021

% Purpose:
% Computes angular velocity of a rigid body/moving reference frame rotating from frame 1 to frame 2 in space-fixed coordinates.

% Inputs: two coordinate frames
% Frame 1 and Frame 2: two 3x3 matrices each representing a coordinate frame in global coordinates, where:
%     The rows represent each individual axis
%     The columns are the magnitudes of the X, Y, and Z (or I, J, K) components of each axis.
%     The origin of both coordinate frames are assumed to be at the same point on the rigid body.
%     The axes within each of the reference frames should be orthogonal. 
%     The axes within each of the reference frames do not need to be unit vectors. That is done here.
% frameRate: in Hz. 1/(time elapsed between the two frames). Optional: Default value is 1 if not specified.
% degOrRad: Specifies preference for outputs to be in radians (0) or degrees (1). Optional: Default value is 0 if not specified.

% Outputs:
% wGlobal: The 1x3 angular velocity vector in global (space-fixed) coordinates. Each element is the angular velocity magnitude about that axis (XYZ/IJK).
%     wGlobal also represents the axis of rotation for the rotation between these two frames.
% wGlobalMag: The magnitude of the wGlobal vector (e.g. amount of rotation/second).

%% The following code can be used as inputs to test this function in a simple test case.
% Rotation about only the first axis (row 1). Axes 2 and 3 both rotate 0.7854 rad (45 degrees)

% Example 1 Input:
% frame1=[1 0 0; 0 1 0; 0 0 1];
% frame2=[1 0 0; [0 1 1]/norm([0 1 1]); [0 -1 1]/norm([0 -1 1])]; % Only axes 2 and 3 rotate.
% frameRate=1;

% Example 1 Output:
% wGlobal=[0.7854 0 0]; % In radians
% wGlobalMag=0.7854; % In radians

%% Exclude all data with any NaN in it from processing.
if any(isnan(frame1),'all') || any(isnan(frame2),'all')
    wGlobal=[NaN NaN NaN];
    wGlobalMag=NaN;
    return;
end

if nargin<3 % frameRate not specified.
    frameRate=1;
end
if nargin<4 % degOrRad not specified.
    degOrRad=0;
end

%% Ensure that the reference frames are comprised of unit vectors.
frame1=frame1./sqrt(frame1(:,1).^2+frame1(:,2).^2+frame1(:,3).^2);
frame2=frame2./sqrt(frame2(:,1).^2+frame2(:,2).^2+frame2(:,3).^2);

%% Check that the reference frames are orthogonal.
% NOTE: There is no fix implemented here if they are not orthogonal, so as to not distort the data.
tol=1e-6;
if ~(dot(frame1(1,:),frame1(2,:))<tol && dot(frame1(1,:),frame1(3,:))<tol && dot(frame1(2,:),frame1(3,:))<tol ... % Frame 1
        && dot(frame2(1,:),frame2(2,:))<tol && dot(frame2(1,:),frame2(3,:))<tol && dot(frame2(2,:),frame2(3,:))<tol) % Frame 2
    warning(strcat(['Reference frames not orthogonal within tol ' num2str(tol)]));
    wGlobal=[NaN NaN NaN];
    wGlobalMag=NaN;
    return;
end

%% Linear displacement vector from frame 1 to frame 2 of each axes' endpoints.
d1=[frame2(1,1)-frame1(1,1) frame2(1,2)-frame1(1,2) frame2(1,3)-frame1(1,3)];
d2=[frame2(2,1)-frame1(2,1) frame2(2,2)-frame1(2,2) frame2(2,3)-frame1(2,3)];
d3=[frame2(3,1)-frame1(3,1) frame2(3,2)-frame1(3,2) frame2(3,3)-frame1(3,3)];

%% Obtain the rotation axis, perpendicular to each of the linear distances traveled by axes 1-3.
% Check if any rotation is occurring at all.
if isequal([d1; d2; d3],zeros(size([d1; d2; d3]))) % No rotation is occurring at all!
    wGlobal=[0 0 0];
    wGlobalMag=0;
    return;
end

% Checks if rotation is occurring purely about one axis. If so, uses the other two to define the axis of rotation.
if isequal(d3,zeros(size(d3))) % Case of rotation about purely one of the axes.
    axOfRot=(cross(d1,d2)/norm(cross(d1,d2)));
elseif isequal(d1,zeros(size(d1)))
    axOfRot=(cross(d3,d2)/norm(cross(d3,d2)));
elseif isequal(d2,zeros(size(d2)))
    axOfRot=(cross(d1,d3)/norm(cross(d1,d3)));
else % Rotation is not purely about one axis. Arbitrarily two axes to cross.
    axOfRot=(cross(d1,d3)/norm(cross(d1,d3)));
end

% Isolate axes by frame. 1st row is current frame, 2nd is prev frame.
ax1=[frame2(1,:); frame1(1,:)]; 
ax2=[frame2(2,:); frame1(2,:)];
ax3=[frame2(3,:); frame1(3,:)];

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
    assert((abs(th1-th2)<tol && abs(th3-th2)<tol && abs(th3-th1)<tol) || any(isnan([th1 th2 th3]))); % Rotation axis not aligned with global.
    axUse=ax1RotPlaneProj;
    thUse=th3;
elseif isequal(d1,zeros(size(d1))) % Rotation axis aligned with ax1.
    assert((abs(th2-th3)<tol));
    axUse=ax2RotPlaneProj;
    thUse=th3;
elseif isequal(d2,zeros(size(d2))) % Rotation axis aligned with ax2.
    assert((abs(th1-th3)<tol));
    axUse=ax1RotPlaneProj;
    thUse=th3;
elseif isequal(d3,zeros(size(d3)))% Rotation axis aligned with ax3.
    assert((abs(th2-th1)<tol));
    axUse=ax1RotPlaneProj;
    thUse=th2;
end

% Determine whether rotation is positive or negative (follows the right hand rule!)
if ~isequal(cross(axUse(2,:),axUse(1,:)),zeros(1,3))
    wDir=cross(axUse(2,:),axUse(1,:))/norm(cross(axUse(2,:),axUse(1,:)));
else % 180 degree (pi radians) turn! Axis of rotation is in either one of the two directions.
    wGlobal=[NaN NaN NaN]; % Placeholder!! This is the result of the above computation anyways!
    wGlobalMag=NaN;
    return;
%     if degOrRad==0 % Radians
% %         wDir=;
%     else % Degrees
% %         wDir=;
%     end
end

% Compute the magnitude of the rotation, using the given frame rate.
wGlobalMag=thUse*frameRate;
if degOrRad==1
    wGlobalMag=rad2deg(wGlobalMag);
end

% Multiply the magnitude of the rotation by the direction of the vector.
wGlobal=wGlobalMag*wDir; % Multiply magnitude (deg or rad/s)*direction of cross product. In space-fixed/global coordinates.

