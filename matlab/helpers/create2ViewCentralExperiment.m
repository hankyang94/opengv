function [v1, v2, t, R] = create2ViewCentralExperiment( noise, pt_number, ...
    fov, max_parallax, outlier_fraction, debug )
%   Author: Heng Yang hankyang@mit.edu
%   Last update: 11/18/2018
%   Adapted from opengv matlab function create2D2DExperimen.m
%   Sample only 3D points that are within FOV of the first camera
%   If debug = true, will plot camera frames, sample points and bearing
%   vectors

%% generate the camera system
% camera simply coincides with the view point
cam_number = 1;
cam_offsets = zeros(3,1);

%% generate random view-points
max_rotation = 0.5;

position1 = zeros(3,1);
rotation1 = eye(3);

position2 = max_parallax * 2.0 * (rand(3,1) - repmat(0.5,3,1));
rotation2 = generateBoundedR(max_rotation);

%% Generate random point-cloud within FOV
minDepth = 1.0;
maxDepth = 5.0;
deg2rad = pi/180.0;

nSamples = 0;
points = zeros(3, pt_number);
points_c2 = zeros(3, pt_number);
pose_c2_c1 = [rotation2, position2; zeros(1,3), 1];
while nSamples < pt_number
    % sample inside a bounded cone within FOV of the first camera
    alpha = -pi + 2*pi*rand(1); % rotation about z
    theta = (-fov/2 + fov*rand(1))*deg2rad; % deviation from z
    direction = [sin(theta)*cos(alpha); sin(theta)*sin(alpha); cos(theta)];
    depth = minDepth + (maxDepth - minDepth) * rand(1);
    point = direction * (depth / direction(end)); % make point's z = depth
    % convert point from first camera to second camera
    point_c2_h = pose_c2_c1 \ [point; 1];
    point_c2 = point_c2_h(1:3);
    % check if point in the cone of second camera 
    theta_2 = atan2(sqrt(point_c2(1)^2 + point_c2(2)^2), point_c2(3));
    if (-fov/2*deg2rad < theta_2 && theta_2 < fov/2*deg2rad)
        nSamples = nSamples + 1;
        points(:, nSamples) = point;
        points_c2(:, nSamples) = point_c2;
    end
end

% normalizedPoints = 2.0*(rand(3,pt_number)-repmat(0.5,3,pt_number));
% norms = sqrt(sum(normalizedPoints.*normalizedPoints));
% directions = normalizedPoints./repmat(norms,3,1);
% points = (maxDepth-minDepth) * normalizedPoints + minDepth * directions;

if debug
%     disp('points')
%     disp(points)
%     disp('points in c2')
%     disp(points_c2)
    figure
    scatter3(points(1,:), points(2,:), points(3,:), 'red', 'filled')
    hold on
    scatter3(position1(1), position1(2), position1(3), 'black', '*')
    quiver3(position1(1), position1(2), position1(3), ...
    0.5, 0, 0, 0, 'green', 'LineWidth', 3)
    quiver3(position1(1), position1(2), position1(3), ...
    0, 0.5, 0, 0, 'blue', 'LineWidth', 3)
    quiver3(position1(1), position1(2), position1(3), ...
    0, 0, 0.5, 0, 'red', 'LineWidth', 3)
    text(position1(1), position1(2), position1(3), 'Camera 1')
    
    c2_coords = [0.5,0,0;0,0.5,0;0,0,0.5];
    c2_coords_in_c1 = pose_c2_c1 * [c2_coords;1,1,1];
%     disp(c2_coords)
%     disp(c2_coords_in_c1)
%     disp(position2)
    scatter3(position2(1), position2(2), position2(3), 'black', '*')
    quiver3(position2(1), position2(2), position2(3), ...
    c2_coords_in_c1(1,1)-position2(1), ...
    c2_coords_in_c1(2,1)-position2(2), ...
    c2_coords_in_c1(3,1)-position2(3), ...
    0, 'green', 'LineWidth', 3)
    quiver3(position2(1), position2(2), position2(3), ...
    c2_coords_in_c1(1,2)-position2(1), ...
    c2_coords_in_c1(2,2)-position2(2), ...
    c2_coords_in_c1(3,2)-position2(3), ...
    0, 'blue', 'LineWidth', 3)
    quiver3(position2(1), position2(2), position2(3), ...
    c2_coords_in_c1(1,3)-position2(1), ...
    c2_coords_in_c1(2,3)-position2(2), ...
    c2_coords_in_c1(3,3)-position2(3),...
    0, 'red', 'lineWidth', 3)
    text(position2(1), position2(2), position2(3), 'Camera 2')
    axis equal
end

%% Now create the correspondences by looping through the cameras

focal_length = 800.0;

v1 = zeros(6,pt_number);
v2 = zeros(6,pt_number);
cam_correspondence = 1;

for i=1:pt_number
    
    cam_offset = cam_offsets(:,cam_correspondence); % zero vector
    
    body_point1 = rotation1' * (points(:,i)-position1);
    body_point2 = rotation2' * (points(:,i)-position2);
    
    % we actually omit the can rotation here by unrotating the bearing
    % vectors already
    bearingVector1 = body_point1 - cam_offset;
    bearingVector2 = body_point2 - cam_offset;
    bearingVector1_norm = norm(bearingVector1);
    bearingVector2_norm = norm(bearingVector2);
    bearingVector1 = bearingVector1/bearingVector1_norm;
    bearingVector2 = bearingVector2/bearingVector2_norm;
    
    % add noise to the bearing vectors here
    bearingVector1_noisy = addNoise(bearingVector1,focal_length,noise);
    bearingVector2_noisy = addNoise(bearingVector2,focal_length,noise);
    
    % store the normalized bearing vectors along with the cameras they are
    % being seen (we create correspondences that always originate from the
    % same camera, you can change this if you want)
    bearingVector1_norm = norm(bearingVector1_noisy);
    bearingVector2_norm = norm(bearingVector2_noisy);
    
    v1(:,i) = [bearingVector1_noisy./bearingVector1_norm; cam_offset];
    v2(:,i) = [bearingVector2_noisy./bearingVector2_norm; cam_offset];
end

if debug
    scale = 5;
    for i=1:pt_number
        v2_rot = rotation2 * v2(1:3,i);
        quiver3(position1(1), position1(2), position1(3),...
            scale*v1(1,i), scale*v1(2,i), scale*v1(3,i), 0, 'yellow')
        quiver3(position2(1), position2(2), position2(3),...
            scale*v2_rot(1), scale*v2_rot(2), scale*v2_rot(3), 0, 'yellow')
    end
end

%% Add outliers
number_outliers = floor(outlier_fraction*pt_number);

if number_outliers > 0
    for i=1:number_outliers

        cam_correspondence = cam_correspondences(1,i);

        cam_offset = cam_offsets(:,cam_correspondence);
        %cam_rotation = cam_rotations(:,(cam_correspondence-1)*3+1:(cam_correspondence-1)*3+3);

        %generate random point
        normalizedPoint = 2.0*(rand(3,1)-repmat(0.5,3,1));
        norm1 = sqrt(sum(normalizedPoint.*normalizedPoint));
        direction = normalizedPoint./norm1;
        point = (maxDepth-minDepth) * normalizedPoint + minDepth * direction;

        body_point2 = rotation2' * (point-position2);

        % store the point (no need to add noise)
        bearingVector2 = body_point2 - cam_offset;

        % store the normalized bearing vectors along with the cameras they are
        % being seen (we create correspondences that always originate from the
        % same camera, you can change this if you want)
        bearingVector2_norm = norm(bearingVector2);

        v2(:,i) = [bearingVector2./bearingVector2_norm; cam_offset];
    end
end

%% compute relative translation and rotation
R = rotation1' * rotation2;
t = rotation1' * (position2 - position1);

%% Cut the cam offset in the single camera case (e.g. central)
if cam_number == 1
    v1 = v1(1:3,:);
    v2 = v2(1:3,:);
end