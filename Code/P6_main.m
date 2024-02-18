clc
clear
close all

load('Variables_6.mat')


% Let's shift the Center 1 to 0,0,0
Cylinder_Surface_World_Shift_Rotate(1,:)=Cylinder_Surface_World_True(1,:)-Center_C1_W_True(1);
Cylinder_Surface_World_Shift_Rotate(2,:)=Cylinder_Surface_World_True(2,:)-Center_C1_W_True(2);
Cylinder_Surface_World_Shift_Rotate(3,:)=Cylinder_Surface_World_True(3,:)-Center_C1_W_True(3);

% Let's rotate everything so that the axis is aligned with the Z axis
axis_direction_vector_True = axis_direction_vector_True / norm(axis_direction_vector_True);
rotation_axis = cross(axis_direction_vector_True, [0, 0, 1]);      % Calculate the rotation axis to Z axis        
rotation_axis = rotation_axis / norm(rotation_axis);

rotation_angle = acos(dot(axis_direction_vector_True, [0, 0, 1])); % Calculate the rotation angle
rotation_matrix = vrrotvec2mat([rotation_axis, rotation_angle]);   % Create the rotation matrix

% Rotate the points
Cylinder_Surface_World_Shift_Rotate(1:3,:) = rotation_matrix * Cylinder_Surface_World_Shift_Rotate(1:3,:);
rotated_axis = rotation_matrix * axis_direction_vector_True;

% Original and Rototranslated Surface in 3D Space--------------------------
figure(1);
xlabel('X');ylabel("Y");zlabel("Z")
scatter3(0,0,0,'r')
text(0,0,0,'camera')
hold on;
axis equal;

% Plot the Original Surface
scatter3(Cylinder_Surface_World_True(1,:),Cylinder_Surface_World_True(2,:),Cylinder_Surface_World_True(3,:),36,[0.6350 0.0780 0.1840]);
hold on;

% Plot the rotated Surface
scatter3(Cylinder_Surface_World_Shift_Rotate(1, :), Cylinder_Surface_World_Shift_Rotate(2, :), Cylinder_Surface_World_Shift_Rotate(3, :),36,[0.4660 0.6740 0.1880]);
hold on;
xlabel('X');
ylabel('Y');
zlabel('Z');



%----------------------------------UNWRAP---------------------------------

% Unwrap by getting everything in cylindrical coordinates but plotting on
% to a plane.

Cylinder_Surface_h_tetha(1,:)=Cylinder_Surface_World_Shift_Rotate(3,:);
Cylinder_Surface_h_tetha(2,:)=atand(Cylinder_Surface_World_Shift_Rotate(2,:)./Cylinder_Surface_World_Shift_Rotate(1,:));


h_min=min(Cylinder_Surface_h_tetha(1,:))
h_max=max(Cylinder_Surface_h_tetha(1,:))
tetha_min=min(Cylinder_Surface_h_tetha(2,:))
tetha_max=max(Cylinder_Surface_h_tetha(2,:))


% There are points from tetha -90 to 90 degrees. Nearly a half circle.
% This half circle has length pi*Radius.
% The height is Distance.
% So radially there are 1.6254 points for every 1 point on the height difference.
% Take this as our aspect ratio of the unwrapped image.
% Coordinates defined by height should have length 1
% Coordinates defined by the angle should have length 1.6
% The final image will be about 4457515 pixels.

img_aspectratio=deg2rad(abs(tetha_max-tetha_min))*Radius/Distance;
unwrap_size_y=round(sqrt(length(Cylinder_Surface_h_tetha)/img_aspectratio));
unwrap_size_x=round(img_aspectratio*unwrap_size_y);

% Map the h and tetha values to 0->hmax-hmin, 0->tethamax-tethamin
Cylinder_Surface_h_tetha(1,:)=Cylinder_Surface_h_tetha(1,:)-h_min;
Cylinder_Surface_h_tetha(2,:)=Cylinder_Surface_h_tetha(2,:)-tetha_min;

% Map the h and tetha values to 1->size_y+1, 1->size_x+1 (everything
% starting from 1,1 is easier for indexing later on.
Cylinder_Surface_h_tetha(1,:)=round(Cylinder_Surface_h_tetha(1,:).*unwrap_size_y./(h_max-h_min))+1;
Cylinder_Surface_h_tetha(2,:)=round(Cylinder_Surface_h_tetha(2,:).*unwrap_size_x./(tetha_max-tetha_min))+1;


IMG_unwrapped=zeros(unwrap_size_y+1, unwrap_size_x+1, 3, 'uint8');
for i=1:length(Cylinder_Surface_h_tetha)
    IMG_unwrapped(Cylinder_Surface_h_tetha(1,i),Cylinder_Surface_h_tetha(2,i),1)=Cylinder_Surface_World_True(4,i);  %R
    IMG_unwrapped(Cylinder_Surface_h_tetha(1,i),Cylinder_Surface_h_tetha(2,i),2)=Cylinder_Surface_World_True(5,i);  %G
    IMG_unwrapped(Cylinder_Surface_h_tetha(1,i),Cylinder_Surface_h_tetha(2,i),3)=Cylinder_Surface_World_True(6,i);  %B
end

figure(2)
imshow(imrotate(circshift(IMG_unwrapped, [0, -unwrap_size_x/2-100]), 180, 'bilinear', 'crop'))
% Manual beautification so that the image isn't cut in half and it's right way up



%Shift everything down to square the pixels.
for i=1:unwrap_size_x
    % Find columns where x value is the same.
    columns_with_x_equals_i = find(Cylinder_Surface_h_tetha(2, :) == i);

    % Extract y values corresponding to the same x value.
    y_values_for_same_x = Cylinder_Surface_h_tetha(1, columns_with_x_equals_i);

    % Find the smallest y=h value
    min_y_value_for_this_x = min(y_values_for_same_x);

    %Shift every pixel down so that the lowest pixels in one row is at h=0
    Cylinder_Surface_h_tetha(1, columns_with_x_equals_i) = Cylinder_Surface_h_tetha(1, columns_with_x_equals_i) - min_y_value_for_this_x +1;
end


IMG_unwrapped=zeros(unwrap_size_y+1, unwrap_size_x+1, 3, 'uint8');
for i=1:length(Cylinder_Surface_h_tetha)
    IMG_unwrapped(Cylinder_Surface_h_tetha(1,i),Cylinder_Surface_h_tetha(2,i),1)=Cylinder_Surface_World_True(4,i);  %R
    IMG_unwrapped(Cylinder_Surface_h_tetha(1,i),Cylinder_Surface_h_tetha(2,i),2)=Cylinder_Surface_World_True(5,i);  %G
    IMG_unwrapped(Cylinder_Surface_h_tetha(1,i),Cylinder_Surface_h_tetha(2,i),3)=Cylinder_Surface_World_True(6,i);  %B
end

% IMG Umwrapped------------------------------------------------------------
figure(3)
imshow(IMG_unwrapped);

% IMG Unwrapped Manually Made More Beautiful-------------------------------
figure(4)
imshow(imrotate(circshift(IMG_unwrapped, [0, -unwrap_size_x/2-100]), 180, 'bilinear', 'crop'));














