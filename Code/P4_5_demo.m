clc
clear
close all

load('Variables_456.mat')


Center_C1_W=Center_C1_W(1:3);
Center_C2_W=Center_C2_W(1:3);

p_a_W=Coord_pi_to_W*H_rect*p_a;
p_b_W=Coord_pi_to_W*H_rect*p_b;

p_a_W=p_a_W./p_a_W(4);
p_b_W=p_b_W./p_b_W(4);

p_a_W=p_a_W(1:3,1);
p_b_W=p_b_W(1:3,1);

axis_direction_vector = Center_C2_W(1:3,1) - Center_C1_W(1:3,1);
axis_direction_vector=axis_direction_vector./norm(axis_direction_vector);

% Define the equation of the planes normal to axis
a = axis_direction_vector(1);
b = axis_direction_vector(2);
c = axis_direction_vector(3);
d1 = -a*Center_C1_W(1)-b*Center_C1_W(2)-c*Center_C1_W(3);       %passing through Center_C1_W
d2 = -a*Center_C2_W(1)-b*Center_C2_W(2)-c*Center_C2_W(3);       %passing through Center_C2_W

t1 = (-d1) / dot(axis_direction_vector, p_a_W);                 % t1*p_a_W is on the plane.
t2 = (-d2) / dot(axis_direction_vector, p_b_W);                 % t2*p_b_W is on the plane.

p_a_W_True = t1 * p_a_W;
p_b_W_True = t2 * p_b_W;

Radius_C1 = norm(p_a_W - Center_C1_W);
Radius_C2 = norm(p_b_W - Center_C2_W);




% Cylinder Axis and Desired Surface in 3D Space----------------------------
figure(1);
xlabel('X');ylabel("Y");zlabel("Z")
scatter3(0,0,0,'r')
text(0,0,0,'camera')
hold on;


% Plot centers, axis
scatter3(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'filled','r')
text(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'C1 Center')
hold on
scatter3(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'filled','r')
text(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'C2 Center')
hold on
plot3([Center_C1_W(1) Center_C2_W(1)],[Center_C1_W(2) Center_C2_W(2)],[Center_C1_W(3) Center_C2_W(3)],'color','b','LineWidth',2) %Draw axis
hold on

%Viewing Rays
plot3([0 Center_C1_W(1)/Center_C1_W(3)],[0 Center_C1_W(2)/Center_C1_W(3)],[0 Center_C1_W(3)/Center_C1_W(3)],'color','y','LineWidth',2)
hold on
plot3([0 Center_C2_W(1)/Center_C2_W(3)],[0 Center_C2_W(2)/Center_C2_W(3)],[0 Center_C2_W(3)/Center_C2_W(3)],'color','y','LineWidth',2)
hold on


%--------------------------------------------------------------------------

figure(2);
xlabel('X');ylabel("Y");zlabel("Z")
scatter3(0,0,0,'r')
text(0,0,0,'camera')
hold on;
axis equal;


% Plot centers, axis
scatter3(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'filled','r')
text(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'C1 Center')
hold on
scatter3(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'filled','r')
text(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'C2 Center')
hold on
plot3([Center_C1_W(1) Center_C2_W(1)],[Center_C1_W(2) Center_C2_W(2)],[Center_C1_W(3) Center_C2_W(3)],'color','b','LineWidth',2) %Draw axis
hold on

%Viewing Rays of Centers
plot3([0 Center_C1_W(1)/Center_C1_W(3)],[0 Center_C1_W(2)/Center_C1_W(3)],[0 Center_C1_W(3)/Center_C1_W(3)],'color','y','LineWidth',2)
hold on
plot3([0 Center_C2_W(1)/Center_C2_W(3)],[0 Center_C2_W(2)/Center_C2_W(3)],[0 Center_C2_W(3)/Center_C2_W(3)],'color','y','LineWidth',2)
hold on

% Points a b
scatter3(p_a_W(1),p_a_W(2),p_a_W(3),'filled','g')
text(p_a_W(1),p_a_W(2),p_a_W(3),'p_a')
hold on
scatter3(p_b_W(1),p_b_W(2),p_b_W(3),'filled','g')
text(p_b_W(1),p_b_W(2),p_b_W(3),'p_b')
hold on

%Viewing Rays of points a and b
plot3([0 p_a_W(1)/p_a_W(3)],[0 p_a_W(2)/p_a_W(3)],[0 p_a_W(3)/p_a_W(3)],'color','y','LineWidth',2)
hold on
plot3([0 p_b_W(1)/p_b_W(3)],[0 p_b_W(2)/p_b_W(3)],[0 p_b_W(3)/p_b_W(3)],'color','y','LineWidth',2)
hold on

% Cross Sectional Planes
[x_plot, y_plot] = meshgrid(-0.4:0.005:0.4, -0.4:0.005:0.4);
z_plot_1 = -(a * x_plot + b * y_plot + d1) / c; % Calculate the corresponding z values for each point using the plane equation
surf(x_plot, y_plot, z_plot_1, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;

z_plot_2 = -(a * x_plot + b * y_plot + d2) / c; % Calculate the corresponding z values for each point using the plane equation
surf(x_plot, y_plot, z_plot_2, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;

% Plot true p_a_W, p_b_W
scatter3(p_a_W_True(1),p_a_W_True(2),p_a_W_True(3),'filled','g')
text(p_a_W_True(1),p_a_W_True(2),p_a_W_True(3),'p_a True')
hold on
scatter3(p_b_W_True(1),p_b_W_True(2),p_b_W_True(3),'filled','g')
text(p_b_W_True(1),p_b_W_True(2),p_b_W_True(3),'p_b True')
hold on


%--------------------------------------------------------------------------

figure(3);
xlabel('X');ylabel("Y");zlabel("Z")
scatter3(0,0,0,'r')
text(0,0,0,'camera')
hold on;
axis equal;


% Plot centers, axis
scatter3(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'filled','r')
text(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'C1 Center')
hold on
scatter3(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'filled','r')
text(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'C2 Center')
hold on
plot3([Center_C1_W(1) Center_C2_W(1)],[Center_C1_W(2) Center_C2_W(2)],[Center_C1_W(3) Center_C2_W(3)],'color','b','LineWidth',2) %Draw axis
hold on

%Viewing Rays of Centers
plot3([0 Center_C1_W(1)/Center_C1_W(3)],[0 Center_C1_W(2)/Center_C1_W(3)],[0 Center_C1_W(3)/Center_C1_W(3)],'color','y','LineWidth',2)
hold on
plot3([0 Center_C2_W(1)/Center_C2_W(3)],[0 Center_C2_W(2)/Center_C2_W(3)],[0 Center_C2_W(3)/Center_C2_W(3)],'color','y','LineWidth',2)
hold on

% Points a b
scatter3(p_a_W(1),p_a_W(2),p_a_W(3),'filled','g')
text(p_a_W(1),p_a_W(2),p_a_W(3),'p_a')
hold on
scatter3(p_b_W(1),p_b_W(2),p_b_W(3),'filled','g')
text(p_b_W(1),p_b_W(2),p_b_W(3),'p_b')
hold on

% Cross Sectional Planes
[x_plot, y_plot] = meshgrid(-0.4:0.005:0.4, -0.4:0.005:0.4);
z_plot_1 = -(a * x_plot + b * y_plot + d1) / c; % Calculate the corresponding z values for each point using the plane equation
surf(x_plot, y_plot, z_plot_1, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;

z_plot_2 = -(a * x_plot + b * y_plot + d2) / c; % Calculate the corresponding z values for each point using the plane equation
surf(x_plot, y_plot, z_plot_2, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;

% Plot true p_a_W, p_b_W
scatter3(p_a_W_True(1),p_a_W_True(2),p_a_W_True(3),'filled','g')
text(p_a_W_True(1),p_a_W_True(2),p_a_W_True(3),'p_a True')
hold on
scatter3(p_b_W_True(1),p_b_W_True(2),p_b_W_True(3),'filled','g')
text(p_b_W_True(1),p_b_W_True(2),p_b_W_True(3),'p_b True')
hold on


% Plot the circles
f_PlotCircle(Radius_C1,Center_C1_W,axis_direction_vector)
f_PlotCircle(Radius_C2,Center_C2_W,axis_direction_vector)



