clc
clear
close all

load('Variables_456.mat')

% DISCLAIMER: THIS CODE TAKES ABOUT 30 MINUTES on an Above Average Computer
% with Intel I7 and 16GB of RAM. This is due to the code given below Part 6

Center_C1_W=Center_C1_W(1:3);
Center_C2_W=Center_C2_W(1:3);

axis_direction_vector = Center_C2_W(1:3,1) - Center_C1_W(1:3,1);
axis_direction_vector=axis_direction_vector./norm(axis_direction_vector);

p_a_W=Coord_pi_to_W*H_rect*p_a;
p_b_W=Coord_pi_to_W*H_rect*p_b;

p_a_W=p_a_W./p_a_W(4);
p_b_W=p_b_W./p_b_W(4);

p_a_W=p_a_W(1:3,1);
p_b_W=p_b_W(1:3,1);

min_err=inf;
for t=0:0.001:1
    Center_C1_W_True=t.*Center_C1_W;

    axis_direction_vector_True = Center_C2_W(1:3,1) - Center_C1_W_True(1:3,1);
    axis_direction_vector_True=axis_direction_vector_True./norm(axis_direction_vector_True);

    % Define the equation of the planes normal to true axis
    a = axis_direction_vector_True(1);
    b = axis_direction_vector_True(2);
    c = axis_direction_vector_True(3);
    d1 = -dot(axis_direction_vector_True,Center_C1_W_True);    %passing through Center_C1_W_True
    d2 = -dot(axis_direction_vector_True,Center_C2_W);         %passing through Center_C2_W
    
    % Projection of the points on the conics to the planes orthogonal to the true axis are their true position
    ta = (-d1) / dot(axis_direction_vector_True, p_a_W);                        % t1*p_a_W is on the plane.
    tb = (-d2) / dot(axis_direction_vector_True, p_b_W);                        % t2*p_b_W is on the plane.

    p_a_W_True = ta * p_a_W;
    p_b_W_True = tb * p_b_W;

    Radius_C1=norm(p_a_W_True - Center_C1_W_True);
    Radius_C2=norm(p_b_W_True - Center_C2_W);

    err=abs(Radius_C1-Radius_C2);
    if err<min_err
        min_err=err;

        t_best=t;
        axis_direction_vector_True_best=axis_direction_vector_True;

        % For projecting to planes orthogonal to axis
        a_best=a;
        b_best=b;
        c_best=c;
        d1_best=d1;
        d2_best=d2;

        ta_best=ta;
        tb_best=tb;
    end
end

min_err
t=t_best;
axis_direction_vector_True=axis_direction_vector_True_best;

a=a_best;
b=b_best;
c=c_best;
d1=d1_best;
d2=d2_best;

ta=ta_best;
tb=tb_best;

Center_C1_W_True=t.* Center_C1_W;

p_a_W_True = ta.* p_a_W;
p_b_W_True = tb.* p_b_W;

%-------------------------------Part 4-------------------------------------
axis_direction_vector_True

%-------------------------------Part 5-------------------------------------
Radius=norm(p_a_W_True - Center_C1_W_True);
Distance=norm(Center_C2_W-Center_C1_W_True);
Radius_over_Distance=Radius/Distance

%------------------------Preamble to Part 6--------------------------------
tic
%Project all points to the 3D Cylinder
options = optimset('Display', 'off');
for i=1:length(Cylinder_Surface_World)
    point=Cylinder_Surface_World(1:3,i);
    fun=@(tcyl) Radius-norm(cross(tcyl*point-Center_C1_W_True,tcyl*point-Center_C2_W))/norm(Center_C2_W-Center_C1_W_True);
    scale_point = fsolve(fun, 1,options);
    point_cyl = scale_point * point;
    Cylinder_Surface_World_True(1:3,i)=point_cyl;
end
toc

Cylinder_Surface_World_True(4:6,:)=Cylinder_Surface_World(4:6,:);
save('Variables_6.mat',"Cylinder_Surface_World_True",'axis_direction_vector_True','Center_C1_W_True','Center_C2_W','Radius','Distance')


% Cylinder Axis and Desired Surface in 3D Space----------------------------
figure(1);
xlabel('X');ylabel("Y");zlabel("Z")
scatter3(0,0,0,'r')
text(0,0,0,'camera')
hold on;
axis equal;

% Initial World Coordinates
scatter3(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'filled','r')
text(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'C1 Center')
hold on
scatter3(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'filled','r')
text(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'C2 Center')
% plot3([Center_C1_W(1) Center_C2_W(1)],[Center_C1_W(2) Center_C2_W(2)],[Center_C1_W(3) Center_C2_W(3)],'color','b','LineWidth',2) %Draw axis
hold on
scatter3(p_a_W(1),p_a_W(2),p_a_W(3),'filled','g')
text(p_a_W(1),p_a_W(2),p_a_W(3),'p_a')
hold on
scatter3(p_b_W(1),p_b_W(2),p_b_W(3),'filled','g')
text(p_b_W(1),p_b_W(2),p_b_W(3),'p_b')
hold on

%Viewing Rays
plot3([0 Center_C1_W(1)/Center_C1_W(3)],[0 Center_C1_W(2)/Center_C1_W(3)],[0 Center_C1_W(3)/Center_C1_W(3)],'color','y','LineWidth',2)
hold on
plot3([0 Center_C2_W(1)/Center_C2_W(3)],[0 Center_C2_W(2)/Center_C2_W(3)],[0 Center_C2_W(3)/Center_C2_W(3)],'color','y','LineWidth',2)
hold on

% True World Coordinates
scatter3(Center_C1_W_True(1),Center_C1_W_True(2),Center_C1_W_True(3),'filled','r')
text(Center_C1_W_True(1),Center_C1_W_True(2),Center_C1_W_True(3),'C1 Center True')
hold on
scatter3(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'filled','r')
text(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'C2 Center True')
hold on
plot3([Center_C1_W_True(1) Center_C2_W(1)],[Center_C1_W_True(2) Center_C2_W(2)],[Center_C1_W_True(3) Center_C2_W(3)],'color','b','LineWidth',2) %Draw axis
hold on
scatter3(p_a_W_True(1),p_a_W_True(2),p_a_W_True(3),'filled','g')
text(p_a_W_True(1),p_a_W_True(2),p_a_W_True(3),'p_a True')
hold on
scatter3(p_b_W_True(1),p_b_W_True(2),p_b_W_True(3),'filled','g')
text(p_b_W_True(1),p_b_W_True(2),p_b_W_True(3),'p_b True')
hold on


% Plot the circle
Radius_C1 = norm(p_a_W_True - Center_C1_W_True);
f_PlotCircle(Radius_C1,Center_C1_W_True,axis_direction_vector_True)
Radius_C2 = norm(p_b_W_True - Center_C2_W);
f_PlotCircle(Radius_C2,Center_C2_W,axis_direction_vector_True)

xlabel('X');ylabel("Y");zlabel("Z");

%--------------------------------------------------------------------------
figure(2);
scatter3(0,0,0,'r')
hold on;
text(0,0,0,'camera')
hold on;
axis equal;

% Initial World Coordinates
scatter3(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'filled','r')
text(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'C1 Center')
hold on
scatter3(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'filled','r')
text(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'C2 Center')
% plot3([Center_C1_W(1) Center_C2_W(1)],[Center_C1_W(2) Center_C2_W(2)],[Center_C1_W(3) Center_C2_W(3)],'color','b','LineWidth',2) %Draw axis
hold on
scatter3(p_a_W(1),p_a_W(2),p_a_W(3),'filled','g')
text(p_a_W(1),p_a_W(2),p_a_W(3),'p_a')
hold on
scatter3(p_b_W(1),p_b_W(2),p_b_W(3),'filled','g')
text(p_b_W(1),p_b_W(2),p_b_W(3),'p_b')
hold on

%Viewing Rays
plot3([0 Center_C1_W(1)/Center_C1_W(3)],[0 Center_C1_W(2)/Center_C1_W(3)],[0 Center_C1_W(3)/Center_C1_W(3)],'color','y','LineWidth',2)
hold on
plot3([0 Center_C2_W(1)/Center_C2_W(3)],[0 Center_C2_W(2)/Center_C2_W(3)],[0 Center_C2_W(3)/Center_C2_W(3)],'color','y','LineWidth',2)
hold on

% True World Coordinates
scatter3(Center_C1_W_True(1),Center_C1_W_True(2),Center_C1_W_True(3),'filled','r')
text(Center_C1_W_True(1),Center_C1_W_True(2),Center_C1_W_True(3),'C1 Center True')
hold on
scatter3(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'filled','r')
text(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'C2 Center True')
hold on
plot3([Center_C1_W_True(1) Center_C2_W(1)],[Center_C1_W_True(2) Center_C2_W(2)],[Center_C1_W_True(3) Center_C2_W(3)],'color','b','LineWidth',2) %Draw axis
hold on
scatter3(p_a_W_True(1),p_a_W_True(2),p_a_W_True(3),'filled','g')
text(p_a_W_True(1),p_a_W_True(2),p_a_W_True(3),'p_a True')
hold on
scatter3(p_b_W_True(1),p_b_W_True(2),p_b_W_True(3),'filled','g')
text(p_b_W_True(1),p_b_W_True(2),p_b_W_True(3),'p_b True')
hold on

% Plot the circle
Radius_C1 = norm(p_a_W_True - Center_C1_W_True);
f_PlotCircle(Radius_C1,Center_C1_W_True,axis_direction_vector_True)
Radius_C2 = norm(p_b_W_True - Center_C2_W);
f_PlotCircle(Radius_C2,Center_C2_W,axis_direction_vector_True)

% Plot points in intial 3D Position
scatter3(Cylinder_Surface_World(1,:),Cylinder_Surface_World(2,:),Cylinder_Surface_World(3,:),1,[0.6350 0.0780 0.1840]);
hold on;

% Plot points mapped on to the 3D Cylindrical Surface
scatter3(Cylinder_Surface_World_True(1,:),Cylinder_Surface_World_True(2,:),Cylinder_Surface_World_True(3,:),1,[0.4660 0.6740 0.1880]);
hold on;


xlabel('X');ylabel("Y");zlabel("Z");








