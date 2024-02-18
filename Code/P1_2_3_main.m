clc
clear
close all


%--------------------------Setup Features----------------------------------
% Load the Image
    IMG = imrotate(imread('PalazzoTe.jpg'), -90, 'bilinear'); % pixel*pixel*3 -> rgb + Image is loaded rotated when downloading for some reason, fix that.
    size_img=size(IMG);

% Plotting Variables
    size_y=size_img(1);
    size_x=size_img(2);
    
    x_range=0:size_x;
    y_range=0:size_y;
    
    x_plot=0:size_x;
    FNT_SZ=29;

% Define the Extracted/Selected Features l1, l2, C1, and C2.
    l1 = [-1.661043177892919e+03;-29.982728842832444;4.718827791588740e+05];
    l2 = [-1.049395509499137e+03;-8.095336787564770e+02;5.440922874469413e+06];
    C1 = [1.106970446841802e-07,-6.770411107233868e-08,-4.730976414356710e-05;-6.770411107233868e-08,2.602117671521042e-07,-4.983472635663797e-04;-4.730976414356710e-05,-4.983472635663797e-04,0.999999498823408];
    C2 = [6.162217728832053e-08,-1.632930635820811e-09,-7.261728892754088e-05;-1.632930635820811e-09,6.720519975155384e-08,-2.579239239372433e-04;-7.261728892754088e-05,-2.579239239372433e-04,0.999999856403943];


    % Depending on the selection/extraction code, these values may not be normalized. Normalize them just in case.
    l1 = l1./l1(3);
    l2 = l2./l2(3);
    C1 = C1./norm(C1);
    C2 = C2./norm(C2);

    [C1_a, C1_b, C1_c, C1_d, C1_e, C1_f] = deal(C1(1), C1(2)*2, C1(5), C1(3)*2, C1(6)*2, C1(9));
    [C2_a, C2_b, C2_c, C2_d, C2_e, C2_f] = deal(C2(1), C2(2)*2, C2(5), C2(3)*2, C2(6)*2, C2(9));




%-------------------------------Part 1-------------------------------------
% Find the intersection points I,J of C1 and C2 (Images of circular points)

% Symbolic variables
    syms x y

% Define the equations for C1 and C2
    C1_eq_syms = C1_a*x^2 + C1_b*x*y + C1_c*y^2 + C1_d*x + C1_e*y + C1_f;
    C2_eq_syms = C2_a*x^2 + C2_b*x*y + C2_c*y^2 + C2_d*x + C2_e*y + C2_f;
    eqns = [C1_eq_syms ==0, C2_eq_syms ==0];
    S_circ = solve(eqns, [x,y]);
    Intersection_of_C1_and_C2=[S_circ.x,S_circ.y]

%  If there is a unique pair of complex conjugates, the image of circular
%  points I and J  are exactly this pair of complex conjugates  (Single Axis
%  Geometry by Fitting Conics by Jiang et al)
    I = [double(S_circ.x(1));double(S_circ.y(1));1]
    J = [double(S_circ.x(2));double(S_circ.y(2));1]
    I_scene=[1;i;0];
    J_scene=[1;-i;0];

% The vanishing line (horizon) passes through I and J.
% The horizon of circumfrences in same or parallel planes: The horizon of
% these planes passes through I, J images.
% (Since I_scene, J_scene are on l_inf and l_inf maps onto h of circumfrences.
    h=cross(I,J);
    h = h./norm(h)  % Normalize for numerical stability




%-------------------------------Part 2-------------------------------------
% Cylinder's axis is the line connecting the centers of the conics.
% Find the centers.
    Center_C1 = [
    (C1_b * C1_e - 2 * C1_c * C1_d) / (4 * C1_a * C1_c - C1_b^2);
    (C1_b * C1_d - 2 * C1_a * C1_e) / (4 * C1_a * C1_c - C1_b^2);
    1];

    Center_C2 = [
    (C2_b * C2_e - 2 * C2_c * C2_d) / (4 * C2_a * C2_c - C2_b^2);
    (C2_b * C2_d - 2 * C2_a * C2_e) / (4 * C2_a * C2_c - C2_b^2);
    1];


% Find the cylinder axis that passes through the centers
    a=cross(Center_C1, Center_C2);          %Cylinder Axis
    a = a./norm(a)

% Find the vanishing point. Vanishing point of the axis is a point where it
% intersects parallel lines on the image. (Where it intersects l1 or l2.)
    V_1 = cross(a,l1);
    V_1 = V_1/V_1(3)



%-------------------------Preamble to Part 3-------------------------------
% We need more vanishing points. Define l3 - l6. By defining points a-d

% Point a: l1 x C1
    l1_eq_syms = l1(1)*x + l1(2)*y + l1(3);
    eqns = [C1_eq_syms ==0, l1_eq_syms ==0];
    S = solve(eqns, [x,y]);

    p_a_1 = [double(S.x(1));double(S.y(1));1];
    p_a_2 = [double(S.x(2));double(S.y(2));1];
    if(p_a_1(2)<p_a_2(2))   %We want the top intersection (smaller y value wrt graphical origin: top left corner)
        p_a=p_a_1;
    else
        p_a=p_a_2;
    end

% Point b: l1 x C2
    eqns = [C2_eq_syms ==0, l1_eq_syms ==0];
    S = solve(eqns, [x,y]);

    p_b_1 = [double(S.x(1));double(S.y(1));1];
    p_b_2 = [double(S.x(2));double(S.y(2));1];
    if(p_b_1(2)<p_b_2(2))   %We want the top intersection (smaller y value wrt graphical origin)
        p_b=p_b_1;
    else
        p_b=p_b_2;
    end

% Point c: l2 x C1
    l2_eq_syms = l2(1)*x + l2(2)*y + l2(3);
    eqns = [C1_eq_syms ==0, l2_eq_syms ==0];
    S = solve(eqns, [x,y]);

    p_c_1 = [double(S.x(1));double(S.y(1));1];
    p_c_2 = [double(S.x(2));double(S.y(2));1];
    if(p_c_1(2)<p_c_2(2))   %We want the top intersection (smaller y value wrt graphical origin)
        p_c=p_c_1;
    else
        p_c=p_c_2;
    end

% Point d: l2 x C2
    eqns = [C2_eq_syms ==0, l2_eq_syms ==0];
    S = solve(eqns, [x,y]);

    p_d_1 = [double(S.x(1));double(S.y(1));1];
    p_d_2 = [double(S.x(2));double(S.y(2));1];
    if(p_d_1(2)<p_d_2(2))   %We want the top intersection (smaller y value wrt graphical origin)
        p_d=p_d_1;
    else
        p_d=p_d_2;
    end

% Line l3 = p_a x p_c
    l3=cross(p_a, p_c);
    l3 = l3./norm(l3);
    l3_eq_syms = l3(1)*x + l3(2)*y + l3(3);

% Line l4 = p_b x p_d
    l4=cross(p_b, p_d);
    l4 = l4./norm(l4);
    l4_eq_syms = l4(1)*x + l4(2)*y + l4(3);

% Line l5 = center_C1 x p_c
    l5=cross(Center_C1, p_c);
    l5 = l5./norm(l5);

% Line l6 = center_C2 x p_d
    l6=cross(Center_C2, p_d);
    l6 = l6./norm(l6);


% We need more vanishing points. Images of parallel lines meet at the vanishing points.
%l3 and l4 are parallel. They meet at V_2
    V_2 = cross(l3,l4);
    V_2 = V_2/V_2(3);

%l5 and l6 are parallel. They meet at V_3
    V_3 = cross(l5,l6);
    V_3 = V_3/V_3(3);



%-------------------------------Part 3-------------------------------------
% Rectify the image wrt. the circular cross sections. Wrt the cylinder caps.
CDCP_img = I*J.' + J*I.'; %Image of the Conic Dual to the Circular Points
CDCP_img=CDCP_img./norm(CDCP_img);

% Be sure that we can get an accurate rectification using the CDCP method:
[V_eigs_CDCP,D_eigs_CDCP]=eigs(CDCP_img);
for i=1:length(D_eigs_CDCP)
    if(D_eigs_CDCP(i,i)<0)
        disp("Rectification Method Using CDCP will be inaccurate!");
    end
end

% No problems above, continue by constructing the rectifying homography.
[U,D,U_t] = svd(CDCP_img);
H_rect = [ sqrt((D(1,1)^-1)) 0 0  ;  0 sqrt((D(2,2)^-1)) 0  ;  0 0 1 ] * U';
H_rect = H_rect / H_rect(3, 3)

% Now let's find the K matrix. To find the K matrix we can find the image
% of the absolute conic w. w of a planar object is related to the rectifing
% homography of that planar object such that: 
% h1'w*h2=0 AND h1'*w*h1 - h2'*w*h2=0. Where [h1 h2 h3] = H_rect^-1
% These equations are derived from the invariant effect of similarity on
% circular points. And H_rect must be a similarity rectifying homography.
% H_rect is such an homography.

% Secondly in an image, w relates the vanishing points of orthogonal
% directions in any plane as such: v1'*w*v2=0
% We know that in the planar image defined by l1, l2, l3, l4 we have 2 such
% vanishing points defined by V_1 and V_2: V_1'*w*V_2=0
% We also know that in the planar image defined by a, l5, l2, l6 we have
% another couple of vanishing points defined by V_1 and V_3: V_1'*w*V_3=0
% Please note that any 2 vanishing points defined by orthogonal parallel
% line sets work as the w gets involved simply due to geometric
% definitions.

% To find w we needed 4 equations and now we have those 4 equations:
% h1'*w*h2=0
% h1'*w*h1 - h2'*w*h2=0
% V_1'*w*V_2=0
% V_1'*w*V_3=0
% Solve for w.

H_rect_inv=inv(H_rect)
h1=H_rect_inv(:,1);
h2=H_rect_inv(:,2);
h3=H_rect_inv(:,3);


syms w11 w13 w23 w33

% Define the equations with the given conditions
eq1 = h1.' * [w11, 0, w13; 0, w11, w23; w13, w23, w33] * h2 == 0;
eq2 = h1.' * [w11, 0, w13; 0, w11, w23; w13, w23, w33] * h1 - h2.' * [w11, 0, w13; 0, w11, w23; w13, w23, w33] * h2 == 0;
eq3 = V_1.' * [w11, 0, w13; 0, w11, w23; w13, w23, w33] * V_2 == 0;
eq4 = V_1.' * [w11, 0, w13; 0, w11, w23; w13, w23, w33] * V_3 == 0;

[A,b] = equationsToMatrix([eq1, eq2, eq3],[w11, w13, w23, w33]);
w_values= null(A);

w=[w_values(1), 0, w_values(2); 0, w_values(1), w_values(3); w_values(2), w_values(3), w_values(4)];


% K can be found from the Cholesky factorization of w:
K=double(inv(chol(w)));
K=K./K(3,3)     %Normalize



%-----------------------Preamble to Part 4,5-------------------------------
Q=inv(K)*H_rect_inv;

i_pi=Q(:,1)
j_pi=Q(:,2)
o_pi=Q(:,3)

Coord_pi_to_W=[i_pi j_pi o_pi;      %Relates the coordinate systems from a planar image (plane pi / rectified image) to the world.
                0    0    1];       %Is the pose of the planar images (Rectified cross sections) which is the same as saying this
                                    %is the pose of the normal to the planar image i.e. this is the pose of the cylinder axis in world.

C1_pi=H_rect_inv'*C1*H_rect_inv;    %Rectify Conics
C2_pi=H_rect_inv'*C2*H_rect_inv;
C1_pi=C1_pi./norm(C1_pi);           %Normalize to Homogenous Coords
C2_pi=C2_pi./norm(C2_pi);
[C1_pi_a, C1_pi_b, C1_pi_c, C1_pi_d, C1_pi_e, C1_pi_f] = deal(C1_pi(1), C1_pi(2)*2, C1_pi(5), C1_pi(3)*2, C1_pi(6)*2, C1_pi(9));
[C2_pi_a, C2_pi_b, C2_pi_c, C2_pi_d, C2_pi_e, C2_pi_f] = deal(C2_pi(1), C2_pi(2)*2, C2_pi(5), C2_pi(3)*2, C2_pi(6)*2, C2_pi(9));


% Pi plane coordinates of conic centers (Rectified image plane)
% We want the world coordinates of conic centers. We know how to get it after getting the pi plane (Rectified planar image) coordinates
% Conicse Change with a special formula we cant get Center_pi from directly applying the homography to the original image centers.
% We must use the center formula with the new conics instead.
    Center_C1_pi= [
    (C1_pi_b * C1_pi_e - 2 * C1_pi_c * C1_pi_d) / (4 * C1_pi_a * C1_pi_c - C1_pi_b^2);
    (C1_pi_b * C1_pi_d - 2 * C1_pi_a * C1_pi_e) / (4 * C1_pi_a * C1_pi_c - C1_pi_b^2);
    1];

    Center_C2_pi= [
    (C2_pi_b * C2_pi_e - 2 * C2_pi_c * C2_pi_d) / (4 * C2_pi_a * C2_pi_c - C2_pi_b^2);
    (C2_pi_b * C2_pi_d - 2 * C2_pi_a * C2_pi_e) / (4 * C2_pi_a * C2_pi_c - C2_pi_b^2);
    1];


% Get the real world coordinates of the conic centers
Center_C1_W=Coord_pi_to_W*Center_C1_pi;
Center_C2_W=Coord_pi_to_W*Center_C2_pi;
Center_C1_W=Center_C1_W./Center_C1_W(4);
Center_C2_W=Center_C2_W./Center_C2_W(4);





%--------------------------Preamble to Part 6------------------------------
%Firstly extract the point coordinates that define the cylindrical surface.
%Cylinder_Surface_Image is 6 row x Extracted_Pixel_Count column matrix
%Holding x,y,1,R,G,B at each row.

% Put all (x,y) pixel combinations to 2 column vectors.
[x_matrix, y_matrix] = meshgrid(x_range(2:end), y_range(2:end)); %Dont take the points with x or y =0, problem in indexing later on.
x_values = reshape(x_matrix, [], 1);
y_values = reshape(y_matrix, [], 1);

l1_result_values=l1(1)*x_values+l1(2)*y_values+l1(3);
l2_result_values=l2(1)*x_values+l2(2)*y_values+l2(3);
l3_result_values=l3(1)*x_values+l3(2)*y_values+l3(3);
l4_result_values=l4(1)*x_values+l4(2)*y_values+l4(3);
C1_result_values=C1_a*x_values.^2 + C1_b*x_values.*y_values + C1_c*y_values.^2 + C1_d*x_values + C1_e*y_values + C1_f;
C2_result_values=C2_a*x_values.^2 + C2_b*x_values.*y_values + C2_c*y_values.^2 + C2_d*x_values + C2_e*y_values + C2_f;

k=1;    %Counter for the extracted points
for i=1:length(x_values)
    x=x_values(i);
    y=y_values(i);
    l1_result=l1_result_values(i);
    l2_result=l2_result_values(i);
    l3_result=l3_result_values(i);
    l4_result=l4_result_values(i);
    C1_result=C1_result_values(i);
    C2_result=C2_result_values(i);

    if l1_result<=0 && l2_result>=0 && C2_result>=0         %Between the generatrix lines and outside C2.
        if (C1_result<=0) || (l3_result>=0 && l4_result<=0) %AND either in C1 or in abcd.
            Cylinder_Surface_Image(1,k)=x;                  %Point x
            Cylinder_Surface_Image(2,k)=y;                  %Point y
            Cylinder_Surface_Image(3,k)=1;                  %Point z=1 (Homogenous)
            Cylinder_Surface_Image(4,k)=IMG(y,x,1);         %R
            Cylinder_Surface_Image(5,k)=IMG(y,x,2);         %G
            Cylinder_Surface_Image(6,k)=IMG(y,x,3);         %B
            k=k+1;
        end
    end

end

Cylinder_Surface_Rectified=H_rect*Cylinder_Surface_Image(1:3,:);    % Rectify
Cylinder_Surface_World=Coord_pi_to_W*Cylinder_Surface_Rectified;    % 3D World Coordinates. Must fix them so they look like X Y Z 1.
Cylinder_Surface_World=Cylinder_Surface_World(1:4,:)./Cylinder_Surface_World(4,:);  % Normalize so that the last element is 1
Cylinder_Surface_World(4:6,:)=Cylinder_Surface_Image(4:6,:);        % Add back in the RGB values.


save('Variables_456.mat', 'Center_C1_W', 'Center_C2_W', 'Coord_pi_to_W', 'Cylinder_Surface_World', 'H_rect', 'H_rect_inv', 'p_a', 'p_b');



%----------------------------Plot Features---------------------------------
% Original Image with Useful Features--------------------------------------
figure(1);
imshow(IMG);
hold on;

% l1 (y=(-ax-c)/b)
    plot( [V_1(1) , - l1(3)/l1(1)] , [V_1(2) , 0], 'LineWidth', 2, 'Color', 'b'); %Plot between V and top of image.
    hold on;
    text(300, 725, 'l1', 'FontSize', FNT_SZ, 'Color', 'b');
    hold on;


% l2
    plot([V_1(1) , x_plot(end)] , [V_1(2) , (-l2(1)*x_plot(end) - l2(3))/l2(2)], 'LineWidth', 2, 'Color', 'b');
    hold on;
    text(3167, 2700, 'l2', 'FontSize', FNT_SZ, 'Color', 'b');
    hold on;

% C1
    C1_eq = @(x_var, y_var) C1_a*x_var.^2 + C1_b*x_var.*y_var + C1_c*y_var.^2 + C1_d*x_var + C1_e*y_var + C1_f;
    graph=ezplot(C1_eq, [0, size_x, 0, size_y]);
    set(graph, 'LineWidth', 2, 'Color', 'b');
    hold on;
    text(2180, 1300, 'C1', 'FontSize', FNT_SZ, 'Color', 'b');
    hold on;

% C2
    C2_eq = @(x_var, y_var) C2_a*x_var.^2 + C2_b*x_var.*y_var + C2_c*y_var.^2 + C2_d*x_var + C2_e*y_var + C2_f;
    graph=ezplot(C2_eq, [0, size_x, 0, size_y]);
    set(graph, 'LineWidth', 2, 'Color', 'b');
    text(1694, 2600, 'C2', 'FontSize', FNT_SZ, 'Color', 'b');
    hold on;


set(gcf, 'PaperPositionMode', 'auto');
frame = getframe(gca);
IMG_features = frame.cdata;


% I', J' of conics.
    plot(I(1), I(2), 'mo', 'MarkerSize', 10);
    text(I(1) + 50, I(2), 'I_{Cross}', 'FontSize', FNT_SZ, 'Color', 'm');

    plot(J(1), J(2), 'mo', 'MarkerSize', 10);
    text(J(1) + 50, J(2), 'J_{Cross}', 'FontSize', FNT_SZ, 'Color', 'm');

% h: Horizon of cross sectional planes.
    y_plot = (-h(1)*x_plot - h(3))/h(2);   %y=(-ax-c)/b
    plot(x_plot,y_plot, 'LineWidth', 2, 'Color', 'c');
    hold on;
    text(3340, 200, 'h', 'FontSize', FNT_SZ, 'Color', 'c');
    hold on;



% Centers of C1 and C2
    plot(Center_C1(1), Center_C1(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    text(Center_C1(1) + 50, Center_C1(2), 'Center C1', 'FontSize', FNT_SZ, 'Color', 'r');

    plot(Center_C2(1), Center_C2(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    text(Center_C2(1) + 50, Center_C2(2), 'Center C2', 'FontSize', FNT_SZ, 'Color', 'r');

% a: Cylinder Axis
    plot([V_1(1) , - a(3)/a(1)] , [V_1(2) , 0], 'LineWidth', 2, 'Color', 'r');
    hold on;

% V: Vanishing Point of cylinder axis. (And parallel lines)
    plot(V_1(1), V_1(2), 'cp', 'MarkerSize', 10);
    text(V_1(1) + 50, V_1(2), 'V_1', 'FontSize', FNT_SZ, 'Color', 'c');




% Point a: l1 x C1
    plot(p_a(1), p_a(2), 'go', 'MarkerSize', 10);
    text(p_a(1) + 50, p_a(2), 'p_a', 'FontSize', FNT_SZ, 'Color', 'g');

% Point b: l1 x C2
    plot(p_b(1), p_b(2), 'go', 'MarkerSize', 10);
    text(p_b(1) + 50, p_b(2), 'p_b', 'FontSize', FNT_SZ, 'Color', 'g');

% Point c: l2 x C1
    plot(p_c(1), p_c(2), 'go', 'MarkerSize', 10);
    text(p_c(1) + 50, p_c(2), 'p_c', 'FontSize', FNT_SZ, 'Color', 'g');

% Point d: l2 x C2
    plot(p_d(1), p_d(2), 'go', 'MarkerSize', 10);
    text(p_d(1) - 200, p_d(2), 'p_d', 'FontSize', FNT_SZ, 'Color', 'g');

% Line l3 = p_a x p_c
    plot([p_a(1) , V_2(1)],[p_a(2) , V_2(2)], 'LineWidth', 2, 'Color', 'g');
    hold on;
    text(1802, 1850, 'l3', 'FontSize', FNT_SZ, 'Color', 'g');
    hold on;

% Line l4 = p_b x p_d
    plot([p_b(1) , V_2(1)],[p_b(2) , V_2(2)], 'LineWidth', 2, 'Color', 'g');
    hold on;
    text(1300, 3300, 'l4', 'FontSize', FNT_SZ, 'Color', 'g');
    hold on;

% Line l5 = center_C1 x p_c
    plot([Center_C1(1) , V_3(1)],[Center_C1(2) , V_3(2)], 'LineWidth', 2, 'Color', 'y');
    hold on;
    text(1300, 2400, 'l5', 'FontSize', FNT_SZ, 'Color', 'y');
    hold on;

% Line l6 = center_C2 x p_d
    plot([Center_C2(1) , V_3(1)],[Center_C2(2) , V_3(2)], 'LineWidth', 2, 'Color', 'y');
    hold on;
    text(850, 3800, 'l6', 'FontSize', FNT_SZ, 'Color', 'y');
    hold on;

% V_2 = l3 x l4
    plot(V_2(1), V_2(2), 'cp', 'MarkerSize', 10);
    text(V_2(1) + 50, V_2(2), 'V_2', 'FontSize', FNT_SZ, 'Color', 'c');

% V_3 = l5 x l6
    plot(V_3(1), V_3(2), 'cp', 'MarkerSize', 10);
    text(V_3(1) + 50, V_3(2), 'V_3', 'FontSize', FNT_SZ, 'Color', 'c');

% Camera's Principal Point
    plot(K(1,3), K(2,3), 'cp', 'MarkerSize', 10);
    text(K(1,3) + 50, K(2,3), '(U_0,V_0)', 'FontSize', FNT_SZ, 'Color', 'c');



% Rectified Image----------------------------------------------------------
figure(2);
title('Image Rectified wrt. Cylinder Caps');
IMG_rectified=imtransform(IMG,maketform( 'projective', H_rect'),'XYScale',1);
imshow(IMG_rectified); %Show the similarity rectified image.
hold on;

% Plotting Variables
    size_img_rectified=size(IMG_rectified);
    size_y_rectified=size_img_rectified(1);
    size_x_rectified=size_img_rectified(2);


% % C1_pi
%     C1_pi_eq = @(x_var, y_var) C1_pi_a*x_var.^2 + C1_pi_b*x_var.*y_var + C1_pi_c*y_var.^2 + C1_pi_d*x_var + C1_pi_e*y_var + C1_pi_f;
%     graph=ezplot(C1_pi_eq, [-size_x_rectified, size_x_rectified, -size_y_rectified, size_y_rectified]);
%     hold on;
%     set(graph, 'LineWidth', 2, 'Color', 'g');
%     hold on;
% 
% % C2_pi
%     C2_pi_eq = @(x_var, y_var) C2_pi_a*x_var.^2 + C2_pi_b*x_var.*y_var + C2_pi_c*y_var.^2 + C2_pi_d*x_var + C2_pi_e*y_var + C2_pi_f;
%     graph=ezplot(C2_pi_eq, [-size_x_rectified, size_x_rectified, -size_y_rectified, size_y_rectified]);
%     set(graph, 'LineWidth', 2, 'Color', 'g');
%     hold on;
% 
% % Centers of C1_pi and C2_pi
%     plot(Center_C1_pi(1), Center_C1_pi(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
%     text(Center_C1_pi(1) + 50, Center_C1_pi(2), 'Center C1', 'FontSize', FNT_SZ, 'Color', 'r');
% 
%     plot(Center_C2_pi(1), Center_C2_pi(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
%     text(Center_C2_pi(1) + 50, Center_C2_pi(2), 'Center C2', 'FontSize', FNT_SZ, 'Color', 'r');


% !!! Conics Must be Manually Shifted to correspond to the automatic shift of the image caused by imtransform!!!
% Please note that we are only shifting the conics and not changing their radii. Therefore it can be seen that the values used in Part E are correct.
% Shift values (Hand Tuned)
shift_x_positive = 6700;
shift_y_positive = 12850;

% Plot shifted C1_pi
    C1_pi_eq_shifted = @(x_var, y_var) C1_pi_a*(x_var - shift_x_positive).^2 + C1_pi_b*(x_var - shift_x_positive).*(y_var - shift_y_positive) + C1_pi_c*(y_var - shift_y_positive).^2 + C1_pi_d*(x_var - shift_x_positive) + C1_pi_e*(y_var - shift_y_positive) + C1_pi_f;
    graph = ezplot(C1_pi_eq_shifted, [-size_x_rectified + shift_x_positive, size_x_rectified + shift_x_positive, -size_y_rectified + shift_y_positive, size_y_rectified + shift_y_positive]);
    hold on;
    set(graph, 'LineWidth', 2, 'Color', 'g');
    hold on;

% Plot shifted C2_pi
    C2_pi_eq_shifted = @(x_var, y_var) C2_pi_a*(x_var - shift_x_positive).^2 + C2_pi_b*(x_var - shift_x_positive).*(y_var - shift_y_positive) + C2_pi_c*(y_var - shift_y_positive).^2 + C2_pi_d*(x_var - shift_x_positive) + C2_pi_e*(y_var - shift_y_positive) + C2_pi_f;
    graph = ezplot(C2_pi_eq_shifted, [-size_x_rectified + shift_x_positive, size_x_rectified + shift_x_positive, -size_y_rectified + shift_y_positive, size_y_rectified + shift_y_positive]);
    set(graph, 'LineWidth', 2, 'Color', 'g');
    hold on;

% Centers of C1_pi and C2_pi
    plot(Center_C1_pi(1)+shift_x_positive, Center_C1_pi(2)+shift_y_positive, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    plot(Center_C2_pi(1)+shift_x_positive, Center_C2_pi(2)+shift_y_positive, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
   






% Cylinder Axis and Desired Surface in 3D Space----------------------------
figure(3);
scatter3(Cylinder_Surface_World(1,:),Cylinder_Surface_World(2,:),Cylinder_Surface_World(3,:));
hold on;
scatter3(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'filled')
text(Center_C1_W(1),Center_C1_W(2),Center_C1_W(3),'C1 Center')
hold on
scatter3(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'filled')
text(Center_C2_W(1),Center_C2_W(2),Center_C2_W(3),'C2 Center')
scatter3(0,0,0,'r')
text(0,0,0,'camera')
plot3([Center_C1_W(1) Center_C2_W(1)],[Center_C1_W(2) Center_C2_W(2)],[Center_C1_W(3) Center_C2_W(3)],'color','b','LineWidth',2) %Draw axis
xlabel('X');ylabel("Y");zlabel("Z")






% Cylinder Surface to be Unwrapped Original--------------------------------
figure(4)
imshow(IMG);
hold on;
scatter(Cylinder_Surface_Image(1,:),Cylinder_Surface_Image(2,:))

% l1 (y=(-ax-c)/b)
    plot( [V_1(1) , - l1(3)/l1(1)] , [V_1(2) , 0], 'LineWidth', 2, 'Color', 'b'); %Plot between V and top of image.
    hold on;

% l2
    plot([V_1(1) , x_plot(end)] , [V_1(2) , (-l2(1)*x_plot(end) - l2(3))/l2(2)], 'LineWidth', 2, 'Color', 'b');
    hold on;

% C1
    graph=ezplot(C1_eq, [0, size_x, 0, size_y]);
    set(graph, 'LineWidth', 2, 'Color', 'b');
    hold on;

% C2
    graph=ezplot(C2_eq, [0, size_x, 0, size_y]);
    set(graph, 'LineWidth', 2, 'Color', 'b');


% Cylinder Surface to be Unwrapped Rectified-------------------------------
figure(5);
title('Image Rectified wrt. Cylinder Caps');
imshow(IMG_rectified); %Show the similarity rectified image.
hold on;

% Rectified Surface
Cylinder_Surface_Rectified=Cylinder_Surface_Rectified(1:3,:)./Cylinder_Surface_Rectified(3,:);
scatter(Cylinder_Surface_Rectified(1,:)+shift_x_positive,Cylinder_Surface_Rectified(2,:)+shift_y_positive)

% Plot shifted C1_pi
graph = ezplot(C1_pi_eq_shifted, [-size_x_rectified + shift_x_positive, size_x_rectified + shift_x_positive, -size_y_rectified + shift_y_positive, size_y_rectified + shift_y_positive]);
hold on;
set(graph, 'LineWidth', 2, 'Color', 'g');
hold on;

% Plot shifted C2_pi
graph = ezplot(C2_pi_eq_shifted, [-size_x_rectified + shift_x_positive, size_x_rectified + shift_x_positive, -size_y_rectified + shift_y_positive, size_y_rectified + shift_y_positive]);
set(graph, 'LineWidth', 2, 'Color', 'g');
hold on;

% Centers of C1_pi and C2_pi
plot(Center_C1_pi(1)+shift_x_positive, Center_C1_pi(2)+shift_y_positive, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(Center_C2_pi(1)+shift_x_positive, Center_C2_pi(2)+shift_y_positive, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    




