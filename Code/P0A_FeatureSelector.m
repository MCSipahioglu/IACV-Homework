clc
clear
close all



% Load the Image
    I = imrotate(imread('PalazzoTe.jpg'), -90, 'bilinear'); % pixel*pixel*3 -> rgb + Image is loaded rotated for some reason, fix that.
    size_img=size(I);

% Plotting Variables
    size_y=size_img(1)
    size_x=size_img(2)

    x_range=0:size_x;
    y_range=0:size_y;
    
    x_plot=0:size_x;
    FNT_SZ=29;


%------------------------Choose Generatrix Lines---------------------------
% Choose 2 points for Generatrix l1. Click on plot, hit enter when you are done.
    figure(1), imshow(I);
    title("Choose 2 Points For l1, Press Enter When Done");
    hold on;
    [x, y] = getpts();

% Save points in the homogeneous coordinates. It is enough to set to 1 the third component
    temp_p1 = [x(1); y(1); 1];
    temp_p2 = [x(2); y(2); 1];

% Draw l1 defined by 2 points
    l1 = cross(temp_p1, temp_p2)    % l1=temp_p1 x temp_p2

    text(temp_p1(1), temp_p1(2), 'l1', 'FontSize', FNT_SZ, 'Color', 'b');
    hold on;
    y_plot = (-l1(1)*x_plot - l1(3))/l1(2);   %y=(-ax-c)/b
    plot(x_plot,y_plot, 'LineWidth', 2, 'Color', 'b');
    hold on;

% Choose 2 more points for Generatrix l2. Click on plot, hit enter when you are done.
    title("Choose 2 Points For l2, Press Enter When Done");
    [x, y] = getpts();
    temp_p1 = [x(1); y(1); 1];
    temp_p2 = [x(2); y(2); 1];


% Draw l2 defined by 2 points
    l2 = cross(temp_p1, temp_p2)    % l2=temp_p1 x temp_p2

    % Draw l2
    text(temp_p1(1), temp_p1(2), 'l2', 'FontSize', FNT_SZ, 'Color', 'b');
    hold on;
    y_plot = (-l2(1)*x_plot - l2(3))/l2(2);   %y=(-ax-c)/b
    plot(x_plot,y_plot, 'LineWidth', 2, 'Color', 'b');
    hold on;


%----------------------------Choose Conics---------------------------------
% Choose 5 points for C1.
    title("Choose 5 Points For C1, Press Enter When Done");
    [x, y]=getpts;

    A=[x.^2 x.*y y.^2 x y ones(size(x))];
    N = null(A);
    cc = N(:, 1);
    [C1_a, C1_b, C1_c, C1_d, C1_e, C1_f] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));
    C1=[C1_a C1_b/2 C1_d/2; C1_b/2 C1_c C1_e/2; C1_d/2 C1_e/2 C1_f]

%Draw C1
    C1_eq = @(x_var, y_var) C1_a*x_var.^2 + C1_b*x_var.*y_var + C1_c*y_var.^2 + C1_d*x_var + C1_e*y_var + C1_f;
    text(x(1), y(1), 'C1', 'FontSize', FNT_SZ, 'Color', 'b');
    hold on;
    h=ezplot(C1_eq, [0, size_x, 0, size_y]);
    set(h, 'LineWidth', 2, 'Color', 'b');
    hold on;


% Choose 5 points for C2.
    title("Choose 5 Points For C2, Press Enter When Done");
    [x, y]=getpts;

    A=[x.^2 x.*y y.^2 x y ones(size(x))];
    N = null(A);
    cc = N(:, 1);
    [C2_a, C2_b, C2_c, C2_d, C2_e, C2_f] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));
    C2=[C2_a C2_b/2 C2_d/2; C2_b/2 C2_c C2_e/2; C2_d/2 C2_e/2 C2_f]

    %Draw C2
    C2_eq = @(x_var, y_var) C2_a*x_var.^2 + C2_b*x_var.*y_var + C2_c*y_var.^2 + C2_d*x_var + C2_e*y_var + C2_f;
    text(x(1), y(1), 'C2', 'FontSize', FNT_SZ, 'Color', 'b');
    hold on;
    h=ezplot(C2_eq, [0, size_x, 0, size_y]);
    set(h, 'LineWidth', 2, 'Color', 'b');

% Print l1, l2, C1, C2
    l1 = l1./l1(3)
    l2 = l2./l2(3)
    C1 = C1./norm(C1)
    C2 = C2./norm(C2)

% Save Feature Selected Image
title("Image with Features:");
saveas(gcf, 'PalazzoTe_FeatureSelected.png');








