clc
clear all
close all

% Rotate the image by -90 degrees using bilinear interpolation
originalImage = imrotate(imread('PalazzoTe.jpg'), -90, 'bilinear');

size_img=size(originalImage);

% Plotting Variables
    size_y=size_img(1);
    size_x=size_img(2);

% Convert the rotated image to grayscale
grayImage = rgb2gray(originalImage);

% Perform edge detection using the Canny edge detector
edges = edge(grayImage, 'Canny');

% Display the rotated image with filtered edges overlaid
imshow(originalImage);
hold on;

% Plot edges with a conical or elliptical appearance and minimum length
[B, L] = bwboundaries(edges, 'noholes');

% Use regionprops on the pixel coordinates
stats = regionprops(edges, 'MajorAxisLength');

% Initialize cell arrays to store selected boundaries and conic coefficients
selectedBoundaries = cell(0, 1);
selectedPointsCell = cell(0, 1);


for k = 1:length(B)
    boundary = B{k};
    
    % Get ellipse properties
    majorAxisLength = stats(k).MajorAxisLength;
    
    if majorAxisLength > 2000
        % Store the boundary in the cell array
        selectedBoundaries = [selectedBoundaries; {boundary}];
        
        % Choose 5 equally spaced x values from within the boundary,
        % excluding the first and last 5%
        x = boundary(:, 2);
        y = boundary(:, 1);
        xRange = max(x) - min(x);
        excludedRange = 0.05 * xRange;
        xExcludeStart = min(x) + excludedRange;
        xExcludeEnd = max(x) - excludedRange;
        
        xExcluded = x(x >= xExcludeStart & x <= xExcludeEnd);
        yExcluded = y(x >= xExcludeStart & x <= xExcludeEnd);
        
        [xUnique, idxUnique] = unique(xExcluded);
        yUnique = yExcluded(idxUnique);
        
        selectedX = linspace(min(xUnique), max(xUnique), 5);
        selectedY = interp1(xUnique, yUnique, selectedX);

        % Save selectedX and selectedY in the cell array
        selectedPointsCell = [selectedPointsCell; {[selectedX; selectedY]}];
        
        % Plot selected points
        plot(selectedX, selectedY, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
        hold on;
    end
end

C1_selected_x=selectedPointsCell{1}(1,:)';
C1_selected_y=selectedPointsCell{1}(2,:)';
C2_selected_x=selectedPointsCell{2}(1,:)';
C2_selected_y=selectedPointsCell{2}(2,:)';




A=[C1_selected_x.^2 C1_selected_x.*C1_selected_y C1_selected_y.^2 C1_selected_x C1_selected_y ones(size(C1_selected_x))];
N = null(A);
cc = N(:, 1);
[C1_a, C1_b, C1_c, C1_d, C1_e, C1_f] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));

A=[C2_selected_x.^2 C2_selected_x.*C2_selected_y C2_selected_y.^2 C2_selected_x C2_selected_y ones(size(C2_selected_x))];
N = null(A);
cc = N(:, 1);
[C2_a, C2_b, C2_c, C2_d, C2_e, C2_f] = deal(cc(1),cc(2),cc(3),cc(4),cc(5),cc(6));

% Plotting C1
C1_eq = @(x_var, y_var) C1_a*x_var.^2 + C1_b*x_var.*y_var + C1_c*y_var.^2 + C1_d*x_var + C1_e*y_var + C1_f;
h_c1 = ezplot(C1_eq, [0, size_x, 0, size_y]);
set(h_c1, 'LineWidth', 2, 'Color', 'b');
hold on;

% Plotting C2
C2_eq = @(x_var, y_var) C2_a*x_var.^2 + C2_b*x_var.*y_var + C2_c*y_var.^2 + C2_d*x_var + C2_e*y_var + C2_f;
h_c2 = ezplot(C2_eq, [0, size_x, 0, size_y]);
set(h_c2, 'LineWidth', 2, 'Color', 'b');
hold on;


% Draw the boundary of the conic-looking edges
for i=1:length(selectedBoundaries)
    boundary=selectedBoundaries{i};
    plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2);
    hold on;
end

% Save Feature Selected Image
title("Image with Extracted Conics:");
saveas(gcf, 'PalazzoTe_FeatureExtracted_Conic.png');



% Plot them together.
figure(2);
imshow(originalImage), hold on;
load('Variables_0B_Line.mat');
plot(1:size(originalImage, 2), y_line1_fit, 'b', 'LineWidth', 2), hold on;  %l1
plot(B_line1(:, 2), B_line1(:, 1), 'r', 'LineWidth', 2), hold on;
plot(1:size(originalImage, 2), y_line2_fit, 'b', 'LineWidth', 2), hold on;  %l2
plot(B_line2(:, 2), B_line2(:, 1), 'r', 'LineWidth', 2), hold on;
h_c1 = ezplot(C1_eq, [0, size_x, 0, size_y]);                               %C1
set(h_c1, 'LineWidth', 2, 'Color', 'b'), hold on;
h_c2 = ezplot(C2_eq, [0, size_x, 0, size_y]);                               %C2
set(h_c2, 'LineWidth', 2, 'Color', 'b'), hold on;

% Plot selected points
plot(C1_selected_x, C1_selected_y, 'bo', 'MarkerSize', 8, 'LineWidth', 2), hold on;
plot(C2_selected_x, C2_selected_y, 'bo', 'MarkerSize', 8, 'LineWidth', 2), hold on;

% Draw the boundary of the conic-looking edges
for i=1:length(selectedBoundaries)
    boundary=selectedBoundaries{i};
    plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2);
    hold on;
end

title("Image with Extracted Lines and Conics:");
saveas(gcf, 'PalazzoTe_FeatureExtracted.png');

