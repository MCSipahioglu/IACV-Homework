clc;
clear all;
close all;

% Rotate the image by -90 degrees using bilinear interpolation
originalImage = imrotate(imread('PalazzoTe.jpg'), -90, 'bilinear');

% Plotting Variables
    size_img=size(originalImage);
    size_y=size_img(1);
    size_x=size_img(2);

% Convert the rotated image to grayscale
grayImage = rgb2gray(originalImage);

% Perform edge detection using the Canny edge detector
edges = edge(grayImage, 'Canny');
[B, L] = bwboundaries(edges, 'noholes');

% Display the rotated image with filtered edges overlaid
imshow(originalImage);
hold on;

% Use regionprops on the pixel coordinates
stats = regionprops(edges, 'Eccentricity', 'Orientation', 'MajorAxisLength');

% Initialize variables to store the longest edges
longestEdge1 = struct('Index', 0, 'Length', 0);
longestEdge2 = struct('Index', 0, 'Length', 0);

% Find the longest edges with specified orientations and eccentricity
for k = 1:length(B)
    boundary = B{k};
    
    % Get eccentricity, orientation, and major axis length
    eccentricity = stats(k).Eccentricity;
    orientation = stats(k).Orientation;
    majorAxisLength = stats(k).MajorAxisLength;

    % Check eccentricity criterion
    if eccentricity > 0.995
        % Check orientation criterion for the first edge
        if orientation > -85 && orientation > 80 && majorAxisLength > longestEdge1.Length
            longestEdge1.Index = k;
            longestEdge1.Length = majorAxisLength;
        end

        % Check orientation criterion for the second edge
        if orientation > 45 && orientation < 60 && majorAxisLength > longestEdge2.Length
            longestEdge2.Index = k;
            longestEdge2.Length = majorAxisLength;
        end
    end
end

% Plot the longest edges and draw the fitted lines
if longestEdge1.Index > 0
    
    % Get coordinates of the longest edge
    x1 = B{longestEdge1.Index}(:, 2);
    y1 = B{longestEdge1.Index}(:, 1);
    % Fit a line to the coordinates
    lineParams1 = polyfit(x1, y1, 1);
    a1 = -lineParams1(1);
    b1 = 1;
    c1 = -lineParams1(2);
    % Draw the fitted line across the image width
    y_line1_fit = (-a1 * (1:size(originalImage, 2)) - c1) / b1;
    plot(1:size(originalImage, 2), y_line1_fit, 'b', 'LineWidth', 2);
    plot(B{longestEdge1.Index}(:, 2), B{longestEdge1.Index}(:, 1), 'r', 'LineWidth', 2);
    B_line1=B{longestEdge1.Index};
end

if longestEdge2.Index > 0
    % Get coordinates of the longest edge
    x2 = B{longestEdge2.Index}(:, 2);
    y2 = B{longestEdge2.Index}(:, 1);
    % Fit a line to the coordinates
    lineParams2 = polyfit(x2, y2, 1);
    a2 = -lineParams2(1);
    b2 = 1;
    c2 = -lineParams2(2);
    % Draw the fitted line across the image width
    y_line2_fit = (-a2 * (1:size(originalImage, 2)) - c2) / b2;
    plot(1:size(originalImage, 2), y_line2_fit, 'b', 'LineWidth', 2);
    plot(B{longestEdge2.Index}(:, 2), B{longestEdge2.Index}(:, 1), 'r', 'LineWidth', 2);
    B_line2=B{longestEdge2.Index};
end

hold off;



% Column vector for the first line equation l1 = [a1; b1; c1]
l1 = [-lineParams1(1); 1; -lineParams1(2)];
l1=l1./l1(3)

% Column vector for the second line equation l2 = [a2; b2; c2]
l2 = [-lineParams2(1); 1; -lineParams2(2)];
l2=l2./l2(3)

% Save Feature Selected Image
title("Image with Extracted Lines:");
saveas(gcf, 'PalazzoTe_FeatureExtracted_Line.png');

save('Variables_0B_Line.mat','y_line1_fit','y_line2_fit','B_line1','B_line2');
