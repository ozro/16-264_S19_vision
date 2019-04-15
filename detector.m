clear all
%% Params
patch_size = 20;
num_lines = 5;
line_dist_thresh = 3;
line_prune_thresh = 1;

%% Image preprocessing
im = imread('plank.png');
% im = imread('plank_blur.jpg');
%im = imgaussfilt(im, 1);

% im = imread('IMG_3761_cropped.jpg');
%im = imresize(im, 0.9);
% im = getMasked(im);

gray = rgb2gray(im);
%% Find corners
C = corner(gray);

%% Find edge candiates from patches
[im_padded, patches, edges, lines_patch, lines_world] = getPatchLines(gray,C,patch_size, num_lines);
eqs = getLineEquations(lines_world);
% plotCorners(2, im_padded, C + patch_size, patch_size, 'c');
% plotLines(1, im_padded, eqs, 'r:')
% plotPatches(1, C, patches, edges, lines_patch, num_lines)

%% Find matched lines
eqs_dist = matchLines(eqs, C + patch_size, im_padded, line_dist_thresh);
% plotLines(1, im_padded, eqs_dist(:, 6:7), 'r:')
plotCorners(1, im_padded, C+patch_size, 0, 'r');
[C_pruned] = pruneVertices(eqs_dist, C + patch_size, line_prune_thresh, im_padded);
plotCorners(1, [], C_pruned(:, 1:2), 2, 'c');

[eqs_matched, matches] = refineMatches(eqs_dist, C_pruned);

plotLines(1, [], eqs_matched(:,6:7), 'g--')
plotMatchedLines(2, im_padded, eqs_matched);

%% Find correct vertices
[vertices, eqs_valid] = getValidVertices(matches, eqs_matched);
plotValidVertices(3, im_padded, vertices, eqs_valid);

%% Plotting
function plotCorners(fig, im_padded, C, patch_size, color)
    figure(fig)
    if(size(im_padded,1) > 0)
        imshow(im_padded);
    end
    hold on
    for i = 1:size(C,1)
        x = C(i,1);
        y = C(i,2);
        plot(x,y, '*', 'Color', color)
        text(x+patch_size-5,y-patch_size-5,"  " + i, 'Color',color)
        rectangle('Position',[x-patch_size, y-patch_size, 2*patch_size+1, 2*patch_size+1], 'EdgeColor', color)
    end
    hold off
end
function plotLines(fig, im_padded, eqs, style)
    figure(fig)
    if(size(im_padded,1) > 0)
        imshow(im_padded);
    end
    pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
    width = pos(3);
    hold on
    for i = 1:size(eqs,3)
        for k = 1:size(eqs, 1)
            eq = eqs(k, :, i);
            [x, y] = getEqPlotPoints(eq, width);
            plot(x,y,style,'LineWidth',.25);
        end
    end
    hold off
end
function plotPatches(fig, C, patches, edges, lines_patch, num_lines)
    figure(fig)
    for i = 1:size(C, 1)
        subplot(3,size(C,1),i)
        imshow(patches(:,:,i)/255)
        subplot(3,size(C,1),i+size(C,1))
        imshow(edges(:,:,i))
        subplot(3,size(C,1),i+size(C,1)*2)
        imshow(edges(:,:,i))
        hold on
        for k = 1:num_lines
           xy = lines_patch(k*2-1:k*2,:,i);
           plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
           plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
           plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
        end
        hold off
    end
end
function plotMatchedLines(fig, im_padded, eqs_matched)
    figure(fig)
    if(size(im_padded,1) > 0)
        imshow(im_padded);
    end
    pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
    width = pos(3);
    hold on
    for k = 1:size(eqs_matched, 1)
        eq = eqs_matched(k, 6:7);
        [x, y] = getEqPlotPoints(eq, width);
        plot(x,y,'g--', 'LineWidth',.5);
        c0 = eqs_matched(k, 8:9);
        c1 = eqs_matched(k, 10:11);
        plot([c0(1) c1(1)], [c0(2) c1(2)], 'y*-', 'LineWidth', 1.5);
    end
    hold off
end
function plotValidVertices(fig, im_padded, vertices, eqs_valid)
    figure(fig)
    if(size(im_padded,1) > 0)
        imshow(im_padded);
    end
    pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
    width = pos(3);
    hold on
    for i = 1:size(eqs_valid, 1)
        eq = eqs_valid(i, 6:7);
        [x, y] = getEqPlotPoints(eq, width);
        plot(x,y,'g--', 'LineWidth',.5);
        c0 = eqs_valid(i, 8:9);
        c1 = eqs_valid(i, 10:11);
        plot([c0(1) c1(1)], [c0(2) c1(2)], 'y*-', 'LineWidth', 1.5);
    end
    plot(vertices(:,1), vertices(:,2), 'co', 'MarkerFaceColor', 'c')
    hold off
end
function [x, y] = getEqPlotPoints(eq, width)
    m = eq(1);
    b = eq(2);
    x0 = 0;
    y0 = m*x0 + b;
    x1 = width;
    y1 = m*x1 + b;
    x = [x0 x1];
    y = [y0 y1];
end

%% Functions
function [vertices, eqs_valid] = getValidVertices(matches, eqs_matched)
    valid = find(sum(matches, 2) == 4);
    vertices = [];
    eqs_valid = zeros(0, 11);
    for i = 1:size(valid)
        v = valid(i);
        mask = (eqs_matched(:, 1) == v) | (eqs_matched(:, 2) == v);
        eqs = eqs_matched(mask,6:7);
        [xi, yi] = getBestIntersection(eqs);
        vertices = [vertices; [xi, yi]];
        eqs_valid = [eqs_valid; eqs_matched(mask, :)];
    end
end
function [eqs_matched, matches] = refineMatches(eqs_dist, C_pruned)
    matches = eye(size(C_pruned, 1));
    eqs_matched = zeros(0, size(eqs_dist, 2));
    % Loop over each possible vertex pairing
    for i = 1:size(C_pruned, 1) 
        for j = i+1:size(C_pruned,1)
            % Get eqs that are connecting the vertices
            rows = (eqs_dist(:,1) == C_pruned(i, 3) & (eqs_dist(:,2) == C_pruned(j, 3)));
            eqs_i = eqs_dist(rows, :);
            if(size(eqs_i, 1) > 0)
                matches(i, j) = matches(i,j) + 1;
                matches(j, i) = matches(j,i) + 1;
                for m = 1:size(eqs_i)
                    eq_avg = eqs_i(m,:);
                    if(sum(~isfinite(eq_avg)) == 0)
                        break
                    end
                end
                if(sum(~isfinite(eq_avg)) ~= 0)
                    break
                end
%                 eq_avg = sum(eqs_i, 1)/size(eqs_i, 1);
                eq_avg(1) = C_pruned(i,3);
                eq_avg(2) = C_pruned(j,3);
                eqs_matched = [eqs_matched; eq_avg];
            end
        end
    end
    eqs_matched = sortrows(eqs_matched);
end
function [C_pruned] = pruneVertices(eqs, C, line_prune_thresh, im_padded)
    invalid = false(size(C, 1), 1);
    % Loop over each possible vertex pairing
    for i = 1:size(C, 1) 
        if(invalid(i))
            continue
        end
        for j = i+1:size(C,1)
            if(invalid(j))
                continue;
            end
            % Get eqs that are connecting the vertices
            rows = (eqs(:,1) == i) & (eqs(:,2) == j);
            eqs_i = eqs(rows, :);
            if(size(eqs_i, 1) > 0)
                %eq_avg = sum(eqs_i, 1)/size(eqs_i, 1);
                for m = 1:size(eqs_i)
                    eq_avg = eqs_i(m,:);
                    if(sum(~isfinite(eq_avg)) == 0)
                        break
                    end
                end
                if(sum(~isfinite(eq_avg)) ~= 0)
                    break
                end
                % Loop over every vertex that aren't the connected ones
                valid_corners = [];
                invalid_corners = [];
                for v = 1:size(C, 1)
                    if(v ~= i && v ~= j && ~invalid(v))
                        m = eq_avg(6);
                        b = eq_avg(7);
                        x = C(v,1);
                        y = C(v,2);
                        dist = getDistToLine(m,b, x, y);    
                        % Check if vertex is on the line
                        if(dist < line_prune_thresh)
                            [xp, yp] = projectPoint(m, b, x, y);
                            x0 = eq_avg(8);
                            y0 = eq_avg(9);
                            x1 = eq_avg(10);
                            y1 = eq_avg(11);
                            dist0 = sqrt((x0-xp)^2 + (y0-yp)^2);
                            dist1 = sqrt((x1-xp)^2 + (y1-yp)^2);
                            distp = sqrt((x1-x0)^2 + (y1-y0)^2);
                            
                            % Check if vertex is between others
                            if abs(((dist0 + dist1)-distp)) < 0.01
                                invalid(v) = 1;
                                invalid_corners = [invalid_corners;C(v,:)];
                            else
                                valid_corners = [valid_corners;C(v,:)];
                            end
                        else
                            valid_corners = [valid_corners;C(v,:)];
                        end
                    end
                end
%                 C_avg = [eq_avg(8:9);eq_avg(10:11)];
%                 plotCorners(2, im_padded, C_avg, 0, 'c');
%                 plotLines(2, [], eq_avg(:,6:7), 'g')
%                 plotLines(2, [], eqs_i(:,6:7), 'b:')
%                 plotCorners(2, [], invalid_corners, 20, 'r');
%                 plotCorners(2, [], valid_corners, 0, 'g');
%                 pause(0.1);
            end
        end
    end
    C_pruned = [C(~invalid, :), find(~invalid)];
end
function [eqs_dist] = matchLines(eqs, C, im_padded, thresh)
    eqs_dist = zeros(0, 11);
    for i=1:size(eqs,3)
        for k=1:size(eqs,1)
            eq = eqs(k, :, i);
            for j=1:size(C,1)
                if(i == j)
                    continue
                end
                if(sum(~isfinite(eq)) > 0)
                    continue
                end
                [dist0, dist1, length, p0, p1] = getLineProjections(eq, C(i,:), C(j,:), im_padded);
                if(dist0 > thresh || dist1> thresh)
                    continue
                end

                if(i < j)
                    eqs_dist = [eqs_dist; [i, j, dist1, dist0, length, eq, p0, p1]];
                else
                    eqs_dist = [eqs_dist; [j, i, dist1, dist0, length, eq, p1, p0]];
                end
            end
        end
    end
    eqs_dist = sortrows(eqs_dist);
end
function [im_padded, patches, edges, lines_patch, lines_world] = getPatchLines(im,C,patch_size, num_lines)
    im_padded = padarray(im, [patch_size patch_size], 0, 'both');
    patches = zeros(patch_size*2+1, patch_size*2+1, size(C,1));
    edges = zeros(patch_size*2+1, patch_size*2+1, size(C,1));
    lines_patch = zeros(num_lines*2, 2, size(C,1));
    lines_world = zeros(num_lines*2, 2, size(C,1));
    
    for i=1:size(C,1)
        x = C(i,1) + patch_size;
        y = C(i,2) + patch_size;
        row = y;
        col = x;
        patches(:,:,i) = im_padded(row-patch_size:row+patch_size, col-patch_size:col+patch_size);
        [xy, edge_patch] = getAllLines(patches(:,:,i), num_lines);
        edges(:,:,i) = edge_patch;
        lines_patch(:,:,i) = xy;
        lines_world(:,:,i) = xy + repmat([C(i,1) C(i,2)], [num_lines*2 1]);
    end
end
function [xy, edge_patch] = getAllLines(patch, num_lines)
    edge_patch = edge(patch/255, 'Canny', 0.1);
    [H,theta,rho] = hough(edge_patch);
    P = houghpeaks(H,num_lines,'threshold',ceil(0.3*max(H(:))));
    lines = houghlines(edge_patch,theta,rho,P,'FillGap',5,'MinLength',7);
    xy = zeros(num_lines*2, 2);
    for k = 1:min([num_lines, size(lines, 2)])
        xy(k*2-1:k*2, :) = [lines(k).point1; lines(k).point2];
    end
end
function [masked] = getMasked(im)
    mask = blueMask(im) | greenMask(im) | redMask(im);
    masked = im;
    mask = repmat(mask, [1 1 3]);
    masked(~mask) = 0;
end
%% Helpers
function [BW,maskedRGBImage] = blueMask(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 25-Mar-2019
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.555;
channel1Max = 0.771;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.313;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.000;
channel3Max = 1.000;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end
function [BW,maskedRGBImage] = greenMask(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 25-Mar-2019
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.181;
channel1Max = 0.338;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.313;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.000;
channel3Max = 1.000;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end
function [BW,maskedRGBImage] = redMask(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 25-Mar-2019
%------------------------------------------------------


% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.930;
channel1Max = 0.018;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.313;
channel2Max = 1.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.000;
channel3Max = 1.000;

% Create mask based on chosen histogram thresholds
sliderBW = ( (I(:,:,1) >= channel1Min) | (I(:,:,1) <= channel1Max) ) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
BW = sliderBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end
function [xi, yi] = getBestIntersection(eqs)
    points = zeros(0, 2);
    for i = 1:size(eqs, 1)
        for j = i+1:size(eqs,1)
            [x, y] = getIntersection(eqs(i, :), eqs(j,:));
            points = [points; [x,y]];
        end
    end
    avg = sum(points)/size(points,1);
    xi = avg(1);
    yi = avg(2);
end
function [xi, yi] = getIntersection(eq1, eq2)
    m1 = eq1(1);
    m2 = eq2(1);
    b1 = eq1(2);
    b2 = eq2(2);
    
    xi = (b2 - b1)/(m1-m2);
    yi = m1 * xi + b1;
end
function [dist] = getDistToLine(m, b, x, y)
    dist = abs(-m*x + y + -b)/(sqrt(m^2 + 1));
end
function [xp, yp] = projectPoint(m, b, x, y)
    mp = -1/m;
    bp = -mp * x + y;
    xp = (bp - b) / (m - mp);
    yp = mp * xp + bp;
end
function [dist0, dist1, length, p0, p1] = getLineProjections(eq, point0, point1, im_padded)
    m = eq(1);
    b = eq(2);
    x0 = point0(1);
    y0 = point0(2);
    x1 = point1(1);
    y1 = point1(2);
    
    dist0 = getDistToLine(m,b, x0, y0);
    dist1 = getDistToLine(m,b, x1, y1);
    
    [xp0, yp0] = projectPoint(m, b, x0, y0);
    [xp1, yp1] = projectPoint(m, b, x1, y1);
    
    length = sqrt((xp1-xp0)^2 + (yp1-yp0)^2);
    
    p0 = [xp0 yp0];
    p1 = [xp1 yp1];
    if(0)%dist0 < 2 && dist1 < 2)
        figure(2)
        imshow(im_padded)
        hold on;
        plot([0, 460], [b, m*460 + b])
        plot([x0, xp0], [y0, yp0], 'g')
        plot([x1, xp1], [y1, yp1], 'g')
        plot(x1, y1, 'r.')
        plot(x0, y0, 'g.')
        plot(xp0, yp0, 'yo')
        plot(xp1, yp1, 'yo')
        text((x0+xp0)/2, (y0+yp0)/2, ""+dist0, "Color", "g")
        text((x1+xp1)/2, (y1+yp1)/2, ""+dist1, "Color", "r")
        text((xp0+xp1)/2, (yp0+yp1)/2, ""+length, "Color", "c")
        hold off;
        xlim([0 460])
        ylim([0 460])
        pause
    end
end
function [eqs] = getLineEquations(lines)
    eqs = zeros(size(lines,1)/2, 2,size(lines,3));
    for i = 1:size(lines, 3)
        for k = 1:size(lines, 1)/2
            xy = lines(k*2-1:k*2, :, i);
            [m, b] = getCoefficients(xy);
            eqs(k, :, i) = [m, b];
        end
    end
end
function [m, b] = getCoefficients(xy)
    m = (xy(2,2)-xy(1,2))/(xy(2,1)-xy(1,1));
    b = xy(2,2) - m*xy(2,1);
end