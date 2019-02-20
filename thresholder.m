clear all
im = imread('IMG_3761.JPG');
im = imresize(im, 0.2);
mask = createMask(im);

se = strel('sphere', 1);
mask = imerode(mask,se);
mask = imdilate(mask,se);
masked = im.*repmat(uint8(mask), [1 1 3]);

edges = detectEdges(masked);
figure(1)
imshowpair(masked,edges,'montage')

lines = detectLines(edges);
figure, imshow(edges), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

function [edges] = detectEdges(im)
    I = rgb2gray(im);
    edges = edge(I,'canny',0.02);
end
function [lines] = detectLines(BW)
    [H,T,R] = hough(BW);
    P  = houghpeaks(H,150,'threshold',ceil(0.3*max(H(:))));
    lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);
end
function [mask] = createMask(RGB)
I = rgb2lab(RGB);

channel1Min = 0.000;
channel1Max = 100.000;
channel2Min = -11.186;
channel2Max = 11.861;
channel3Min = -24.879;
channel3Max = 14.744;
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);
mask = ~sliderBW;
end
