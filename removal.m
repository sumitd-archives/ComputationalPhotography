function [] = removal(filename , winsize)
tic
gray_img = rgb2gray(imread(filename));
target_object = roipoly(gray_img);
idx = find(target_object);
gray_img(idx) = 0;
%incompleteImg = im2double(imread(filename));
incompleteImg = im2double(gray_img);
Image = incompleteImg;
[rows, cols] = size(Image); 
windowlessSize = [(rows - winsize + 1) (cols - winsize + 1)];
halfWindow = (winsize-1)/2;
filled = repmat(1 , size(incompleteImg));
unfilledIndices = find(incompleteImg == 0);
nUnfilled = size(unfilledIndices);
nUnfilled = nUnfilled(1);
filled(unfilledIndices) = 0;
patches = im2col(incompleteImg(:, :, 1), [winsize winsize], 'sliding');
[patch_size patch_count] = size(patches);

%remove patches with zero pixel at center
columns = 1:patch_count;
row = repmat((winsize^2 + 1) / 2, [1 patch_count]);
idx = sub2ind(size(patches) , row , columns);
zero_idx = find(patches(idx) == 0);
patches(:,zero_idx) = [];
columns(:,zero_idx) = [];

gaussMask = fspecial('gaussian',winsize, winsize/6.4);

nfilled = 0;
while nfilled < nUnfilled    
    %progress = false;
    
    [pixelRows, pixelCols] = GetUnfilledPixels(filled , winsize);
     
     for i = 1:length(pixelRows)
        disp(sprintf('Pixels filled: %d / %d\n', nfilled, nUnfilled)); 
        pixelRow = pixelRows(i);
        pixelCol = pixelCols(i);
   
        rowRange = pixelRow-halfWindow:pixelRow+halfWindow;
        colRange = pixelCol - halfWindow:pixelCol + halfWindow;

        deadRows = rowRange < 1 | rowRange > rows;
        deadCols = colRange < 1 | colRange > cols; 

        if sum(deadRows) + sum(deadCols) > 0 
            safeRows = rowRange(~deadRows); 
            safeCols = colRange(~deadCols); 

            template = zeros(winsize, winsize); 
            template(~deadRows, ~deadCols, :) = Image(safeRows, safeCols, :); 

            validMask = repmat(false, [winsize winsize]); 
            validMask(~deadRows, ~deadCols) = filled(safeRows, safeCols); 
        else
            template = Image(rowRange, colRange, :);
            validMask = filled(rowRange, colRange); 

        end

        [bestMatches, SSD] = FindMatches(template, validMask, gaussMask, patches);

        matchIdx = RandomPick(bestMatches);
        if (matchIdx == 0)
            continue;
        end
        
        [matchRow, matchCol] = ind2sub(windowlessSize, columns(matchIdx)); 
             
             %match coords are at corner of window and need to be offset
             matchRow = matchRow + halfWindow;
             matchCol = matchCol + halfWindow;  

             Image(pixelRow, pixelCol, :) = incompleteImg(matchRow, matchCol, :);

             filled(pixelRow, pixelCol) = 1;   
             nfilled = nfilled + 1; 
             progress = true;
     end
end
toc
time = toc
figure;
imshow(Image);

function [pixelRows, pixelCols] = GetUnfilledPixels(filled, winsize) 
    size(filled) 
    border = bwmorph(filled,'dilate')-filled;
    
    [pixelRows, pixelCols] = find(border);
    len = length(pixelRows); 
     
     %sort by number of neighbors     
     filledSums = colfilt(filled, [winsize winsize], 'sliding', @sum); 
     numFilledNeighbors = filledSums( sub2ind(size(filled), pixelRows, pixelCols) ); 
     [sorted, sortIndex] = sort(numFilledNeighbors, 1, 'descend');
     
     pixelRows = pixelRows(sortIndex); 
     pixelCols = pixelCols(sortIndex); 


%% Pick a random pixel from valid patches
function idx = RandomPick(matches)
    indices = find(matches);
    if (length(indices) == 0)
        idx = 0;
    else
        idx = indices(ceil(rand() * length(indices)));
    end
    %idx = ceil(rand() * length(in));
      
     
%% Find candidate patches that match template
function [pixelList, SSD] = FindMatches (template, validMask, gaussMask, patches)
%ErrThreshold = 0.3; 
[pixels_per_patch, npatches] = size(patches); 
indices = find(validMask);
totalWeight = sum(sum(gaussMask(indices)));
mask = (gaussMask .* validMask) / totalWeight;
mask_vec = mask(:)'; 
new_patch = reshape(template, [pixels_per_patch 1]); 
new_patch_copies = repmat(new_patch, [1 npatches]);
dist =  mask_vec * (new_patch_copies - patches).^2; 
SSD = dist;
%pixelList = SSD <= min(SSD) * (1+ErrThreshold);
pixelList = SSD == min(SSD);














