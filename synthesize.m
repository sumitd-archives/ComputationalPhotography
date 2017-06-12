function [Image, time] = synthesize(filename, winsize, newRows, newCols)
tic
%MaxErrThreshold = 0.1;

img = im2double(imread(filename)); 

inputImg =  img;  

[rows, cols, channels] = size(inputImg); 
windowlessSize = [(rows - winsize + 1) (cols - winsize + 1)];

halfWindow = (winsize - 1) / 2;

npixels = newRows * newCols; 
Image = zeros(newRows, newCols, 1); 

patches = im2col(inputImg(:, :, 1), [winsize winsize], 'sliding'); 

%initialize new texture with a random 3x3 patch from the sample
randRow = ceil(rand() * (rows - 2)); 
randCol = ceil(rand() * (cols - 2));

seedSize = 3; 
seedRows = ceil(newRows/2):ceil(newRows/2)+seedSize-1;
seedCols = ceil(newCols/2):ceil(newCols/2)+seedSize-1;
Image(seedRows, seedCols, :) = inputImg(randRow:randRow+seedSize-1, randCol:randCol+seedSize-1, :);

nfilled = seedSize * seedSize; 
filled = repmat(false, [newRows newCols]); 
filled(seedRows, seedCols) = repmat(true, [3 3]); 

gaussMask = fspecial('gaussian',winsize, winsize/6.4);

nskipped = 0; 

while nfilled < npixels    
    progress = false;
    
    [pixelRows, pixelCols] = GetUnfilledPixels(filled, winsize);
     
     for i = 1:length(pixelRows)
        pixelRow = pixelRows(i);
        pixelCol = pixelCols(i);
        
        rowRange = pixelRow-halfWindow:pixelRow+halfWindow;
        colRange =  pixelCol - halfWindow:pixelCol + halfWindow;

        deadRows = rowRange < 1 | rowRange > newRows;
        deadCols = colRange < 1 | colRange > newCols; 


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

       %[bestMatches, SSD] = FindMatches(template, validMask, gaussMask, red_patches, green_patches, blue_patches);
       [bestMatches, SSD] = FindNearestPatch(template, validMask, gaussMask, patches);


        matchIdx = RandomPick(bestMatches);
        matchError = SSD(matchIdx); 

         %if matchError < MaxErrThreshold 
             [matchRow, matchCol] = ind2sub(windowlessSize, matchIdx); 
             
             %match coords are at corner of window and need to be offset
             matchRow = matchRow + halfWindow;
             matchCol = matchCol + halfWindow;  

             Image(pixelRow, pixelCol, :) = inputImg(matchRow, matchCol, :);

             filled(pixelRow, pixelCol) = true;   
             nfilled = nfilled + 1; 
             progress = true;
         %else
         %    nskipped = nskipped + 1; 
         %end
    end
    
    
    disp(sprintf('Pixels filled: %d / %d', nfilled, npixels)); 
          
  %      MaxErrThreshold = MaxErrThreshold * 1.1;
 end

 %figure;
 %imshow(Image);
 toc
time = toc; 

function [pixelRows, pixelCols] = GetUnfilledPixels(filled, winsize) 
    border = bwmorph(filled,'dilate')-filled;
    
    [Rows, Cols] = find(border);
    len = length(Rows); 
     
     filledSums = colfilt(filled, [winsize winsize], 'sliding', @sum); 
     filledNeighborsCnt = filledSums( sub2ind(size(filled), Rows, Cols) ); 
     [sorted, sorted_Index] = sort(filledNeighborsCnt, 1, 'descend');
     
     pixelRows = Rows(sorted_Index); 
     pixelCols = Cols(sorted_Index); 

     
function idx = RandomPick(matches)
    indices = find(matches);
    idx = indices(ceil(rand() * length(indices)));
        
     
function [pixelList, SSD] = FindNearestPatch (template, validMask, gaussMask, patches)
ErrThreshold = 0.3; 

[pixels_per_patch, npatches] = size(patches); 

totalWeight = sum(sum(gaussMask(validMask)));

mask = (gaussMask .* validMask) / totalWeight;
mask_row = mask(:)'; 
 
new_patch = reshape(template, [pixels_per_patch 1]); 
new_patch_copies = repmat(new_patch, [1 npatches]);

dist =  mask_row * (new_patch_copies - patches).^2; 

SSD = dist;
pixelList = SSD <= min(SSD) * (1+ErrThreshold);
