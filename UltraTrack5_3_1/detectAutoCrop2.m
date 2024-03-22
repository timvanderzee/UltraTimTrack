function boundToCrop = detectAutoCrop2(frames,FrameRate)
    try
        % Determine the size of the frames matrix.
        [vidHeight, vidWidth, numberOfFrames] = size(frames);

        % Initialize a binary matrix.
        binT = logical(zeros(vidHeight, vidWidth));

        % Calculate the mean gray level for each frame.
        meanGrayLevels = squeeze(mean(mean(frames, 1), 2));

        % Initialize the adaptive background.
        alpha = 0.5;
        Background = frames(:, :, 1);

        for frame = 2 : round(FrameRate/4) :  numberOfFrames
            % Change background slightly at each frame.
            Background = (1 - alpha) * frames(:, :, frame) + alpha * Background;

            % Calculate the difference between this frame and the background.
            differenceImage = frames(:, :, frame) - uint8(Background);

            % Threshold with Otsu method.
            %grayImage = rgb2gray(differenceImage);
            grayImage = (differenceImage);
            thresholdLevel = graythresh(grayImage);
            binaryImage = im2bw(grayImage, thresholdLevel);

            % Add binary image to the sum.
            binT = binT + binaryImage;
        end

    catch ME
        % Handle errors.
        strErrorMessage = sprintf('Error!!!');
        disp(strErrorMessage);
        return;
    end

    % Apply a median filter, [10 10] neighbor pixels to the binary matrix.
    filteredMatrix = medfilt2(binT, [10, 10], 'zeros');
    %maybe check the filled area to be sure that  i don't small bastards
    %around
    
    % Find connected components in the filtered matrix.
    stats = regionprops(filteredMatrix, 'BoundingBox');

    % Extract the bounding box information.
    boundingBoxes = cat(1, stats.BoundingBox);

    % Calculate the overall bounding box that encompasses all smaller bounding boxes.
    boundToCrop = [min(boundingBoxes(:, 1)), min(boundingBoxes(:, 2)), ...
                   max(boundingBoxes(:, 1) + boundingBoxes(:, 3)) - min(boundingBoxes(:, 1)), ...
                   max(boundingBoxes(:, 2) + boundingBoxes(:, 4)) - min(boundingBoxes(:, 2))];
end
