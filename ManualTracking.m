clear;
clc;
% Add required functions to the Matlab path.
addpath(genpath(cd));

[fname, path] = uigetfile("*.mp4", "Select videos to track", "MultiSelect", "on");

if ischar(fname)
    fname = {fname};
end

frameImg = [];
hFig = figure('Name', 'Frame Viewer', 'NumberTitle', 'off');
% Initialize FasData as a struct with empty arrays
FasData.FLength = [];
FasData.FAngle = [];
FasData.pts = {};
FasData.digitizedFrames = [];

ApoData.Angle = [];
ApoData.pts = {};
ApoData.digitizedFrames = [];


for n_file = 1:numel(fname)
    vidObj = VideoReader([path fname{n_file}]);
    %n_fr = 1;

    step = (vidObj.NumFrames/100);
    frames_to_track = 1: step : vidObj.NumFrames;
    % Create a random array of integers for the frames to track manually
    rnd_fr = randperm(length(frames_to_track)); %digitez n_frames / 10

    ZeroPadL = 200 * ones(vidObj.Height, vidObj.Width, 'uint8');
    ZeroPadR = 200 * ones(vidObj.Height, vidObj.Width, 'uint8');

    for ii = 1:numel(rnd_fr)
        % Read the frame
        frameImg = rgb2gray(read(vidObj, rnd_fr(ii)));
        currentImage = [ZeroPadL, frameImg, ZeroPadR];

        imshow(currentImage);
        hold on;
        fas = drawline('Color', 'red');
        apo = drawline('Color', 'green');

        % Store fas and FasData in guidata
        data.fas = fas;
        data.FasData = FasData;
        data.apo = apo;
        data.ApoData = ApoData;

        data.currentFrame = rnd_fr(ii); % Store the current frame


        guidata(hFig, data);

        % Create "Accept Line" button
        hButton = uicontrol('Style', 'pushbutton', 'String', 'Accept Line', ...
            'Position', [20 20 100 40], 'Callback', @acceptLineCallback);

        % Wait for button press
        uiwait(gcf);
        delete(hButton); % Remove the button after it's pressed

        % Update FasData from guidata
        data = guidata(hFig);
        FasData = data.FasData;
        ApoData = data.ApoData;
    end
    FasData.FileName = fname{n_file};
end

close(hFig); % Close the figure when done
uisave({'FasData','ApoData'},'Manual_Tracking.mat');

%% Function to handle button click
function acceptLineCallback(~, ~)
% Retrieve fas and FasData from guidata
data = guidata(gcbf);
fas = data.fas;
FasData = data.FasData;
ApoData = data.ApoData;
% Get the position of the line
pos = fas.Position;
% Calculate the length of the line
length = sqrt((pos(2,1) - pos(1,1))^2 + (pos(2,2) - pos(1,2))^2);
% Calculate the angle of the line (in degrees)
angle = atan2d((pos(2,2) - pos(1,2)), (pos(2,1) - pos(1,1)));

% Append the data to the FasData struct
FasData.FLength(end+1) = length;
FasData.FAngle(end+1) = angle;
FasData.pts{end+1} = pos;
FasData.digitizedFrames(end+1) = data.currentFrame;

% Apo calculation
apo = data.apo;
%apoData = data.apoData;

% Get the position of the line
pos = apo.Position;
% Calculate the length of the line
%length = sqrt((pos(2,1) - pos(1,1))^2 + (pos(2,2) - pos(1,2))^2);
% Calculate the angle of the line (in degrees)
angle = atan2d((pos(2,2) - pos(1,2)), (pos(2,1) - pos(1,1)));

% Append the data to the FasData struct
%ApoData.FLength(end+1) = length;
ApoData.Angle(end+1) = angle;
ApoData.pts{end+1} = pos;
ApoData.digitizedFrames(end+1) = data.currentFrame;

% Store updated FasData in guidata
data.FasData = FasData;
data.ApoData = ApoData;
guidata(gcbf, data);

% Resume UI execution
uiresume(gcbf);
end

