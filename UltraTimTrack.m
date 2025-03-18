function varargout = UltraTimTrack(varargin)
% ULTRATIMTRACK M-file for UltraTimTrack.fig
%      ULTRATIMTRACK, by itself, creates a new ULTRATIMTRACK or raises the existing
%      singleton*.
%      H = ULTRATIMTRACK returns the handle to a new ULTRATIMTRACK or the handle to
%      the existing singleton*.
%
%      ULTRATIMTRACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ULTRATIMTRACK.M with the given input arguments.
%
%      ULTRATIMTRACK('Property','Value',...) creates a new ULTRATIMTRACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UltraTimTrack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UltraTimTrack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UltraTimTrack

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UltraTimTrack_OpeningFcn, ...
    'gui_OutputFcn',  @UltraTimTrack_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT

% --- Executes just before UltraTimTrack is made visible.
function UltraTimTrack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UltraTimTrack (see VARARGIN)

% Check if Parallel Computing Toolbox is installed
% chkParallelToolBox(); %if exists it runs infinitely till Ultratimtrack is closed

% Choose default command line output for UltraTimTrack
handles.output = hObject;

%add automatically all files and subfolders dynamically
filename = [mfilename,'.m'];
fullpath = which(filename);
mainfoldername = erase(fullpath,filename);
addpath(genpath(mainfoldername));

load('TimTrack_parms.mat','parms')
handles.parms = parms;

handles.BlockSize = [21 71]; %initialize it
% check to make sure there is a settings mat-file present, if not then make
% one in the directory where the m-file sits.
[program_directory, ~, ~] = fileparts(mfilename('fullpath'));

if ismac || isunix
    if isempty(dir([program_directory '/ultrasound_tracking_settings.mat']))
        ImageDepth = 50.7;
        % Needed for KL algorithm, but not KLT algorithm
        %Sigma = 3;
        %S_Step = 3;
        Position = get(gcf,'Position');
        Default_Directory = cd;
        save([program_directory '/ultrasound_tracking_settings.mat'], 'ImageDepth',...
            'Position', 'Default_Directory');
    else
        load ultrasound_tracking_settings.mat
    end
end
%
if ispc
    if isempty(dir([program_directory '\ultrasound_tracking_settings.mat']))
        ImageDepth = 50.7;
        % Needed for KL algorithm, but not KLT algorithm
        %Sigma = 3;
        %S_Step = 3;
        Position = get(gcf,'Position');
        Default_Directory = cd;
        save([program_directory '\ultrasound_tracking_settings.mat'], 'ImageDepth',...
            'Position', 'Default_Directory');
    else
        load ultrasound_tracking_settings.mat
    end

end

% load the settings mat-file and set default settings
handles.ID = ImageDepth;
set(handles.ImDepthEdit,'String',num2str(ImageDepth));
% Needed for KL algorithm, but not KLT algorithm
%handles.SIGMA = Sigma;
%handles.S_STEP = S_Step;
set(0,'RecursionLimit',3000)
set(gcf,'DoubleBuffer','on','Position',Position);
%cd(Default_Directory)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UltraTimTrack wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function menu_load_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% First clean up some variables from any previously loaded files
if isfield(handles,'movObj')
    %handles = rmfield(handles,'mov');
    handles = rmfield(handles,'movObj');
end
if isfield(handles,'BIm')
    handles=rmfield(handles,'BIm');
    handles=rmfield(handles,'Bheader');
end

if isfield(handles,'ImStack')
    handles=rmfield(handles,'ImStack');
end
if isfield(handles,'points')
        handles=rmfield(handles,'points');
end
if isfield(handles,'Region')

    handles = rmfield(handles,'Region');
end

% if ~isempty(get(handles.keyframe_list,'String'))
%     set(handles.keyframe_list,'String',[])
% end

if isfield(handles,'crop_rect')
    handles.crop_rect = [];
end


version = ver;
% do some checking to make sure the image processing toolbox is available
if isempty(strfind([version.Name],'Image Processing Toolbox'))
    error('This application requires the Matlab Image Processing Toolbox to be installed')
end

% determine the release date and use the appropriate file loading function
% newer versions > 2010b use VideoReader, older versions use mmreader
if datenum(version(1).Date) >= 734353
    reader = 'VideoReader';
else reader = 'mmreader';
end

% if post 2010a version of Matlab then look at the available video formats
% for the list
if datenum(version(1).Date) >= 734202
    %determine the available video formats
    file_formats = eval([reader '.getFileFormats']);
    format_list{1,1} = [];
    format_list{1,2} = 'All Video Files (';
    for i = 1:length(file_formats)
        format_list{1,1} = [format_list{1,1} '*.' file_formats(i).get.Extension ';'];
        format_list{1,2} = [format_list{1,2} file_formats(i).get.Description ', '];
        format_list{i+1,1} = ['*.' file_formats(i).get.Extension];
        format_list{i+1,2} = ['*.' file_formats(i).get.Extension ' - ' file_formats(i).get.Description];
    end

    format_list{1,2} = [format_list{1,2}(1:end-2) ')'];

    newlist = strcat(format_list{1,1},'*.b32;*.b8;*.mat');% add some options that we can use

    %load the avi file
    [handles.fname, handles.pname] = uigetfile(newlist, 'Pick a movie file');

else % pre 2010a mmreader function cannot use the getFileFormats method so only allow AVI files to be selected
    [handles.fname, handles.pname] = uigetfile('*.avi', 'Pick a movie file');
end

if isequal(handles.fname,0) || isequal(handles.pname,0)
    return;
end

mb = waitbar(0,'Loading Video....');
cd (handles.pname)
[~,~,Ext]=fileparts([handles.pname handles.fname]);

if strcmp(Ext,'.b32')||strcmp(Ext,'.b8')

    [BIm,Bheader] = RPread([handles.pname handles.fname]); %read in the .b32 file

    handles.ImStack = fliplr(BIm / max(max(max(BIm))));
    handles.NumFrames = size(BIm,3); % get the number of image frames
    handles.vidHeight = (Bheader.bl(2)-1)-(Bheader.ul(2)+1); % get the height in pixels of the images
    handles.vidWidth = (Bheader.br(1)-1)-(Bheader.bl(1)+1); % get the width in pixels
    handles.FrameRate = Bheader.dr;

elseif strcmp(Ext,'.mat')

    load([handles.pname handles.fname])
    handles.ImStack = TVDdata.Im;
    handles.vidHeight = double(TVDdata.Height);
    handles.vidWidth = double(TVDdata.Width);
    handles.NumFrames = double(TVDdata.Fnum);
    % sort the time data - the last timestamp is always duplicated and so
    % needs to be adjusted
    handles.TimeStamps = double(TVDdata.Time)/1000;
    handles.FrameRate = round(1/(handles.TimeStamps(end-1)/(handles.NumFrames-1)));
    handles.TimeStamps(end) = handles.TimeStamps(end-1)+(1/handles.FrameRate);

elseif strcmp(Ext,'.png')

    imgRGB = imread([handles.pname handles.fname]);
    imgCrop = imcrop(imgRGB, [290 44 977 779]);
    TVDdata.Im = im2double(im2gray(imgCrop));
    handles.ImStack = TVDdata.Im;
    handles.ImStack(:,:,2) = TVDdata.Im;
    handles.vidHeight = size(handles.ImStack,1);
    handles.vidWidth = size(handles.ImStack,2);
    handles.NumFrames = 2;
    % sort the time data - the last timestamp is always duplicated and so
    % needs to be adjusted
    handles.TimeStamps = [1; 2];
    handles.FrameRate = 1;
    %handles.Time = [0 1];

else

    switch reader


        case 'VideoReader'

            handles.movObj = VideoReader([handles.pname handles.fname]);

            % get info
            handles.vidHeight = handles.movObj.Height;
            handles.vidWidth = handles.movObj.Width;
            handles.NumFrames = handles.movObj.NumFrames;
            handles.FrameRate = handles.movObj.FrameRate;

            i=1;

            % start clean
            if isfield(handles,'ImTrack')
                handles = rmfield(handles, 'ImTrack');
            end
            if isfield(handles,'ImStack')
                handles = rmfield(handles, 'ImStack');
            end
            if isfield(handles,'ImStackOr')
                handles = rmfield(handles, 'ImStackOr');
            end

            handles.ImStack     = zeros(handles.vidHeight, handles.vidWidth, handles.NumFrames,'uint8');

            while hasFrame(handles.movObj)
                waitbar(handles.movObj.CurrentTime/handles.movObj.Duration,mb)
                if regexp(handles.movObj.VideoFormat,'RGB')
                    handles.ImStack(:,:,i) = im2gray(readFrame(handles.movObj));
                else
                    handles.ImStack(:,:,i) = readFrame(handles.movObj);
                end

                handles.ImBrightness(i) = mean(handles.ImStack(:,:,i),'all');
                i=i+1;
            end


        case 'mmreader' % pre R2010b uses mmreader
            handles.movObj = eval([reader '([handles.pname handles.fname])']);
            handles.NumFrames = handles.movObj.NumberOfFrames;
            handles.vidHeight = handles.movObj.Height;
            handles.vidWidth = handles.movObj.Width;
            handles.FrameRate = handles.movObj.NumberOfFrames/handles.movObj.Duration;

            for i = 1:handles.NumFrames
                cdata = read(handles.movObj,i);
                waitbar(i/handles.NumFrames,mb)
                % create image from movie frame
                if regexp(handles.movObj.VideoFormat,'RGB')
                    handles.ImStack(:,:,i) = rgb2gray(cdata);
                else
                    handles.ImStack(:,:,i) = cdata;
                end
                clear cdata
            end

    end

end

% check whether a mat file exists in the location with the same name with
% setting of the video (ImageDepth)
[path,name,~] = fileparts([handles.pname handles.fname]); %more elegant, people may not have necessarly mp4
if exist([path '/' name '.mat'],"file")
    TVD = load([path '/' name '.mat']);
    if isfield(TVD, 'TVDdata')
        if isfield(TVD.TVDdata,'cmPerPixY') %check whether the field exists and update scalar
            %ImageDepth = round(TVD.TVDdata.cmPerPixY*10,3); %round to 3 digits
            ImageDepth = round(cast(TVD.TVDdata.Height * TVD.TVDdata.cmPerPixY,'single'),3)*10;
            handles.ID = ImageDepth;
            set(handles.ImDepthEdit,'String',num2str(ImageDepth));
            
            % Update handles structure
            guidata(hObject, handles);
        end
    end
    clearvars path name
end

% display the path and name of the file in the filename text box
set(handles.filename,'String',[handles.pname handles.fname])

% set the string in the frame_number box to the current frame value (1)
set(handles.frame_number,'String',num2str(1))

% allows Cut_frames_before_Callback to work
handles.start_frame = 1;

% set the limits on the slider - use number of frames to set maximum (min =
% 1)
set(handles.frame_slider,'Min',1);
set(handles.frame_slider,'Max',handles.NumFrames);
set(handles.frame_slider,'Value',1);
set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 10/handles.NumFrames]);
set(handles.frame_rate,'String',handles.FrameRate(1))
set(handles.vid_width,'String',handles.vidWidth(1))
set(handles.vid_height,'String',handles.vidHeight(1))

% set the image croppable area to the maximum area
if ~isfield(handles,'crop_rect')||isempty(handles.crop_rect)

    if strcmp(Ext,'.b32')||strcmp(Ext,'.b8')

        handles.crop_rect = [Bheader.ul(1),Bheader.ul(2),...
            handles.vidWidth handles.vidHeight];

    else

        handles.crop_rect = [1 1 handles.vidWidth handles.vidHeight];

    end
end

cd(handles.pname)

% make a timeline which corresponds to the ultrasound frame data
if strcmp(Ext,'.mat')
    handles.Time = handles.TimeStamps;
elseif exist('TVD','var') && isfield(TVD, 'TVDdata')
    if isfield(TVD.TVDdata,'Time') %check if also time exists as Telemed has uncostant framerate
        handles.Time = TVD.TVDdata.Time; %note that the last timestamp is repeated in Echo Wave II recordings
    end
else
    handles.Time = (double(1/handles.FrameRate):double(1/handles.FrameRate):double(handles.NumFrames/handles.FrameRate))';
end

handles.KeyframeInd = logical(0.0);

if exist('TrackingData','var')
    chkload = questdlg('Do you want to load previous tracking?','Tracking data detected','Yes');
    if strcmp(chkload,'Yes')

        handles.Region = TrackingData.Region;
        handles.start_frame = TrackingData.start_frame;
        handles.NumFrames = TrackingData.NumFrames;
        set(handles.frame_slider,'Min',1);
        set(handles.frame_slider,'Max',handles.NumFrames);
        set(handles.frame_slider,'Value',1);
        set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 5/handles.NumFrames]);
        % set the string in the frame_number to 1
        set(handles.frame_number,'String',1);
    end
end

cd(handles.pname)

% update the image axes using show_image function (bottom)
clear_fascicle_Callback(hObject, eventdata, handles);

handles.ImStackOr = handles.ImStack;
if handles.flipimage.Value == 1%check based on flip tick box value
     handles = do_flip(hObject, eventdata, handles);
end
guidata(hObject, handles);
show_image(hObject,handles);
waitbar(1,mb)
close(mb)

% --- Outputs from this function are returned to the command line.
function varargout = UltraTimTrack_OutputFcn(~, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cut_frames_before.
function cut_frames_before_Callback(hObject, eventdata, handles)
% hObject    handle to cut_frames_before (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'movObj')||isfield(handles,'BIm')||isfield(handles,'ImStack')

    frame_no = round(get(handles.frame_slider,'Value'));

    handles.start_frame = frame_no + handles.start_frame;

    handles.NumFrames = handles.NumFrames-handles.start_frame+1;

    handles.Time = handles.Time(handles.start_frame:handles.start_frame + handles.NumFrames-1);

    set(handles.frame_slider,'Min',1);
    set(handles.frame_slider,'Max',handles.NumFrames);
    set(handles.frame_slider,'Value',1);
    set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 5/handles.NumFrames]);

    % set the string in the frame_number to 1
    set(handles.frame_number,'String',1);

    % update the image axes using show_image function (bottom)
    show_image(hObject,handles);
end

% --- Executes on button press in cut_frames_after.
function cut_frames_after_Callback(hObject, eventdata, handles)
% hObject    handle to cut_frames_after (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'movObj')||isfield(handles,'BIm')||isfield(handles,'ImStack')
    frame_no = round(get(handles.frame_slider,'Value'));

    handles.NumFrames = frame_no;

    set(handles.frame_slider,'Min',1);
    set(handles.frame_slider,'Max',handles.NumFrames);
    set(handles.frame_slider,'Value',handles.NumFrames);
    set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 5/handles.NumFrames]);

    % set the string in the frame_number to 1
    set(handles.frame_number,'String',frame_no);

    % update the image axes using show_image function (bottom)
    show_image(hObject,handles);
end

% --- Executes on slider movement.
function frame_slider_Callback(hObject, eventdata, handles)
% hObject    handle to frame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% get the current value from the slider (round to ensure it is integer)
frame_no = round(get(handles.frame_slider,'Value'));

% set the string in the frame_number box to the current frame value
set(handles.frame_number,'String',num2str(frame_no));

% update the image axes using show_image function (bottom)
show_image(hObject,handles);

% --- Executes during object creation, after setting all properties.
function frame_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.


if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function frame_number_Callback(hObject, eventdata, handles)
% hObject    handle to frame_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_number as text
%        str2double(get(hObject,'String')) returns contents of frame_number as a double

frame_no = str2num(get(handles.frame_number,'String'));
set(handles.frame_slider,'Value',round(frame_no));

% update the image axes using show_image function (bottom)
show_image(hObject,handles);

% --- Executes during object creation, after setting all properties.
function frame_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -----------------------------------------------
% -----------------------------------------------
function goto_next_frame(hObject, handles)

if get(handles.autorun_but,'BackgroundColor') == [0 1 0]

    % get the current value from the slider (round to ensure it is integer)
    frame_no = round(get(handles.frame_slider,'Value'))+1;

    if frame_no < get(handles.frame_slider,'Max')
        % set the slider
        set(handles.frame_slider,'Value',frame_no);
        % set the string in the frame_number box to the current frame value
        set(handles.frame_number,'String',num2str(frame_no));
    else  handles.run = 0;
        set(handles.autorun_but,'BackgroundColor',[1 0 0]);
    end

    % update the image axes using show_image function (bottom)
    show_image(hObject,handles);

end

% --- Executes on button press in clear_fascicle.
function[handles] = clear_fascicle_Callback(hObject, eventdata, handles)
% hObject    handle to clear_fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')
    %remove opticflow and state estimated
    if isfield(handles,'Region')
        handles = rmfield(handles,'Region');
    end
    %remove timtrack
    if isfield(handles,'geofeatures')
        handles = rmfield(handles, 'geofeatures');
    end
    %clear up also points, otherwise the are annoying see them in a
    %different region potentially
    if isfield(handles,'points')
        handles = rmfield(handles,'points');
    end
    % reset tracking

    if isfield(handles, 'ImTrack')
        handles = rmfield(handles, 'ImTrack');
    end

    % set current frame to 1
    set(handles.frame_slider,'Value',1);
    set(handles.frame_number,'String',1);

    cla(handles.length_plot); %clean fascicle length data
    cla(handles.mat_plot);%clean fascicle angle data
    cla(handles.axes1); %clean image data

    show_image(hObject,handles);
end

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Image_Callback(hObject, eventdata, handles)
% hObject    handle to Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tracking_Callback(hObject, eventdata, handles)
% hObject    handle to Tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Fascicle_Callback(hObject, eventdata, handles)
% hObject    handle to Fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Settings_Callback(hObject, eventdata, handles)
% hObject    handle to Settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_settings_Callback(hObject, eventdata, handles)
% hObject    handle to save_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ImageDepth = handles.ID;
% Needed for KL algorithm, but not KLT algorithm
%Sigma = handles.SIGMA;
%S_Step = handles.S_STEP;
Position = get(gcf,'Position');
Default_Directory = cd;

[directory, ~, ~] = fileparts(mfilename('fullpath'));

if ismac || isunix
    save([directory '/ultrasound_tracking_settings.mat'], 'ImageDepth', ...
        'Position', 'Default_Directory');
end

if ispc
    save([directory '\ultrasound_tracking_settings.mat'], 'ImageDepth',...
        'Position', 'Default_Directory');
end

msgbox('Settings saved')

% --------------------------------------------------------------------
function[handles] = menu_clear_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to menu_clear_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    if isfield(handles,'Region')
        handles = rmfield(handles,'Region');
    end

    % set current frame to 1
    set(handles.frame_slider,'Value',1);
    set(handles.frame_number,'String',1);

    cla(handles.length_plot)
    cla(handles.mat_plot)
    cla(handles.axes1) %clean image data

    show_image(hObject,handles);

end

% --------------------------------------------------------------------
function[handles] = AutoCrop_Callback(hObject, eventdata, handles)
% hObject    handle to AutoCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Im = handles.ImStack;

% relative change with time
dIm = sum(abs(diff(Im,1, 3)),3);
ndIm = dIm / max(dIm(:));

% threshold to binarize
th = .05;
ndIm(ndIm>th) = 1;
ndIm(ndIm<=th)= 0;

% find connected components in the filtered matrix.
BW2 = bwareaopen(ndIm,50);
stats = regionprops(BW2, 'BoundingBox');

% extract the bounding box information.
boundingBoxes = cat(1, stats.BoundingBox);

% calculate the overall bounding box that encompasses all smaller bounding boxes.
B = round([min(boundingBoxes(:, 1)), min(boundingBoxes(:, 2)), ...
    max(boundingBoxes(:, 1) + boundingBoxes(:, 3)) - min(boundingBoxes(:, 1)), ...
    max(boundingBoxes(:, 2) + boundingBoxes(:, 4)) - min(boundingBoxes(:, 2))]);

handles.ImStack = Im(B(2):(B(2)+B(4)-1), B(1):(B(1)+B(3)-1),:);

handles.vidHeight = size(handles.ImStack,1);
handles.vidWidth = size(handles.ImStack,2);

% clear all tracking
handles = menu_clear_tracking_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% update the image axes using show_image function (bottom)
show_image(hObject,handles);

% --------------------------------------------------------------------
function menu_crop_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_crop_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.ImStack = handles.ImStackOr;
if isfield(handles,'ImStack')

    % define current axes
    set(handles.axes1)
    axes(handles.axes1)

    % use imcrop tool to determine croppable area
    [~,handles.crop_rect] = imcrop;

    handles.crop_rect = round(handles.crop_rect);
    handles.vidHeight = handles.crop_rect(4)+1;
    handles.vidWidth = handles.crop_rect(3)+1;

    % save a copy and overwrite
    ImStackOld = handles.ImStack;
    handles.ImStack     = zeros(handles.vidHeight, handles.vidWidth, handles.NumFrames,'uint8');

    %Crop all images before updating
    for ii = handles.start_frame : handles.NumFrames
        handles.ImStack(:,:,ii) = imcrop(ImStackOld(:,:,ii),handles.crop_rect);
    end

    clearvars tmp
    % Clean axis from original image and tight axis on the cropped image
    cla
    % update the image axes using show_image function (bottom)
    show_image(hObject,handles);
    axis tight

    guidata(hObject, handles);

end

% --------------------------------------------------------------------
function menu_reset_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_reset_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    % set the image croppable area to the maximum area
    %handles.crop_rect = [1 1 handles.vidWidth handles.vidHeight];
    handles.ImStack = handles.ImStackOr; %restore original image
    % update the image axes using show_image function (bottom)
    show_image(hObject,handles);

end


% --------------------------------------------------------------------
function menu_set_depth_Callback(hObject, eventdata, handles)
% hObject    handle to menu_set_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    if ~isfield(handles, 'ID')
        handles.ID = 60.4;
    end

    A = inputdlg('Enter image depth', 'Image Depth', 1, {num2str(handles.ID)});
    handles.ID = str2double(A{1});

    set(handles.ImDepthEdit, 'String',A{1})
    % Update handles structure
    guidata(hObject, handles);

end


% --------------------------------------------------------------------
function [handles] = menu_load_fascicle_Callback(hObject, eventdata, handles,varargin)
% hObject    handle to menu_load_fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if the 'file' input is provided 
if nargin < 4
    % 'file' input is not provided, handle accordingly
    fileFas = []; % You can set a default file or leave it empty
else
    % 'file' input is provided in case of process folder
    fileFas = varargin{1};
end

%make trasparent
set(handles.S, 'EdgeAlpha',0,'FaceAlpha',0.1,'InteractionsAllowed','none')
set(handles.D, 'EdgeAlpha',0,'FaceAlpha',0.1,'InteractionsAllowed','none')


if isfield(handles,'ImStack')

   % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'))+handles.start_frame-1;
    if isempty(fileFas)
        [fname, pname] = uigetfile('*.mat','Load tracking MAT file');
        load([pname fname]);
    else
        load(fileFas);
    end
    
    %just be sure such necessary data exists
    if exist('Fdat','var') && isfield(Fdat,'Region') && isfield(Fdat.Region,'Fascicle')

        for i = 1:length(Fdat.Region)

            for k = 1:length(Fdat.Region(i).Fascicle)

                handles.Region(i).Fascicle(k).fas_x{frame_no} = Fdat.Region(i).Fascicle(k).fas_x;
                handles.Region(i).Fascicle(k).fas_y{frame_no} = Fdat.Region(i).Fascicle(k).fas_y;

                %in case some wants to load a single fascicle from a series
                %tracked
                handles.Region(i).Fascicle(k).current_xy(1,1) = handles.Region(i).Fascicle(k).fas_x{frame_no}(1);
                handles.Region(i).Fascicle(k).current_xy(1,2) = handles.Region(i).Fascicle(k).fas_y{frame_no}(1);
                handles.Region(i).Fascicle(k).current_xy(2,1) = handles.Region(i).Fascicle(k).fas_x{frame_no}(2);
                handles.Region(i).Fascicle(k).current_xy(2,2) = handles.Region(i).Fascicle(k).fas_y{frame_no}(2);
                
                %Add ROI data
                handles.Region(i).ROIx{frame_no} = Fdat.Region(i).ROIx;
                handles.Region(i).ROIy{frame_no} = Fdat.Region(i).ROIy;

                try
                    handles.Region(i).deep_x{frame_no} =  Fdat.Region(i).deep_x ;
                    handles.Region(i).deep_y{frame_no} =  Fdat.Region(i).deep_y ;
                    handles.Region(i).sup_x{frame_no} =   Fdat.Region(i).sup_x ;
                    handles.Region(i).sup_y{frame_no} =  Fdat.Region(i).sup_y ;
                catch
                    fprintf('Apo pts defined by ROI because not specify in the loaded fascicle file!\n')
                    handles.Region(i).deep_x{frame_no} =  Fdat.Region(i).ROIx([2,3]);
                    handles.Region(i).deep_y{frame_no}=  Fdat.Region(i).ROIy([2,3]);
                    handles.Region(i).sup_x{frame_no} =  Fdat.Region(i).ROIx([end,end-1]);
                    handles.Region(i).sup_y{frame_no} =  Fdat.Region(i).ROIy([end,end-1]);
                end
                % handles.Region(i).fas_pen(frame_no,k) = atan2d(abs(diff(handles.Region(i).Fascicle(k).fas_y{frame_no})),abs(diff(handles.Region(i).Fascicle(k).fas_x{frame_no})));
                % scalar = handles.ID;%/handles.vidHeight;
                % handles.Region(i).fas_length(frame_no,k) = scalar*sqrt(diff(handles.Region(i).Fascicle(k).fas_y{frame_no}).^2 + diff(handles.Region(i).Fascicle(k).fas_x{frame_no}).^2);
                % calculate fascicle length and penn
                handles = calc_fascicle_length_and_pennation(handles,frame_no);

                if ~isfield(handles.Region(i).Fascicle(k),'analysed_frames')
                    handles.Region(i).Fascicle(k).analysed_frames = frame_no;
                else
                    handles.Region(i).Fascicle(k).analysed_frames = sort([handles.Region(i).Fascicle(k).analysed_frames frame_no]);
                end
            end


        end

      
    else

        warndlg('Fascicle data not found','ERROR')

    end

    show_data(hObject,handles);
    show_image(hObject,handles);
end

% --------------------------------------------------------------------
function menu_save_fascicle_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'movObj')||isfield(handles,'BIm')||isfield(handles,'ImStack')

    % find current frame number from slider in relation to starting frame
    frame_no = round(get(handles.frame_slider,'Value'))+handles.start_frame-1;
    try
        for i = 1:length(handles.Region)
            for k = 1:length(handles.Region(i).Fascicle)
    
                if isfield(handles.Region(i).Fascicle(k),'fas_x')
    
                    Fdat.Region(i).Fascicle(k).fas_x = handles.Region(i).Fascicle(k).fas_x{frame_no};
                    Fdat.Region(i).Fascicle(k).fas_y = handles.Region(i).Fascicle(k).fas_y{frame_no};
    
                end
    
            end
            % Save apo data (OK)
            Fdat.Region(i).sup_x = handles.Region(i).sup_x{frame_no};
            Fdat.Region(i).sup_y = handles.Region(i).sup_y{frame_no};
            Fdat.Region(i).deep_x = handles.Region(i).deep_x{frame_no};
            Fdat.Region(i).deep_y = handles.Region(i).deep_y{frame_no};
            % Save fascicle length and pennation
            Fdat.Region(i).fas_length = handles.Region(i).fas_length(frame_no);
            Fdat.Region(i).fas_pen = handles.Region(i).fas_pen(frame_no);
            % Also save the ROI data
            Fdat.Region(i).ROIx = handles.Region(i).ROIx{frame_no};
            Fdat.Region(i).ROIy = handles.Region(i).ROIy{frame_no};
            %Fdat.Region(i).ROI = handles.Region(i).ROIp{frame_no};
        end
    
        % save to file
        [fname, pname] = uiputfile('*.mat', 'Save fascicle to MAT file');
        save([pname fname],'Fdat');
    catch
        warndlg('There is no fascicle to save at this frame!','ERROR')
    end
end

% --------------------------------------------------------------------
function menu_save_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')
    [fileout, pathout, FI] = uiputfile('*.tif', 'Save video as');

    if FI > 0
        IM_out = getframe(handles.axes1);
        imwrite(IM_out.cdata,[pathout fileout], 'TIFF');
    end

end

% --------------------------------------------------------------------
function save_video_Callback(hObject, eventdata, handles)
% hObject    handle to save_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    filename = [handles.pname, handles.fname(1:end-4), '_tracked_Q=',strrep(num2str(handles.Q),'.','')];
    vidObj = VideoWriter(filename,'MPEG-4');
    vidObj.FrameRate = handles.FrameRate;
    open(vidObj);

    h = waitbar(0,['Saving frame 1/', num2str(handles.NumFrames)],'Name','Saving to video file...');
    i = 1;
    j = 1;

    for f = handles.start_frame:1:get(handles.frame_slider,'Max')

        if isfield(handles, 'Region')
            % create tracked frame
            d = round(size(handles.ImStack,2)/2);

            ZeroPadL = 200*ones(size(handles.ImStack,1), ceil(size(handles.ImStack,2)/2),'uint8');
            ZeroPadR = 200*ones(size(handles.ImStack,1), floor(size(handles.ImStack,2)/2),'uint8');

            % add padding
            currentImage = [ZeroPadL, handles.ImStack(:,:,f), ZeroPadR];

            % add fascicle

            currentImage = insertShape(currentImage,'line',[handles.Region(i).Fascicle(j).fas_x{f}(1)+d, handles.Region(i).Fascicle(j).fas_y{f}(1), ...
                handles.Region(i).Fascicle(j).fas_x{f}(2)+d,handles.Region(i).Fascicle(j).fas_y{f}(2)], 'LineWidth',5, 'Color','red');

            currentImage = insertMarker(currentImage,[handles.Region(i).Fascicle(j).fas_x{f}(1)+d, handles.Region(i).Fascicle(j).fas_y{f}(1);...
                handles.Region(i).Fascicle(j).fas_x{f}(2)+d, handles.Region(i).Fascicle(j).fas_y{f}(2)], 'o', 'Color','red','size',5);

            % add aponeurosis
            currentImage = insertShape(currentImage,'line',[handles.Region(i).sup_x{f}(1)+d, handles.Region(i).sup_y{f}(1), ...
                handles.Region(i).sup_x{f}(2)+d,handles.Region(i).sup_y{f}(2)], 'LineWidth',5, 'Color','blue');

            currentImage = insertShape(currentImage,'line',[handles.Region(i).deep_x{f}(1)+d, handles.Region(i).deep_y{f}(1), ...
                handles.Region(i).deep_x{f}(2)+d,handles.Region(i).deep_y{f}(2)], 'LineWidth',5, 'Color','green');

            % add ROI
            currentImage = insertShape(currentImage,'Polygon',[handles.Region(i).ROIx{f}(1)+d, handles.Region(i).ROIy{f}(1), ...
                handles.Region(i).ROIx{f}(2)+d, handles.Region(i).ROIy{f}(2),handles.Region(i).ROIx{f}(3)+d, handles.Region(i).ROIy{f}(3),...
                handles.Region(i).ROIx{f}(4)+d, handles.Region(i).ROIy{f}(4),handles.Region(i).ROIx{f}(5)+d, handles.Region(i).ROIy{f}(5)],'LineWidth',1, 'Color','red');

            % save
            ImTrack = currentImage;
        else
            ImTrack = handles.ImStack(:,:,f);
        end

        if isfield(handles, 'S')
            if isvalid(handles.S)
                % show region
                spos = ceil([handles.S.Position(1:2) ceil(size(handles.ImStack,2)/2) handles.S.Position(4)]);
                dpos = ceil([handles.D.Position(1:2) ceil(size(handles.ImStack,2)/2) handles.D.Position(4)]);

                ImTrack(spos(2):(spos(2)+spos(4)),spos(1):(spos(1)+spos(3)),3) = 230;
                ImTrack(dpos(2):(dpos(2)+dpos(4)),dpos(1):(dpos(1)+dpos(3)),2) = 230;
            end
        end

        F = ImTrack;
        writeVideo(vidObj,F)

        frac_progress = f/handles.NumFrames;
        waitbar(frac_progress,h, ['Processing frame ', num2str(f), '/', num2str(get(handles.frame_slider,'Max'))])

    end
    %     end
    close(vidObj)
    close(h)

end


% --------------------------------------------------------------------
function menu_process_all_Callback(hObject, eventdata, handles)
% hObject    handle to menu_process_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

process_all_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Save_As_Txt_Callback(hObject, eventdata, handles)
% hObject    handle to Save_As_Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'Region')

    for i = 1:length(handles.Region)


        if isfield(handles,'movObj')||isfield(handles,'BIm')||isfield(handles,'ImStack') && isfield(handles.Region(i),'fas_length')

            time = handles.Time;

            %determine any non-zero entries in fascicle length array
            nz = logical(handles.Region(i).fas_length(:,1) ~= 0);

            T = time(nz)';

            if isfield(handles.Region(i),'fas_length_corr') && ~isempty(handles.Region(i).fas_length_corr)
                R(i).FL = handles.Region(i).fas_length_corr(nz,:)';
                R(i).PEN = handles.Region(i).fas_pen_corr(nz,:)';

                for j = 1:size(handles.Region(i).fas_length_corr,2)
                    if sum(handles.Region(i).fas_length_corr(nz,j)) == 0

                        R(i).FL(j,:) = handles.Region(i).fas_length(nz,j)';
                        R(i).PEN(j,:) = handles.Region(i).fas_pen(nz,j)';

                    end
                end

            else R(i).FL = handles.Region(i).fas_length(nz,:)';
                R(i).PEN = handles.Region(i).fas_pen(nz,:)';
            end

            col_head{i} = ['Time\tFascicle Length R' num2str(i) '_F1\tFascicle Angle R' num2str(i) '_F1\t'];
            col_type{i} = '%3.4f\t%3.6f\t%3.6f';

            data_out{i} = [T' R(i).FL(1,:)' R(i).PEN(1,:)'];

            if size(R(i).FL',2) > 1
                for k = 2:size(R(i).FL',2)

                    data_out{i} = [data_out{i} R(i).FL(k,:)' R(i).PEN(k,:)'];
                    regionlabel = num2str(i);
                    faslabel = num2str(k);
                    col_head{i} = [col_head{i} 'Fascicle Length R' regionlabel '_F' faslabel '\tPennation Angle R' regionlabel '_F' faslabel '\t'];
                    col_type{i} = [col_type{i} '\t%3.6f\t%3.6f'];
                end
            end
        end
    end

    header = [];
    type = [];
    dout = [];
    for i = 1:length(handles.Region)
        header = [header col_head{i}];
        type = [type col_type{i}];
        dout = [dout data_out{i}];
    end


    header = [header '\n'];
    %     type = [type '\n'];

    fileout_suggest = [handles.pname handles.fname(1:end-3) 'txt'];
    [fileout, pathout] = uiputfile(fileout_suggest,'Save fascicle data as...');
    cd(pathout);
    % open a new text file for writing
    fid = fopen(fileout,'w');
    fprintf(fid,header);
    dlmwrite(fileout,dout,'-append','delimiter','\t')
    %     fprintf(fid,type,dout);

    fclose(fid);
end


% --------------------------------------------------------------------
function Save_As_Mat_Callback(hObject, eventdata, handles)
% hObject    handle to Save_As_Mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'Region')

    for i = 1:length(handles.Region)

        if isfield(handles,'movObj')||isfield(handles,'BIm')||isfield(handles,'ImStack') && isfield(handles.Region(i),'fas_length')

            time = handles.Time;
            %determine any non-zero entries in fascicle length array
            nz = logical(handles.Region(i).fas_length(handles.start_frame:end,1) ~= 0);

            T = time(nz)';

            TrackingData.res = handles.ID;
            TrackingData.start_frame = handles.start_frame;
            TrackingData.NumFrames = handles.NumFrames;
            %info about tracking for replication purposes
            TrackingData.ProcessingTime = handles.ProcessingTime; %two
            TrackingData.BlockSize = handles.BlockSize;
            TrackingData.Parallel = handles.do_parfor.Value;
            TrackingData.info = "Processing [TimTrack; Opticflow], %%\nBlockSize [width; height], %%\nGains [Apo, Position, Angle]";
            TrackingData.S = handles.S;
            TrackingData.D = handles.D;

            if isfield(handles,'R')
                Fdat.R = handles.R;
            end

            if isfield(handles,'Frequency')
                Fdat.Frequency = handles.Frequency;
            end

            if isfield(handles, 'geofeatures')
                Fdat.geofeatures = handles.geofeatures;
            end

            Fdat.Region(i) = handles.Region(i);
            Fdat.Region(i).FL = handles.Region(i).fas_length(nz,:)';
            Fdat.Region(i).PEN = handles.Region(i).fas_pen(nz,:)';
            if isfield(handles.Region,'fas_ang') %this exists only when estimator runs
                Fdat.Region(i).ANG = handles.Region(i).fas_ang(nz,:)';
            end
            Fdat.Region(i).Time = time(nz);
        end
    end
end

filename = [handles.pname, handles.fname(1:end-4), '_tracked_Q=',strrep(num2str(handles.Q),'.','')];
save(filename,'TrackingData','Fdat');

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if strcmp(eventdata.Key,'rightarrow')
    % get the current value from the slider (round to ensure it is integer)
    % and add one to make the new frame number
    frame_no = round(get(handles.frame_slider,'Value'))+1;

    if frame_no < get(handles.frame_slider,'Max')
        % set the slider
        set(handles.frame_slider,'Value',frame_no);
        % set the string in the frame_number box to the current frame value
        set(handles.frame_number,'String',num2str(frame_no));

        show_image(hObject,handles);
    end
end

if strcmp(eventdata.Key,'leftarrow')
    % get the current value from the slider (round to ensure it is integer)
    % and subtract 1 to make the new frame number
    frame_no = round(get(handles.frame_slider,'Value'))-1;

    if frame_no > 1
        % set the slider
        set(handles.frame_slider,'Value',frame_no);
        % set the string in the frame_number box to the current frame value
        set(handles.frame_number,'String',num2str(frame_no));

        show_image(hObject,handles);
    end
end


% --- Executes on button press in PlayButton.
function PlayButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'Interruptible','on');
% show_data(hObject, handles);

frame_no = round(get(handles.frame_slider,'Value'));
end_frame = get(handles.frame_slider,'Max');
handles.stop = 0;
guidata(hObject,handles);

for f = frame_no:end_frame
    stop = get(handles.StopButton,'Value');

    if stop
        set(handles.StopButton,'Value',0.0);
        return
    else
        set(handles.frame_slider,'Value',f);
        set(handles.frame_number,'String',num2str(f));

        % update image
        handles = show_image(hObject,handles);
        drawnow
    end
end

% --- Executes on button press in StopButton.
function StopButton_Callback(hObject, eventdata, handles)
% hObject    handle to StopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.StopButton,'Value',1.0)


% --------------------------------------------------------------------
function Load_All_Tracked_Frames_Callback(hObject, eventdata, handles)
% hObject    handle to Load_All_Tracked_Frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%load the avi file

menu_clear_tracking_Callback(hObject, eventdata, handles); % clear any current tracking

%select file
[fname, pname] = uigetfile('*.mat', 'Pick a .MAT file');
load([pname fname],'TrackingData','Fdat');

handles.Region = Fdat.Region;
handles.start_frame = TrackingData.start_frame;
handles.NumFrames = TrackingData.NumFrames;
handles.S = TrackingData.S;
handles.D = TrackingData.D;
handles.geofeatures = Fdat.geofeatures;
set(handles.frame_slider,'Min',1);
set(handles.frame_slider,'Max',handles.NumFrames);
set(handles.frame_slider,'Value',1);
set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 5/handles.NumFrames]);
% set the string in the frame_number to 1
set(handles.frame_number,'String',1);

show_image(hObject,handles);
show_data(hObject, handles); %update plots


% --- Executes on button press in zoominvideo.
function zoominvideo_Callback(hObject, eventdata, handles)
% hObject    handle to zoominvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% function too zoom the image by adjusting the axes' limits
axes(handles.axes1)
cx = get(gca,'XLim');
cy = get(gca,'YLim');

xrange = cx(2)-cx(1);
yrange = cy(2)-cy(1);
px = xrange/10;
py = yrange/10;

nx = [cx(1)+(0.5*px),cx(2)-(0.5*px)];
ny = [cy(1)+(0.5*py),cy(2)-(0.5*py)];


set(gca,'XLim',nx,'YLim',ny);
handles.Vidax_X = nx;
handles.Vidax_Y = ny;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in zoomoutvideo.
function zoomoutvideo_Callback(hObject, eventdata, handles)
% hObject    handle to zoomoutvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
cx = get(gca,'XLim');
cy = get(gca,'YLim');

xrange = cx(2)-cx(1);
yrange = cy(2)-cy(1);
px = xrange/10;
py = yrange/10;

nx = [cx(1)-(0.5*px),cx(2)+(0.5*px)];
ny = [cy(1)-(0.5*py),cy(2)+(0.5*py)];

set(gca,'XLim',nx,'YLim',ny);
handles.Vidax_X = nx;
handles.Vidax_Y = ny;
% Update handles structure
guidata(hObject, handles);

% --- Executes on value changed
function ImDepthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ImDepthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImDepthEdit as text
%        str2double(get(hObject,'String')) returns contents of ImDepthEdit as a double
scalarOld = handles.ID;
handles.ID = str2double(get(handles.ImDepthEdit,'String'));

for i = 1:length(handles.Region)
    if isfield(handles.Region(i),'Fascicle')
        if isfield(handles.Region(i),'fas_length')

            if ~isempty(handles.Region(i).fas_length)
                FL = handles.Region(i).fas_length;
                FL = FL ./ scalarOld;
                FL = FL .* handles.ID;
                handles.Region(i).fas_length = FL;
            end
        end
        show_data(hObject, handles); %update plots
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ImDepthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImDepthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function show_data(hObject, handles)
% difference with show_image is that this is called once (not per frame)

if isfield(handles, 'Region')
    for i = 1:length(handles.Region)
        if isfield(handles.Region(i),'fas_length')
            if ~isempty(handles.Region(i).fas_length)

                dt = 1/handles.FrameRate;
                FL = handles.Region(i).fas_length;
                PEN = handles.Region(i).fas_pen;
                time = 0:dt:((handles.NumFrames+handles.start_frame-2)*dt);
                
                if ~isfield(handles.Region(i), 'fas_length_manual')
                    handles.Region(i).fas_length_manual = nan(size(time));
                    handles.Region(i).fas_pen_manual = nan(size(time));
                end

                FLm = handles.Region(i).fas_length_manual;
                PENm = handles.Region(i).fas_pen_manual;

                axes(handles.length_plot); 
                hold off;
                plot(handles.length_plot,time,FL,'r', time, FLm, 'mx','linewidth',2);
                set(handles.length_plot,'ylim',[min(FL)*0.85 max(FL)*1.15],'xlim', [0 max(time)],'box','off'); %set axis 15% difference of min and and value,easier to read
                xlabel('Time (s)');
                ylabel('Fascicle Length (mm)');

                axes(handles.mat_plot);
                plot(handles.mat_plot,time,PEN,'r', time, PENm, 'mx','linewidth',2);
                set(handles.mat_plot,'ylim',[min(PEN)*0.85 max(PEN)*1.15], 'xlim', [0 max(time)],'box','off'); %set axis 15% difference of min and and value,easier to read
                xlabel('Time (s)');
                ylabel('Fascicle angle (deg)');

            end
        end
    end
end

%---------------------------------------------------------
% Function to show image with appropriate image processing
%---------------------------------------------------------
function handles = show_image(hObject, handles)

if isfield(handles,'ImStack')

    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value')) + handles.start_frame - 1;

    i = 1;
    j = 1;

    if isfield(handles, 'Region')
        if length(handles.Region(i).Fascicle(j).fas_x) >= frame_no
            f = frame_no;

            fasx = handles.Region(i).Fascicle(j).fas_x{f};
            fasy = handles.Region(i).Fascicle(j).fas_y{f};

            % create analyzed frame
            d = round(size(handles.ImStack,2)/2);

            ZeroPadL = 200*ones(size(handles.ImStack,1), ceil(size(handles.ImStack,2)/2),'uint8');
            ZeroPadR = 200*ones(size(handles.ImStack,1), floor(size(handles.ImStack,2)/2),'uint8');

            % add padding
            currentImage = [ZeroPadL, handles.ImStack(:,:,frame_no), ZeroPadR];

            % add fascicle
            currentImage = insertShape(currentImage,'line',[fasx(1)+d, fasy(1), ...
                fasx(2)+d,fasy(2)], 'LineWidth',5, 'Color','red');

            currentImage = insertMarker(currentImage,[fasx(1)+d, fasy(1);...
                fasx(2)+d, fasy(2)], 'o', 'Color','red','size',5);

            if isfield(handles.Region(i).Fascicle(j), 'fas_y_end') && ~isempty(handles.Region(i).Fascicle(j).fas_y_end{frame_no})
                fasx_end = [handles.Region(i).Fascicle(j).fas_x_end{f}];
                fasy_end = [handles.Region(i).Fascicle(j).fas_y_end{f}];

                % add fascicle
                currentImage = insertShape(currentImage,'line',[fasx_end(1)+d, fasy_end(1), ...
                    fasx_end(2)+d,fasy_end(2)], 'LineWidth',5, 'Color','red');

                currentImage = insertMarker(currentImage,[fasx_end(1)+d, fasy_end(1);...
                    fasx_end(2)+d, fasy_end(2)], 'o', 'Color','red','size',5);
            end


            if isfield(handles.Region(i).Fascicle(j), 'fas_x_manual') && length(handles.Region(i).Fascicle(j).fas_x_manual) >= frame_no
                if ~isempty(handles.Region(i).Fascicle(j).fas_x_manual{frame_no})

                    fasx_manual = [handles.Region(i).Fascicle(j).fas_x_manual{f}];
                    fasy_manual = [handles.Region(i).Fascicle(j).fas_y_manual{f}];

                    % add fascicle
                    currentImage = insertShape(currentImage,'line',[fasx_manual(1)+d, fasy_manual(1), ...
                        fasx_manual(2)+d,fasy_manual(2)], 'LineWidth',5, 'Color','magenta');

                    currentImage = insertMarker(currentImage,[fasx_manual(1)+d, fasy_manual(1);...
                        fasx_manual(2)+d, fasy_manual(2)], 'o', 'Color','magenta','size',5);
                end
            end

            % add aponeurosis
            currentImage = insertShape(currentImage,'line',[handles.Region(i).sup_x{f}(1)+d, handles.Region(i).sup_y{f}(1), ...
                handles.Region(i).sup_x{f}(2)+d,handles.Region(i).sup_y{f}(2)], 'LineWidth',5, 'Color','blue');

            currentImage = insertShape(currentImage,'line',[handles.Region(i).deep_x{f}(1)+d, handles.Region(i).deep_y{f}(1), ...
                handles.Region(i).deep_x{f}(2)+d,handles.Region(i).deep_y{f}(2)], 'LineWidth',5, 'Color','green');

            if isfield(handles,'points')
                if ~isempty(handles.points{f})

                    currentImage = insertMarker(currentImage,[handles.points{f}(:,1)+d, handles.points{f}(:,2)], '+', 'Color','red','size',2);
                    currentImage = insertText(currentImage, [10 10], ['Number of feature points: ' ,num2str(length(handles.points{f}))],'BoxColor','white');
                end
            end

            % add ROI
            currentImage = insertShape(currentImage,'Polygon',[handles.Region(i).ROIx{f}(1)+d, handles.Region(i).ROIy{f}(1), ...
                handles.Region(i).ROIx{f}(2)+d, handles.Region(i).ROIy{f}(2),handles.Region(i).ROIx{f}(3)+d, handles.Region(i).ROIy{f}(3),...
                handles.Region(i).ROIx{f}(4)+d, handles.Region(i).ROIy{f}(4),handles.Region(i).ROIx{f}(5)+d, handles.Region(i).ROIy{f}(5)],'LineWidth',1, 'Color','red');

            % save
            ImTrack = currentImage;

            % show region
            set(handles.S, 'Position', [1 handles.parms.apo.super.cut(1)*handles.vidHeight floor(handles.vidWidth/2) diff(handles.parms.apo.super.cut)*handles.vidHeight])
            set(handles.D, 'Position', [1 handles.parms.apo.deep.cut(1)*handles.vidHeight floor(handles.vidWidth/2) diff(handles.parms.apo.deep.cut)*handles.vidHeight])

        end
    end


    % show the image
    if isfield(handles,'ImStack')

        if exist('ImTrack','var')
            Im = ImTrack;
        else
            Im = handles.ImStack(:,:,frame_no);
        end

        if ~isfield(handles, 'image') || ~isvalid(handles.image)
            axes(handles.axes1)
            handles.image = image(Im);
            colormap(gray(256));
            axis off;
            axis equal

        else
            set(handles.image, 'CData',Im);
        end
    end

    % if we're showing original, show the aponeurosis regions
    if ~isfield(handles,'ImTrack') && (~isfield(handles, 'S') || ~isvalid(handles.S))
        handles.S = images.roi.Rectangle(handles.axes1,'position', [1 handles.parms.apo.super.cut(1)*handles.vidHeight handles.vidWidth diff(handles.parms.apo.super.cut)*handles.vidHeight],'color','blue');
        handles.D = images.roi.Rectangle(handles.axes1,'position', [1 handles.parms.apo.deep.cut(1)*handles.vidHeight handles.vidWidth diff(handles.parms.apo.deep.cut)*handles.vidHeight],'color','green');
    end


    % Update Im and NIm
    handles.Im = Im;
    handles.NIm = handles.Im;

    % remove previous vertical lines
    children = get(handles.length_plot, 'children');
    if length(children) > 2
        delete(children(1));
    end

    children = get(handles.mat_plot, 'children');
    if length(children) > 2
        delete(children(1));
    end


    if isfield(handles,'Region')

        FL = handles.Region(1).fas_length(:);
        PEN = handles.Region(1).fas_pen(:);

        % add new vertical lines
        dt = 1/handles.FrameRate;
        time = 0:dt:((handles.NumFrames+handles.start_frame-2)*dt);
        line(handles.length_plot, 'xdata', time(frame_no) * ones(1,2), 'ydata', [.85*min(FL) 1.15*max(FL)],'color',[0 0 0]);
        line(handles.mat_plot, 'xdata', time(frame_no) * ones(1,2), 'ydata', [.85*min(PEN) 1.15*max(PEN)],'color', [0 0 0]);

    end

    % Update handles structure
    guidata(hObject, handles);

end

% --- Executes on button press in process_all.
function[handles] = process_all_Callback(hObject, eventdata, handles)
% hObject    handle to process_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% detect first frame if required or last frame if backward tracking
% & set the slider and frame number for clarity

if handles.trackbck_chkBox.Value == 0
    frame_no = 1;%handles.start_frame;
    set(handles.frame_slider,'Value',frame_no);
    set(handles.frame_number,'String',num2str(frame_no));

else

    frame_no = handles.NumFrames;
    set(handles.frame_slider,'Value',frame_no);
    set(handles.frame_number,'String',num2str(handles.NumFrames));
end 

% Update handles structure
guidata(hObject, handles);
%update image so people have a feedback regarding forward/and or backward
show_image(hObject, handles);

if ~isfield(handles,'Region') || isnan(handles.Region(1).fas_length(frame_no+handles.start_frame-1))%fas length is in the correct index now
    handles = Auto_Detect_Callback(hObject, eventdata, handles);
end
% Run TimTrack
if contains(handles.ROItype, 'Hough')
    handles = process_all_TimTrack(hObject, eventdata, handles);
end

%if isfield(handles, 'geofeatures') %remove because this impose to RUN with
%hough
% Run UltraTrack (note: includes state estimation on ROI)
try
    if handles.trackbck_chkBox.Value == 1
            handles = process_all_UltraTrack_backwards(hObject, eventdata, handles);
    else
            handles = process_all_UltraTrack(hObject, eventdata, handles);
    end

    if contains(handles.ROItype, 'Hough')
         % % State estimation
         handles = do_state_estimation(hObject, eventdata, handles);
    end
catch    
    %close %close the waitbar <-- not the best because the wbar is a local var into process_all_Ultratrack so it's just because we are fast
    fprintf('Nothing tracked\n');
end

% update the image axes using show_image function (bottom)
show_data(hObject, handles);

% Update handles structure
guidata(hObject, handles);

%%%% Ultratrack (KLT optic flow)
function[handles] = process_all_UltraTrack(hObject, eventdata, handles)

%need to flip here for optic flow?
% if handles.trackbck_chkBox.Value == 1
%     handles.geofeatures = flip(handles.geofeatures);
%     handles.ImStack = flip(handles.ImStack,3);
% end

im1 = handles.ImStack(:,:,handles.start_frame+1);
h = waitbar(0,['Processing frame 1/', num2str(handles.NumFrames)],'Name','Running UltraTrack...'); 

tstart = tic;

frames = handles.start_frame:(handles.start_frame + handles.NumFrames-1);
% numIterations = length(frames);

for i = 1:length(handles.Region)
    for f = frames

        % extract image
        im = handles.ImStack(:,:,f);

        % make a copy
        I_fmasked = im;

        if isfield(handles,'geofeatures') && strcmp(handles.ROItype, 'Hough - local')
            M = zeros(size(im1,1), size(im1,2), handles.parms.fas.npeaks);

            for j = 1:handles.parms.fas.npeaks
                x1 = handles.geofeatures(f).x(j,1) * handles.imresize_fac;
                y1 = handles.geofeatures(f).y(j,1) * handles.imresize_fac;

                x2 = handles.geofeatures(f).x(j,2) * handles.imresize_fac;
                y2 = handles.geofeatures(f).y(j,2) * handles.imresize_fac;

                dy = 5; %what is 5?

                ROIx = [x1 x1 x2 x2 x1];
                ROIy = [y1-dy y1+dy y2+dy y2-dy y1-dy]';

                if sum(isfinite(ROIx)) == length(ROIx) && sum(isfinite(ROIy)) == length(ROIy)
                    M(:,:,j) = poly2mask(ROIx,ROIy, size(im1,1), size(im1,2));
                end
            end

            % mask
            fmask = sum(M,3);
            fmask(fmask>1) = 1;

            I_fmasked(fmask~=1) = 0;
        end


        % get current ROI
        if contains(handles.ROItype, 'Hough')
            I_amasked = im;

            n = handles.vidWidth;
            ROIx = [1 1 n n 1]';
            super_apo = handles.geofeatures(f).super_pos';
            deep_apo = handles.geofeatures(f).deep_pos';
            thickness = deep_apo - super_apo;

            r = .1;
            ROIy_fcor = round([super_apo(1)+thickness(1)*r; deep_apo-thickness*r; super_apo([2,1])+thickness([2,1])*r]);

            % mask
            fmask = poly2mask(ROIx, ROIy_fcor, size(im,1), size(im,2));
            I_fmasked(fmask~=1) = 0;

            m = handles.vidHeight;

            s = handles.parms.apo.super.cut;
            ROIys = [s(1) s(2) s(2) s(1) s(1)] * m;

            d = handles.parms.apo.deep.cut;
            ROIyd = [d(1) d(2) d(2) d(1) d(1)] * m;

            dmask =  poly2mask(ROIx, ROIyd, size(im,1), size(im,2));
            smask =  poly2mask(ROIx, ROIys, size(im,1), size(im,2));

            amask = dmask + smask;
            amask(amask > 1) = 1;
            I_amasked(amask~=1) = 0;
        end


        for j = 1:length(handles.Region(i).Fascicle)

            if f == handles.start_frame
                % detect points
                fpoints = detectMinEigenFeatures(I_fmasked,'FilterSize',11, 'MinQuality', 0.005);

                if contains(handles.ROItype, 'Hough')
                    fpoints = fpoints.selectStrongest(300);
                end

                % get location
                fpoints = double(fpoints.Location);

                % must be in ROI
                ROIy = handles.Region(i).ROIy{f};
                ROIx = handles.Region(i).ROIx{f};
                inPoints = inpolygon(fpoints(:,1), fpoints(:,2), ROIx, ROIy);
                fpoints = fpoints(inPoints,:);

                % define fascicle tracker
                fpointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',handles.BlockSize);
                initialize(fpointTracker,fpoints,im);

                if contains(handles.ROItype, 'Hough')
                    % define aponeurosis tracker
                    apoints = detectMinEigenFeatures(I_amasked,'FilterSize',11, 'MinQuality', 0.005);
                    apoints =  double(apoints.Location);
                    apointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',handles.BlockSize);
                    initialize(apointTracker,apoints,im);
                end
            else

                if ~exist('fpointTracker','var')

                    fpoints = detectMinEigenFeatures(I_fmasked,'FilterSize',11, 'MinQuality', 0.005);

                    if contains(handles.ROItype, 'Hough')
                        fpoints = fpoints.selectStrongest(300);
                    end

                    % get location
                    fpoints = double(fpoints.Location);

                    % must be in ROI
                    ROIy = handles.Region(i).ROIy{f};
                    ROIx = handles.Region(i).ROIx{f};
                    inPoints = inpolygon(fpoints(:,1), fpoints(:,2), ROIx, ROIy);
                    fpoints = fpoints(inPoints,:);

                    % define fascicle tracker
                    fpointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',handles.BlockSize);
                    initialize(fpointTracker,fpoints,im);

                    if contains(handles.ROItype, 'Hough')
                        % define aponeurosis tracker
                        apoints = detectMinEigenFeatures(I_amasked,'FilterSize',11, 'MinQuality', 0.005);
                        apoints =  double(apoints.Location);
                        apointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',handles.BlockSize);
                        initialize(apointTracker,apoints,im);
                    end

                end
                % Compute the flow and new roi
                [fpointsNew, isFound] = step(fpointTracker, im);
                [wf,~] = estimateGeometricTransform2D(fpoints(isFound,:), fpointsNew(isFound,:), 'affine', 'MaxDistance',50);
                handles.Region(i).warp(:,:,f-1) = wf;

                if contains(handles.ROItype, 'Hough')
                    % Compute the flow and new roi
                    [apointsNew, isFound] = step(apointTracker, im);
                    [wa,~] = estimateGeometricTransform2D(apoints(isFound,:), apointsNew(isFound,:), 'affine', 'MaxDistance',50);
                    handles.Region(i).awarp(:,:,f-1) = wa;
                end

                % apply the warp to fascicles
                fas_prev = [handles.Region(i).Fascicle(j).fas_x{f-1} handles.Region(i).Fascicle(j).fas_y{f-1}];
                fas_new = transformPointsForward(wf, fas_prev);

                % save
                handles.Region(i).Fascicle(j).fas_x{f} = fas_new(:,1);
                handles.Region(i).Fascicle(j).fas_y{f} = fas_new(:,2);

                % make a copy
                handles.Region(i).Fascicle(j).fas_x_original{f} = handles.Region(i).Fascicle(j).fas_x{f};
                handles.Region(i).Fascicle(j).fas_y_original{f} = handles.Region(i).Fascicle(j).fas_y{f};

                % in the new version, aponeurosis has its own warp (in the old version it doesnt)
                if contains(handles.ROItype, 'Hough')

                    % apply warp to aponeurosis
                    super_prev = [handles.Region(i).sup_x{f-1} handles.Region(i).sup_y{f-1}];
                    super_new = transformPointsForward(wa, super_prev);

                    deep_prev = [handles.Region(i).deep_x{f-1} handles.Region(i).deep_y{f-1}];
                    deep_new = transformPointsForward(wa, deep_prev);

                    % save
                    handles.Region(i).sup_x{f} = [1 n]';
                    handles.Region(i).sup_y{f} = super_new(:,2);
                    handles.Region(i).deep_x{f} = [1 n]';
                    handles.Region(i).deep_y{f} = deep_new(:,2);

                    % get ROI
                    ROIx = [1 1 n n 1]';
                    super_apo = handles.geofeatures(f).super_pos';
                    deep_apo = handles.geofeatures(f).deep_pos';

                    ROIy = [super_apo(1); deep_apo; super_apo([2,1])];
                    ROIy_fcor = round([super_apo(1)+thickness(1)*r; deep_apo-thickness*r; super_apo([2,1])+thickness([2,1])*r]);

                else
                    % apply warp to ROI
                    ROIpos = transformPointsForward(wf, [handles.Region(i).ROIx{f-1} handles.Region(i).ROIy{f-1}]);

                    ROIx = ROIpos(:,1);
                    ROIy = ROIpos(:,2);

                    ROIx(ROIx > handles.vidWidth) = handles.vidWidth;
                    ROIy(ROIy > handles.vidHeight) = handles.vidHeight;
                    ROIx(ROIx < 1) = 1;
                    ROIy(ROIy < 1) = 1;

                    handles.Region(i).sup_x{f}  = ROIx([1,4]);
                    handles.Region(i).sup_y{f}  = ROIy([1,4]);
                    handles.Region(i).deep_x{f} = ROIx([2,3]);
                    handles.Region(i).deep_y{f} = ROIy([2,3]);
                end

                % save ROI
                handles.Region(i).ROIx{f} = ROIx;
                handles.Region(i).ROIy{f} = ROIy;

                % calculate the length and pennation for the current frame
                handles = calc_fascicle_length_and_pennation(handles,f);

                % update the points
                fpoints = fpointsNew;

                % if drops below 100, define new points
                if length(fpoints) < 100 || ~strcmp(handles.ROItype(1:5), 'Hough')

                    % detect points
                    fpoints = detectMinEigenFeatures(I_fmasked,'FilterSize',11, 'MinQuality', 0.005);

                    if strcmp(handles.ROItype(1:5), 'Hough')
                        fpoints = fpoints.selectStrongest(300);
                    end

                    fpoints = double(fpoints.Location);
                end

                % must be in ROI
                if strcmp(handles.ROItype(1:5), 'Hough')
                    inPoints = inpolygon(fpoints(:,1),fpoints(:,2), ROIx, ROIy_fcor);
                else
                    inPoints = inpolygon(fpoints(:,1),fpoints(:,2), ROIx, ROIy);
                end

                fpoints = fpoints(inPoints,:);

                % set tracker
                setPoints(fpointTracker, fpoints);

                % update the points
                if strcmp(handles.ROItype(1:5), 'Hough')
                    apoints = apointsNew;

                    if length(apoints) < 500
                        % detect points
                        apoints = detectMinEigenFeatures(I_amasked,'FilterSize',11, 'MinQuality', 0.005);
                        apoints = double(apoints.Location);
                    end

                    m = handles.vidHeight;

                    s = handles.parms.apo.super.cut;
                    ROIys = [s(1) s(2) s(2) s(1) s(1)]' * m;

                    d = handles.parms.apo.deep.cut;
                    ROIyd = [d(1) d(2) d(2) d(1) d(1)]' * m;

                    % must be in ROI
                    dinPoints = inpolygon(apoints(:,1),apoints(:,2), ROIx, ROIyd);
                    sinPoints = inpolygon(apoints(:,1),apoints(:,2), ROIx, ROIys);
                    apoints = apoints(dinPoints | sinPoints,:);

                    % set tracker
                    setPoints(apointTracker, apoints);
                end
            end

            % save the points
            handles.points{f} = fpoints;

            if strcmp(handles.ROItype(1:5), 'Hough')
                handles.apoints{f} = apoints;
            end

        end

        frac_progress = ((f-handles.start_frame)+(get(handles.frame_slider,'Max')*(i-1))) / (get(handles.frame_slider,'Max')*length(handles.Region));
        waitbar(frac_progress,h, ['Processing frame ', num2str((f-handles.start_frame+1)), '/', num2str(get(handles.frame_slider,'Max'))])
    end

end
close(h)
handles.ProcessingTime(2) = toc(tstart);

%%%%%%% TimTrack process
function[handles] = process_all_TimTrack(hObject, eventdata, handles)

% run TimTrack on all frames
frames = (handles.start_frame+1):(handles.start_frame + handles.NumFrames-1);
numIterations = length(frames);

parms = handles.parms;
parms.extrapolation = 1;

n = handles.vidWidth;


if isfield(handles,'ImStack')
    im2 = imresize(handles.ImStack, 1/handles.imresize_fac);

    % call once to get the correct fascicle region
    auto_ultrasound(im2(:,:,handles.start_frame), parms);

    % call again a bunch of times to get estimate the total duration
    for i = 1:min([size(im2,3), 5])
        tstart = tic;
        auto_ultrasound(im2(:,:,1), parms);
        dt = toc(tstart);
    end

    est_duration = dt * numIterations;

    % prompt to optionally change processing based on computational time
    answer = 'Undefined';
    if (est_duration > 60) && handles.do_parfor.Value == 0
        answer = questdlg(['Estimated TimTrack duration: ', num2str(est_duration), ' s, would you like to use parallel pool?'], 'Type of computation', 'Yes','No','Cancel','Yes');
    end

    if strcmp(answer, 'Yes')
        handles.do_parfor.Value = 1;
    elseif strcmp(answer, 'No')
        handles.do_parfor.Value = 0;
    elseif strcmp(answer,'Cancel')
        return
    end

    if ~strcmp(answer,'Cancel')
        for i = 1:length(handles.Region)
            % TimTrack (parfor or for)
            if handles.do_parfor.Value
                % Then construct a ParforProgressbar object:
                WaitMessage = parfor_wait(numIterations,'Waitbar', true,'Title','Running TimTrack...');

                tstart = tic;
                parfor f = frames
                    geofeatures(f) = auto_ultrasound(im2(:,:,f), parms);
                    WaitMessage.Send; %update waitbar parfor
                end
                WaitMessage.Destroy(); %update waitbar parfor
                handles.ProcessingTime(1) = toc(tstart);
                %

            else
                %single thread for loop
                tstart = tic;
                hwb = waitbar(0,'','Name','Running TimTrack...');
                for f = frames

                    geofeatures(f) = auto_ultrasound(im2(:,:,f), parms);
                    waitbar((f-handles.start_frame) / numIterations, hwb, sprintf('Processing frame %d/%d', (f-handles.start_frame), numIterations)); %maybe numIterations +1 on bar, but not sure for cutting frames

                end
                close(hwb)
                handles.ProcessingTime(1) = toc(tstart);
            end

            % Adjust the parameter of geofeatures
            for kk = frames
                geofeatures(kk).fas_coef(2)     = geofeatures(kk).fas_coef(2) * handles.imresize_fac;
                geofeatures(kk).super_coef(2)   = geofeatures(kk).super_coef(2) * handles.imresize_fac;
                geofeatures(kk).deep_coef(2)    = geofeatures(kk).deep_coef(2) * handles.imresize_fac;
                geofeatures(kk).thickness       = geofeatures(kk).thickness * handles.imresize_fac;
                geofeatures(kk).faslen          = geofeatures(kk).faslen * handles.imresize_fac;

                % get vertical locations at image boundaries
                geofeatures(kk).super_pos = polyval(geofeatures(kk).super_coef, [1 n]);
                geofeatures(kk).deep_pos = polyval(geofeatures(kk).deep_coef, [1 n]);

                % calculate fascicle
                Deep_intersect_x = round((geofeatures(kk).deep_coef(2) - geofeatures(kk).fas_coef(2))   ./ (geofeatures(kk).fas_coef(1) - geofeatures(kk).deep_coef(1)));
                Super_intersect_x = round((geofeatures(kk).super_coef(2) - geofeatures(kk).fas_coef(2)) ./ (geofeatures(kk).fas_coef(1) - geofeatures(kk).super_coef(1)));
                Super_intersect_y = polyval(geofeatures(kk).super_coef, Super_intersect_x);
                Deep_intersect_y = polyval(geofeatures(kk).deep_coef, Deep_intersect_x);

                handles.Region(i).sup_x{kk} = [1 n]';
                handles.Region(i).sup_y{kk} = polyval(geofeatures(kk).super_coef, [1 n]');

                handles.Region(i).deep_x{kk} = [1 n]';
                handles.Region(i).deep_y{kk} = polyval(geofeatures(kk).deep_coef, [1 n]');

                handles.Region(i).ROIx{kk} = [1 1 n n 1]';
                handles.Region(i).ROIy{kk} = [polyval(geofeatures(kk).super_coef, 1); polyval(geofeatures(kk).deep_coef, [1 n]'); polyval(geofeatures(kk).super_coef, [n 1]')];

                handles.Region.Fascicle.fas_x{kk} = [Deep_intersect_x Super_intersect_x]';
                handles.Region.Fascicle.fas_y{kk} = [Deep_intersect_y Super_intersect_y]';

                handles.Region.Fascicle.fas_x_original{kk} = handles.Region.Fascicle.fas_x{kk};
                handles.Region.Fascicle.fas_y_original{kk} = handles.Region.Fascicle.fas_y{kk};

            end

            handles.geofeatures = geofeatures;
        end
        
        % get the first from auto_detect
        geofeatures(handles.start_frame) = handles.geofeatures(handles.start_frame);
        
        handles.geofeatures = geofeatures;
    end
end

function[handles] = estimate_variance(hObject, eventdata, handles)

geofeatures = handles.geofeatures;

% frames = handles.start_frame:handles.start_frame + handles.NumFrames-1;
% numIterations = length(frames);

x = nan(handles.NumFrames+handles.start_frame-1,5);

for f = handles.start_frame+1:handles.NumFrames+handles.start_frame-1
    x(f,1) = geofeatures(f).alpha;
    %     x(f,2) = geofeatures(f).thickness;
    x(f,2) = geofeatures(f).super_pos(1);
    x(f,3) = geofeatures(f).super_pos(2);
    x(f,4) = geofeatures(f).deep_pos(1);
    x(f,5) = geofeatures(f).deep_pos(1);
end

Wn = 1.5*handles.fc_lpf / (.5 * handles.FrameRate);
Wn(Wn>=1) = 1-1e-6;
Wn(Wn<=0) = 1e-6;
[b,a] = butter(2, Wn, 'high');

xcya = isnan(x(:,1));
x(xcya,:) = [];

y = nan(size(x));
for i = 1:size(x,2)
    y(:,i) = filtfilt(b,a,x(:,i));
end

handles.R = var(y);

function[handles] = do_state_estimation(hObject, eventdata, handles)
% hObject    handle to do_state_estimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get Qmax
if ~isnan(handles.Q)
    handles = estimate_variance(hObject, eventdata, handles);

    if handles.trackbck_chkBox.Value == 0
    %%here they are for normal forward  Optic flow tracking
    % forward state estimation
    for f = (handles.start_frame+1):length(handles.geofeatures)
        for i = 1:length(handles.Region)
            for j = 1:length(handles.Region(i).Fascicle)
                % state estimation
                handles = apo_state_estimator(handles,f,f-1);
            end
        end
    end

    % forward state estimation
    for f = (handles.start_frame+1):length(handles.geofeatures)
        for i = 1:length(handles.Region)
            for j = 1:length(handles.Region(i).Fascicle)
                % state estimation
                handles = state_estimator(handles,f,f-1);
            end
        end
    end

%     Rauch-Tung-Striebel backwards filter
    for f = (length(handles.geofeatures)-1):-1:handles.start_frame
        for i = 1:length(handles.Region)
            for j = 1:length(handles.Region(i).Fascicle)
                handles = state_smoothener(handles,f,f+1);
            end
        end
    end

    else %%% backwards estimator for backwards opticflow tracking
        %%here they are for normal forward  Optic flow tracking
        % forward state estimation
        
        for f = (handles.NumFrames+handles.start_frame-1) :-1 : (handles.start_frame)
            for i = 1:length(handles.Region)
                for j = 1:length(handles.Region(i).Fascicle)
                        %frame_geo = (f - (handles.NumFrames))+1; %old
                        %solution with no cutting including
                        frame_geo = (handles.NumFrames + handles.start_frame - 1) - (f - handles.start_frame); % Direct calculation
                        % state estimation
                        handles = apo_state_estimator(handles,f-1,f,frame_geo);
                    
                end
            end
        end
    
        % forward state estimation
        for f = (handles.NumFrames+handles.start_frame-1) :-1 :(handles.start_frame)
            for i = 1:length(handles.Region)
                for j = 1:length(handles.Region(i).Fascicle)
                   %frame_geo = -(f - (handles.NumFrames));
                    frame_geo = (handles.NumFrames + handles.start_frame - 1) - (f - handles.start_frame); % Direct calculation

                    % state estimation
                    handles = state_estimator(handles,f-1,f,frame_geo);
                end
            end
        end
    
    %     Rauch-Tung-Striebel backwards filter
        for f = handles.start_frame +1 :  handles.NumFrames+handles.start_frame-1
            for i = 1:length(handles.Region)
                for j = 1:length(handles.Region(i).Fascicle)
                    handles = state_smoothener(handles,f,f-1); 
                end
            end
        end
        
    end
    show_image(hObject,handles);
    show_data(hObject, handles);
    guidata(hObject, handles);
end

function[Q] = getQ(handles, dx)

% Optical flow error is proportional to flow
Q = handles.Q  * dx^2;

function [K] = run_kalman_filter(k)
% this assumes we already have the aposteriori state estimate (k.x_minus),
% the measurement (k.y) and the process- and measurement noise covariances (k.R and k.Qvalue)

% adjust kalman gain based on measurement variance
K.K = k.P_minus / (k.P_minus + k.R);

% check for weird gains
if (K.K < 0) || (K.K > 1)
    disp('warning: kalman gain outside 0-1 interval');

    K.K(K.K<0) = 0;
    K.K(K.K>1) = 1;
end

% update state
K.x_plus = k.x_minus + K.K * (k.y - k.x_minus);

% update variance
K.P_plus = (1-K.K) * k.P_minus;

%State estimator smoothener
function[handles] = state_smoothener(handles,frame_no,prev_frame_no)


i = 1;
j = 1;

Pcorr = handles.Region(i).Fascicle(j).fas_p{frame_no};
Ppred = handles.Region(i).Fascicle(j).fas_p_minus{prev_frame_no};

xcorr = [handles.Region(i).Fascicle(j).fas_x{frame_no}(2) handles.Region(i).Fascicle(j).alpha{frame_no}];
xpred = handles.Region(i).Fascicle(j).X_minus{prev_frame_no};
xsmooth = [handles.Region(i).Fascicle(j).fas_x{prev_frame_no}(2) handles.Region(i).Fascicle(j).alpha{prev_frame_no}];

Psmooth = ones(1,2);

for m = 1:2
    A           = Pcorr(m)/Ppred(m);
    A(isnan(A)) = 1;

    xsmooth(m)     = xcorr(m) + A*(xsmooth(m) - xpred(m));
    Psmooth(m)     = Pcorr(m) + A*(Psmooth(m) - Ppred(m))*A;
end

fasx2_smooth = xsmooth(1);
alpha_smooth = xsmooth(2);
% fasx2_smooth = xcorrm(1);

% recalc length and pennation
% fit the current aponeurosis
super_apo   = [handles.Region(i).sup_x{frame_no} handles.Region(i).sup_y{frame_no}];
deep_apo    = [handles.Region(i).deep_x{frame_no} handles.Region(i).deep_y{frame_no}];
super_coef  = polyfit(super_apo(:,1), super_apo(:,2), 1);
deep_coef   = polyfit(deep_apo(:,1), deep_apo(:,2), 1);

% get the vertical point from the estimated aponeurosis
% fasy2_smooth = super_coef(2) + fasx2_smooth*super_coef(1);
if handles.trackbck_chkBox.Value == 0
    fasy2 = handles.Region(i).Fascicle(j).fas_y{handles.start_frame}(2);
else
    fasy2 = handles.Region(i).Fascicle(j).fas_y{handles.NumFrames}(2);
end
% get the deep attachment point from the superficial point and the angle
fas_coef(1) = -tand(alpha_smooth);
fas_coef(2) =  fasy2 - fas_coef(1) * fasx2_smooth;

fasx1_smooth = (fas_coef(2) - deep_coef(2)) / (deep_coef(1) - fas_coef(1));
fasy1_smooth = deep_coef(2) + fasx1_smooth*deep_coef(1);

fasy2_smooth = fasy2;

fasx2_smooth_end = (fas_coef(2) - super_coef(2)) / (super_coef(1) - fas_coef(1));
fasy2_smooth_end = super_coef(2) + fasx2_smooth_end*super_coef(1);

handles.Region(i).Fascicle(j).fas_x{frame_no}   = [fasx1_smooth; fasx2_smooth];
handles.Region(i).Fascicle(j).fas_y{frame_no}   = [fasy1_smooth; fasy2_smooth];

handles.Region(i).Fascicle(j).fas_x_end{frame_no}   = [fasx1_smooth; fasx2_smooth_end];
handles.Region(i).Fascicle(j).fas_y_end{frame_no}   = [fasy1_smooth; fasy2_smooth_end];

handles.Region(i).Fascicle(j).alpha{frame_no}       = alpha_smooth;
handles.Region(i).Fascicle(j).fas_p{frame_no}       = Psmooth;

handles.Region(i).Fascicle(j).A{frame_no} = A;

% calculate the length and pennation for the current frame
handles = calc_fascicle_length_and_pennation(handles,frame_no);


function[handles] = apo_state_estimator(handles,frame_no,prev_frame_no, varargin)

% Check if the 'geofeatures frame' input is provided (when backwards
% trackingw as perform)
if nargin < 4
    % 'frame n' input is not provided, handle accordingly
    frame_no_geo = frame_no; % You can set a default because normal tracking 
else
    %
    frame_no_geo = varargin{1};
end
i = 1;
n = handles.vidWidth;

%% Apply warp
super_prev = [handles.Region(i).sup_x{prev_frame_no} handles.Region(i).sup_y{prev_frame_no}];
deep_prev = [handles.Region(i).deep_x{prev_frame_no} handles.Region(i).deep_y{prev_frame_no}];

w = handles.Region(i).awarp(:,:,prev_frame_no);
super_new = transformPointsForward(w, super_prev);
deep_new = transformPointsForward(w, deep_prev);

apo_prev = [super_prev; deep_prev];
apo_new = [super_new; deep_new];

apo_new_y = apo_new(:,2);
apo_prev_y = apo_prev(:,2);

apo_plus = nan(size(apo_new_y));
P_plus = nan(size(apo_new_y));

apo_y = [handles.geofeatures(frame_no_geo).super_pos'; handles.geofeatures(frame_no_geo).deep_pos'];

% check whether manual tracking exists
if isfield(handles.Region(i), 'sup_x_manual')
    if length(handles.Region(i).sup_x_manual) >= frame_no
        if ~isempty(handles.Region(i).sup_x_manual{frame_no})
            apo_y2 = [handles.Region(i).sup_y_manual{frame_no}; handles.Region(i).deep_y_manual{frame_no}];
        end
    end
end

Rs = handles.R(2:end) * .01;
if handles.trackbck_chkBox.Value == 0
    handles.Region(i).apo_p{handles.start_frame} = Rs; %why start frame without+1?
else
    handles.Region(i).apo_p{handles.NumFrames + handles.start_frame-1} = Rs; %why end frame?
end

% loop over points (4)
for kk = 1:numel(apo_new_y)
    dapo = abs(apo_new_y(kk) - apo_prev_y(kk));

    % A priori state estimate
    x_minus = apo_new_y(kk);

    % State estimation superficial aponeurosis attachment
    % previous estimate covariance
    P_prev = handles.Region(i).apo_p{prev_frame_no}(kk);

    % get the process noise and measurement noise covariance
    s.Q = getQ(handles, dapo);
    R(1) = Rs(kk);

    % a priori estimate from optical flow
    s.x_minus = x_minus;

    % 'measurement' is from TimTrack
    y(1) = apo_y(kk);

    if exist('apo_y2','var')
        y(2) = apo_y2(kk);
        R(2) = handles.R_manual;
    end

    % previous state covariance
    s.P_prev = P_prev;

    % a posteriori variance estimate
    s.P_minus = s.P_prev + s.Q;

    % run kalman filter
    for ii = 1:length(y)

        if ii == 2
            % update
            s.x_minus = S.x_plus;
            s.P_minus = S.P_plus;
        end

        s.y = y(ii);
        s.R = R(ii);

        S = run_kalman_filter(s);

    end

    apo_gain(kk) = S.K;
    apo_plus(kk) = S.x_plus;
    P_plus(kk) = S.P_plus;

    if isinf(handles.Q)
        apo_plus(kk) = apo_y(kk);
        P_plus(kk) = s.R;
    end

end

handles.Region(i).apo_p{frame_no} = P_plus;

handles.Region(i).apo_gain{frame_no} = apo_gain';

handles.Region(i).sup_x{frame_no} = [1 n]';
handles.Region(i).sup_y{frame_no} = apo_plus(1:2);
handles.Region(i).deep_x{frame_no} = [1 n]';
handles.Region(i).deep_y{frame_no} = apo_plus(3:4);

function[handles] = state_estimator(handles,frame_no,prev_frame_no,varargin)
% Check if the 'geofeatures frame' input is provided (when backwards
% trackingw as perform)
if nargin < 4
    % 'frame n' input is not provided, handle accordingly
    frame_no_geo = frame_no; % You can set a default because normal tracking 
else
    %
    frame_no_geo = varargin{1};
end

i = 1; j = 1;

% if the second frame, we need to calculate the first frame
if frame_no == (handles.start_frame + 1) && handles.trackbck_chkBox.Value == 0% && ~isfield(handles.Region,'fas_ang')

    alpha0 = nan(1,handles.NS);

    for k = 1:handles.NS % number of starting frames
        alpha0(k) = atan2d(-diff(handles.Region(i).Fascicle(j).fas_y_original{k+handles.start_frame-1}), diff(handles.Region(i).Fascicle(j).fas_x_original{k+handles.start_frame-1}));
%         alpha0(k) = handles.geofeatures(k+handles.start_frame-1).alpha;
    end

    if isempty(handles.Region.Fascicle.fas_x_original{handles.start_frame})
        handles.Region(i).Fascicle(j).X_plus{handles.start_frame} = [handles.Region(i).Fascicle(j).fas_x{handles.start_frame}(2) mean(alpha0)];
    else
        handles.Region(i).Fascicle(j).X_plus{handles.start_frame} = [handles.Region(i).Fascicle(j).fas_x_original{handles.start_frame}(2) mean(alpha0)];
    end
    handles.Region(i).Fascicle(j).fas_p{handles.start_frame} = [0 var(alpha0)];

    % if manual is available for first frame, overrule
    if isfield(handles.Region(i).Fascicle(j), 'fas_x_manual')
        if ~isempty(handles.Region(i).Fascicle(j).fas_x_manual)
            if ~isempty(handles.Region(i).Fascicle(j).fas_x_manual{handles.start_frame})
                handles.Region(i).Fascicle(j).X_plus{handles.start_frame} = [handles.Region(i).Fascicle(j).fas_x_manual{handles.start_frame}(2) handles.Region(i).fas_ang_manual(handles.start_frame)];
                handles.Region(i).Fascicle(j).fas_p{handles.start_frame} = [0 0];
            end
        end
    end

    % a priori is the same as a positeriori
    handles.Region(i).Fascicle(j).fas_p_minus{handles.start_frame} = handles.Region(i).Fascicle(j).fas_p{handles.start_frame};
    handles.Region(i).Fascicle(j).X_minus{handles.start_frame} = handles.Region(i).Fascicle(j).X_plus{handles.start_frame};

    handles = update_Fascicle(handles,handles.start_frame);

    %in case of backwards tracking here
elseif frame_no == (handles.NumFrames -1 + handles.start_frame - 1) && handles.trackbck_chkBox.Value == 1
    alpha0 = nan(1,handles.NS);

    for k = 1:handles.NS % number of starting frames
        alpha0(k) = atan2d(-diff(handles.Region(i).Fascicle(j).fas_y_original{frame_no-k}), diff(handles.Region(i).Fascicle(j).fas_x_original{frame_no-k}));
%         alpha0(k) = handles.geofeatures(k+handles.start_frame-1).alpha;
    end

    if isempty(handles.Region.Fascicle.fas_x_original{frame_no+1})
        handles.Region(i).Fascicle(j).X_plus{frame_no+1} = [handles.Region(i).Fascicle(j).fas_x{frame_no+1}(2) mean(alpha0)];
    else
        handles.Region(i).Fascicle(j).X_plus{frame_no+1} = [handles.Region(i).Fascicle(j).fas_x_original{frame_no+1}(2) mean(alpha0)];
    end
    handles.Region(i).Fascicle(j).fas_p{frame_no+1} = [0 var(alpha0)];

    % if manual is available for first frame, overrule
    if isfield(handles.Region(i).Fascicle(j), 'fas_x_manual')
        if ~isempty(handles.Region(i).Fascicle(j).fas_x_manual)
            if ~isempty(handles.Region(i).Fascicle(j).fas_x_manual{frame_no+1})
                handles.Region(i).Fascicle(j).X_plus{frame_no+1} = [handles.Region(i).Fascicle(j).fas_x_manual{frame_no+1}(2) handles.Region(i).fas_ang_manual(frame_no+1)];
                handles.Region(i).Fascicle(j).fas_p{frame_no+1} = [0 0];
            end
        end
    end

    % a priori is the same as a positeriori
    handles.Region(i).Fascicle(j).fas_p_minus{frame_no+1} = handles.Region(i).Fascicle(j).fas_p{frame_no};
    handles.Region(i).Fascicle(j).X_minus{frame_no+1} = handles.Region(i).Fascicle(j).X_plus{frame_no};

    handles = update_Fascicle(handles,frame_no+1);

%     if isfield(handles.Region,'fas_ang')
% 
%         handles.Region(1).Fascicle(1).alpha{handles.start_frame} = handles.Region(1).fas_ang(handles.start_frame);
%         handles.Region(i).Fascicle(j).X_plus{handles.start_frame} = [handles.Region(i).Fascicle(j).fas_x{handles.start_frame}(2) handles.Region(1).fas_ang(handles.start_frame)];
%         handles.Region(i).Fascicle(j).fas_p{handles.start_frame} = [0 var(handles.Region(1).fas_ang(handles.start_frame))];
%         handles.Region(i).Fascicle(j).fas_p_minus{handles.start_frame} = handles.Region(i).Fascicle(j).fas_p{handles.start_frame};
%         handles.Region(i).Fascicle(j).X_minus{handles.start_frame} = handles.Region(i).Fascicle(j).X_plus{handles.start_frame};
% 
%     end

end

%% Apply warp
fas_prev = [handles.Region(i).Fascicle(j).fas_x{prev_frame_no} handles.Region(i).Fascicle(j).fas_y{prev_frame_no}];
alpha_prev = handles.Region(i).Fascicle(j).alpha{prev_frame_no};

w = handles.Region(i).warp(:,:,prev_frame_no);
fas_new = transformPointsForward(w, fas_prev);

% Estimate the change in fascicle angle from the change in points
dalpha = abs(atan2d(diff(fas_new(:,2)), diff(fas_new(:,1)))) - abs(atan2d(diff(fas_prev(:,2)), diff(fas_prev(:,1))));
alpha_new = alpha_prev + dalpha;

% A priori state estimate
x_minus = [fas_new(2,1) alpha_new];

%% State estimation superficial aponeurosis attachment
% previous estimate covariance
P_prev = handles.Region(i).Fascicle(j).fas_p{prev_frame_no};

% get the process noise and measurement noise covariance
dx = sqrt((fas_new(2,1)-fas_prev(2,1)).^2 + (fas_new(2,2)-fas_prev(2,2)).^2);
s.Q = getQ(handles, dx);
R(1) = handles.X;

% a priori estimate from optical flow
s.x_minus = x_minus(1);

% 'measurement', here is the first value
y(1) = handles.Region(i).Fascicle(j).fas_x{handles.start_frame}(2);

% if there is a manual estimate, add a second measurement
if isfield(handles.Region(i).Fascicle(j), 'fas_x_manual')
    if length(handles.Region(i).Fascicle(j).fas_x_manual) >= frame_no
        if ~isempty(handles.Region(i).Fascicle(j).fas_x_manual{frame_no})
            y(2) = handles.Region(i).Fascicle(j).fas_x_manual{frame_no}(2);
            R(2) = handles.R_manual;
        end
    end
end

% previous state covariance
s.P_prev = P_prev(1);

% a posteriori variance estimate
s.P_minus = s.P_prev + s.Q;

% run kalman filter
for ii = 1:length(y)

    if ii == 2
        % update
        s.x_minus = S.x_plus;
        s.P_minus = S.P_plus;
    end

    s.y = y(ii);
    s.R = R(ii);

    S = run_kalman_filter(s);

end

% ascribe
fasx2_plus = S.x_plus;
fasx2_minus = s.x_minus;
supP_plus = S.P_plus;
supP_minus = s.P_minus;

if isinf(handles.Q)
    fasx2_plus = s.y;
    supP_plus = s.R;
elseif isinf(handles.X)
    fasx2_plus = s.x_minus;
    supP_plus = s.Q;
end

%% Fascicle angle estimate
% get the process noise and measurement noise covariance
dx = abs(dalpha);
dx(dx<0.005) = 0;
f.Q = getQ(handles, dx);
R(1) = handles.R(1);

% apriori estimate from optical flow
f.x_minus = x_minus(2);

% measurement from Hough transform
y(1) = handles.geofeatures(frame_no_geo).alpha;

% if there is a manual estimate, add a second measurement
if isfield(handles.Region(i).Fascicle(j), 'fas_x_manual')
    if length(handles.Region(i).Fascicle(j).fas_x_manual) >= frame_no
        if ~isempty(handles.Region(i).Fascicle(j).fas_x_manual{frame_no})
            y(2) = handles.Region(i).fas_ang_manual(frame_no);
            R(2) = handles.R_manual / (pi*handles.Region(i).fas_length_manual(frame_no)*handles.vidHeight/handles.ID) * 180;
        end
    end
end

% previous state covariance
f.P_prev = P_prev(2);

% a posteriori variance estimate
f.P_minus = f.P_prev + f.Q;

% run kalman filter
for ii = 1:length(y)

    if ii == 2
        % update
        f.x_minus = F.x_plus;
        f.P_minus = F.P_plus;
    end

    f.y = y(ii);
    f.R = R(ii);

    F = run_kalman_filter(f);

end

% ascribe
alpha_minus = f.x_minus;
alpha_plus = F.x_plus;
Kgain = F.K;
fasP_plus = F.P_plus;
fasP_minus = f.P_minus;

if isinf(handles.Q)
    alpha_plus = f.y;
    Kgain = 1;
    fasP_plus = f.R;
end

% state estimate
handles.Region(i).Fascicle(j).X_plus{frame_no}  = [fasx2_plus alpha_plus];
handles.Region(i).Fascicle(j).X_minus{frame_no} = [fasx2_minus alpha_minus];

% state covariance
handles.Region(i).Fascicle(j).fas_p{frame_no} = [supP_plus fasP_plus];
handles.Region(i).Fascicle(j).fas_p_minus{frame_no} = [supP_minus fasP_minus];

% kalman xshiftcor for fascicle
handles.Region(i).Fascicle(j).K(frame_no) = Kgain;

% update fascicle
handles = update_Fascicle(handles,frame_no);

function[handles] = update_Fascicle(handles,frame_no)
% gets fascicle tracking estimates from the state and tracked aponeuroses
i = 1;
j = 1;

% if frame_no == 5
%     keyboard
% end

% get the state
fasx2_plus = handles.Region(i).Fascicle(j).X_plus{frame_no}(1);
alpha_plus = handles.Region(i).Fascicle(j).X_plus{frame_no}(2);

% fit the current aponeurosis
super_apo   = [handles.Region(i).sup_x{frame_no} handles.Region(i).sup_y{frame_no}];
deep_apo    = [handles.Region(i).deep_x{frame_no} handles.Region(i).deep_y{frame_no}];
super_coef  = polyfit(super_apo(:,1), super_apo(:,2), 1);
deep_coef   = polyfit(deep_apo(:,1), deep_apo(:,2), 1);

% get the vertical point from the estimated aponeurosis
% fasy2_plus = super_coef(2) + fasx2_plus*super_coef(1);

% vertical location is fixed
fasy2 = handles.Region(i).Fascicle(j).fas_y{handles.start_frame}(2);

% get the deep attachment point from the superficial point and the angle
fas_coef(1) = -tand(alpha_plus);
fas_coef(2) =  fasy2 - fas_coef(1) * fasx2_plus;

% deep
fasx1_end = (fas_coef(2) - deep_coef(2)) / (deep_coef(1) - fas_coef(1));
fasy1_end = deep_coef(2) + fasx1_end*deep_coef(1);

% superficial
fasx2_end = (fas_coef(2) - super_coef(2)) / (super_coef(1) - fas_coef(1));
fasy2_end = super_coef(2) + fasx2_end*super_coef(1);

%% update
% state and dependent variables
handles.Region(i).Fascicle(j).fas_x{frame_no}   = [fasx1_end; fasx2_plus];
handles.Region(i).Fascicle(j).fas_y{frame_no}   = [fasy1_end; fasy2];

handles.Region(i).Fascicle(j).fas_x_end{frame_no}   = [fasx1_end; fasx2_end];
handles.Region(i).Fascicle(j).fas_y_end{frame_no}   = [fasy1_end; fasy2_end];

handles.Region(i).Fascicle(j).alpha{frame_no}   = alpha_plus;

handles = calc_fascicle_length_and_pennation(handles,frame_no);


function[handles] = calc_fascicle_length_and_pennation(handles,frame_no)
i = 1;
j = 1;
%if to check and init because in case of cut before/after then load a
%fascicle, it crashes in plotting because Time is long the entire cut video
%but not the featuere to plot (y data)
if ~isfield(handles.Region,'fas_length') || ~isfield(handles.Region,'fas_pen') || length(handles.Region(i).fas_length) < handles.NumFrames
    handles.Region(i).fas_pen = nan((handles.NumFrames+handles.start_frame-1),j);
    handles.Region(i).fas_ang = nan((handles.NumFrames+handles.start_frame-1),j);
    handles.Region(i).fas_length = nan((handles.NumFrames+handles.start_frame-1),j);
end
% fit the current aponeurosis
% ROI         = [handles.Region(i).ROIx{frame_no} handles.Region(i).ROIy{frame_no}];
% deep_apo    = ROI([2,3],:);

deep_apo =  [handles.Region(i).deep_x{frame_no} handles.Region(i).deep_y{frame_no}];

% deep aponeurosis angle
gamma = atan2d(-diff(deep_apo(:,2)), diff(deep_apo(:,1)));

% calculate the length and pennation for the current frame
%if length(handles.Region(i).Fascicle(j).alpha) >= frame_no
if isfield(handles.Region(i).Fascicle(j),'alpha') && length(handles.Region(i).Fascicle(j).alpha) >= frame_no && ~isempty(handles.Region(i).Fascicle(j).alpha{frame_no})
    handles.Region(i).fas_pen(frame_no,j) = handles.Region(i).Fascicle(j).alpha{frame_no} - gamma;
    handles.Region(i).fas_ang(frame_no,j) = handles.Region(i).Fascicle(j).alpha{frame_no};
else
    handles.Region(i).fas_pen(frame_no,j) = atan2d(-diff(handles.Region(i).Fascicle(j).fas_y{frame_no}), diff(handles.Region(i).Fascicle(j).fas_x{frame_no}));
end

if isfield(handles.Region(i).Fascicle(j), 'fas_y_end') && ~isempty(handles.Region(i).Fascicle(j).fas_y_end{frame_no})

    fasx = handles.Region(i).Fascicle(j).fas_x_end{frame_no};
    fasy = handles.Region(i).Fascicle(j).fas_y_end{frame_no};
else
    fasx = handles.Region(i).Fascicle(j).fas_x{frame_no};
    fasy = handles.Region(i).Fascicle(j).fas_y{frame_no};
end

handles.Region(i).fas_length(frame_no,j) = (handles.ID/handles.vidHeight)*sqrt(diff(fasx).^2 + diff(fasy).^2);

if isfield(handles.Region(i).Fascicle(j), 'fas_x_manual')
    if length(handles.Region(i).Fascicle(j).fas_x_manual) >= frame_no
        if ~isempty(handles.Region(i).Fascicle(j).fas_x_manual{frame_no})
            fasx = handles.Region(i).Fascicle(j).fas_x_manual{frame_no};
            fasy = handles.Region(i).Fascicle(j).fas_y_manual{frame_no};

            handles.Region(i).fas_length_manual(frame_no,j) = (handles.ID/handles.vidHeight)*sqrt(diff(fasx).^2 + diff(fasy).^2);
            handles.Region(i).fas_ang_manual(frame_no,j) = atan2d(-diff(fasy), diff(fasx));
            handles.Region(i).fas_pen_manual(frame_no,j) = atan2d(-diff(fasy), diff(fasx)) - gamma;
        end
    end
end


% --- Executes on button press in Auto_Detect.
function [handles] = Auto_Detect_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_Detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
% initialize
N = handles.NumFrames + handles.start_frame - 1;
handles.Region(1).fas_length    = nan(N,1);
handles.Region(1).fas_pen       = nan(N,1);
handles.Region(1).fas_length_manual    = nan(N,1);
handles.Region(1).fas_pen_manual       = nan(N,1);

w = handles.vidWidth;

if ~isvalid(handles.S)
    handles.S = images.roi.Rectangle(handles.axes1,'position', [1 handles.parms.apo.super.cut(1)*handles.vidHeight w diff(handles.parms.apo.super.cut)*handles.vidHeight],'color','blue');
end

if ~isvalid(handles.D)
    handles.D = images.roi.Rectangle(handles.axes1,'position', [1 handles.parms.apo.deep.cut(1)*handles.vidHeight w diff(handles.parms.apo.deep.cut)*handles.vidHeight],'color','green');
end

%% Aponeurosis detection
axes(handles.axes1);

handles.parms.apo.super.cut = [handles.S.Position(2) handles.S.Position(2)+handles.S.Position(4)] / handles.vidHeight;
handles.parms.apo.deep.cut = [handles.D.Position(2) handles.D.Position(2)+handles.D.Position(4)] / handles.vidHeight;

set(handles.S, 'EdgeAlpha',0,'FaceAlpha',0.1,'InteractionsAllowed','none')
set(handles.D, 'EdgeAlpha',0,'FaceAlpha',0.1,'InteractionsAllowed','none')

% max range
% parms.fas.range = 90 - [-90 89];

% find the first frame
frame_no = handles.start_frame + round(get(handles.frame_slider,'Value')) - 1;

% % detect orientation
data = imresize(handles.ImStack(:,:,frame_no), 1/handles.imresize_fac);

% run TimTrack
handles.parms.fas.redo_ROI = 1;
[geofeatures, ~, parms] = auto_ultrasound(data, handles.parms);

% save parms
parms.fas.redo_ROI = 0;
handles.parms = parms;

% scale
geofeatures.thickness = geofeatures.thickness * handles.imresize_fac;
geofeatures.super_coef(2) = geofeatures.super_coef(2)     .* [handles.imresize_fac];
geofeatures.deep_coef(2) = geofeatures.deep_coef(2)       .* [handles.imresize_fac];
geofeatures.fas_coef(2) = geofeatures.fas_coef(2)         .* [handles.imresize_fac];

n = handles.vidWidth;
geofeatures.super_pos = polyval(geofeatures.super_coef, [1 n]);
geofeatures.deep_pos = polyval(geofeatures.deep_coef, [1 n]);

% save geofeatures
handles.geofeatures = geofeatures;

Deep_intersect_x = round((geofeatures.deep_coef(2) - geofeatures.fas_coef(2))   ./ (geofeatures.fas_coef(1) - geofeatures.deep_coef(1)));
Super_intersect_x = round((geofeatures.super_coef(2) - geofeatures.fas_coef(2)) ./ (geofeatures.fas_coef(1) - geofeatures.super_coef(1)));
Super_intersect_y = polyval(geofeatures.super_coef, Super_intersect_x);
Deep_intersect_y = polyval(geofeatures.deep_coef, Deep_intersect_x);

i = 1;
handles.Region(i).sup_x{frame_no} = [1 n]';
handles.Region(i).sup_y{frame_no} = polyval(geofeatures.super_coef, [1 n]');

handles.Region(i).deep_x{frame_no} = [1 n]';
handles.Region(i).deep_y{frame_no} = polyval(geofeatures.deep_coef, [1 n]');

handles.Region(i).ROIx{frame_no} = [1 1 n n 1]';
handles.Region(i).ROIy{frame_no} = [polyval(geofeatures.super_coef, 1); polyval(geofeatures.deep_coef, [1 n]'); polyval(geofeatures.super_coef, [n 1]')];

handles.Region(i).Fascicle.fas_x{frame_no} = [Deep_intersect_x Super_intersect_x]';
handles.Region(i).Fascicle.fas_y{frame_no} = [Deep_intersect_y Super_intersect_y]';

handles.Region.Fascicle.fas_x_original{frame_no} = handles.Region.Fascicle.fas_x{1};
handles.Region.Fascicle.fas_y_original{frame_no} = handles.Region.Fascicle.fas_y{1};

[handles] = calc_fascicle_length_and_pennation(handles,frame_no);

% Update handles structure
guidata(hObject, handles);

show_image(hObject,handles);
show_data(hObject, handles)
catch
    warndlg('No video loaded')
end

function [handles] = extract_estimates(hObject, eventdata, handles)

i = 1;
j = 1;
frame_no = handles.start_frame + round(get(handles.frame_slider,'Value')) - 1;
d = round(size(handles.ImStack,2)/2);

handles.Region(i).sup_x_manual{frame_no} = handles.s.Position(:,1)-d;
handles.Region(i).sup_y_manual{frame_no} = handles.s.Position(:,2);

handles.Region(i).deep_x_manual{frame_no} = handles.d.Position(:,1)-d;
handles.Region(i).deep_y_manual{frame_no} = handles.d.Position(:,2);

handles.Region(i).ROIx_manual{frame_no} = [handles.s.Position(1,1); handles.d.Position([1 2],1); handles.s.Position([2 1],1)]-d;
handles.Region(i).ROIy_manual{frame_no} = [handles.s.Position(1,2); handles.d.Position([1 2],2); handles.s.Position([2 1],2)];

handles.Region(i).Fascicle(j).fas_x_manual{frame_no} = handles.h.Position(:,1)-d;
handles.Region(i).Fascicle(j).fas_y_manual{frame_no} = handles.h.Position(:,2);

handles = calc_fascicle_length_and_pennation(handles,frame_no);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in do_parfor.
function do_parfor_Callback(hObject, eventdata, handles)
% hObject    handle to do_parfor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of do_parfor
handles.do_parfor = get(hObject,'Value');

% --- Executes on button press in flipimage.
function flipimage_Callback(hObject, eventdata, handles)
% hObject    handle to flipimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flipimage
handles = do_flip(hObject, eventdata, handles);
%end
guidata(hObject, handles);
show_image(hObject,handles);


% ----- Function to perform flipping
function [handles] = do_flip(hObject, eventdata, handles)
% hObject    handle to do_state_estimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles.flip = ~handles.flip;%change value of flip

if isfield(handles,"Region")
    updateX = @(fas_x) flip(handles.vidWidth - fas_x); %only here we need correction as axis starts from 1 (plotting)

    for i = 1:numel(handles.Region)
        %adjust ROI coordinates and logic mask image
        handles.Region(i).ROIx = cellfun(updateX, handles.Region(i).ROIx, 'UniformOutput', false);
        handles.Region(i).ROIy = cellfun(@flip, handles.Region(i).ROIy, 'UniformOutput', false);

        if isfield(handles.Region(i),"ROI")
            handles.Region(i).ROI = cellfun(@fliplr, handles.Region(i).ROI, 'UniformOutput', false);
        end
        if isfield(handles.Region(i),"deep_x") && isfield(handles.Region(i),"deep_y")
            handles.Region(i).sup_y = cellfun(@flip, handles.Region(i).sup_y , 'UniformOutput', false);
            handles.Region(i).deep_y = cellfun(@flip, handles.Region(i).deep_y , 'UniformOutput', false);

        end
        %adjust each fascicle's pts
        for j = 1:numel(handles.Region(i).Fascicle)
            handles.Region(i).Fascicle(j).fas_x = cellfun(updateX, handles.Region(i).Fascicle(j).fas_x, 'UniformOutput', false);
            handles.Region(i).Fascicle(j).fas_y = cellfun(@flip, handles.Region(i).Fascicle(j).fas_y, 'UniformOutput', false);
            if isfield(handles.Region(i).Fascicle(j),'fas_x_end') %if estimator ran
                handles.Region(i).Fascicle(j).fas_x_end = cellfun(updateX, handles.Region(i).Fascicle(j).fas_x_end, 'UniformOutput', false);
                handles.Region(i).Fascicle(j).fas_y_end = cellfun(@flip, handles.Region(i).Fascicle(j).fas_y_end, 'UniformOutput', false);
            end
        end
    end

end

% If statement not necessary, if tick flip else flip back, so everytime flipimage
% changes which depends on the callback, flip the image
%if handles.flipimage

if isfield(handles, 'ImStack')
    handles.ImStack = flip(handles.ImStack, 2);
end



% --- Function to check whether ParallelToolbox exists and run it
function chkParallelToolBox()
% First checks if Parallel Computing Toolbox exists and then
% Returns num of workers available for executing parallel computing

available_toolboxes = ver;
isexist_ParallelToolBox = false;

% Check if Parallel Computing Toolbox exists
for index = 1:size(available_toolboxes, 2)
    if convertCharsToStrings(available_toolboxes(index).Name) == "Parallel Computing Toolbox"
        isexist_ParallelToolBox = true;
        break;
    end
end

% If Parallel Computing Toolbox is available, start it if not already running
if isexist_ParallelToolBox
    if isempty(gcp('nocreate')) %if parallel is not running yet, start it
        myCluster = parcluster('Processes');

        % Function to start the Parallel Computing Toolbox pool
        disp('Parallel Computing Toolbox exists and it will be opened...');
        % Open parallel pool with no idle timeout
        parpool(myCluster, 'IdleTimeout', Inf);
    else
        disp('Parpool is already running!')
    end

end


% --- Executes when user attempts to close figure1 (i.e., using the close icon
% to close the GUI)
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if Parallel Computing Toolbox is running and shut it down
delete(gcp('nocreate'));
% Hint: delete(hObject) closes the figure
delete(hObject);


% --------------------------------------------------------------------
function spatial_cal_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'ImStack')

    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'));

    % select two point for calculating the calibration mmperpx
    set(handles.axes1);
    [~,distanceInPixels]=ginputYellow(2);

    %create a simple dlg to type the real value
    %Ask the user for the real-world distance.
    userPrompt = {'Enter real distance in mm'};
    dialogTitle = 'Calibration';
    def = {''};
    answer = inputdlg(userPrompt, dialogTitle, 1, def);

    while isnan(str2double(answer{1}))  %check if it's non numeric
        answer = inputdlg(userPrompt, dialogTitle, 1, def);
    end
    %get the answer
    dist_mm = str2double(answer{1});

    %calculate mmperpx factor
    calibration_value = dist_mm /  abs(round(diff(distanceInPixels)));
    calibration_value = round(calibration_value,3); %round otherwise the conversion crash in the calculation (don't ask why)
    %update handles and GUI
    set(handles.ImDepthEdit,"String",string(calibration_value)); %this should automatically update the plots with Fascicle data
    ImDepthEdit_Callback(handles.ImDepthEdit, [], handles); % Manually call the callback function because it doesn't do automatically
end

% --- Executes on button press in Process_folder.
function Process_folder_Callback(hObject, eventdata, handles)
% hObject    handle to Process_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dir_data = uigetdir(cd,'Select folder with video(s)');

if dir_data == 0 %no folder selected, just return
    return
end
files = dir(dir_data);

%get video format list
file_formats = VideoReader.getFileFormats;
video_formats = cell(numel(file_formats),1);

for i = 1 : numel(file_formats)
    video_formats{i} = ['.' file_formats(i).get.Extension];
end

ind_toRemove = [];
for n_file = 1 : numel(files)
    [~,~,ext] = fileparts(files(n_file).name);
    %save indexes of thefile that are not videos
    if sum(strcmp(video_formats,ext)) == 0
        ind_toRemove = [ind_toRemove n_file];
    end
end

files(ind_toRemove) = []; %keep videos
% clearvars -except files eventdata hObject handles subFolders


for k = 1:numel(files) %foreach file

    handles.fname = files(k).name;
    handles.pname = [files(k).folder,'/'];
    handles.start_frame = 1;
    cd(handles.pname)
    handles.movObj = VideoReader([handles.pname handles.fname]);
    mb = waitbar(0,'Loading Video....');

    % get info
    handles.vidHeight = handles.movObj.Height;
    handles.vidWidth = handles.movObj.Width;
    handles.NumFrames = handles.movObj.NumFrames;
    handles.FrameRate = handles.movObj.FrameRate;

    i=1; %simple counter for frames reading

    % start clean
    %if isfield(handles,'ImTrack')
    %    handles = rmfield(handles, 'ImTrack');
    %end
    if isfield(handles,'ImStack')
        handles = rmfield(handles, 'ImStack');
    end
    if isfield(handles,'ImStackOr')
        handles = rmfield(handles, 'ImStackOr');
    end

    handles.ImStack = zeros(handles.vidHeight, handles.vidWidth, handles.NumFrames,'uint8');

    while hasFrame(handles.movObj)
        waitbar(handles.movObj.CurrentTime/handles.movObj.Duration,mb)
        if regexp(handles.movObj.VideoFormat,'RGB')
            handles.ImStack(:,:,i) = im2gray(readFrame(handles.movObj));
        else
            handles.ImStack(:,:,i) = readFrame(handles.movObj);
        end
        handles.ImBrightness(i) = mean(handles.ImStack(:,:,i),'all');
        i=i+1;
    end

    waitbar(1,mb)
    close(mb)

    cd(handles.pname)
    % check whether a mat file exists in the location with the same name with
    % setting of the video (ImageDepth)
    [path,name,~] = fileparts([handles.pname handles.fname]); %more elegant, people may not have necessarly mp4
    if exist([path '/' name '.mat'],"file")
        TVD = load([path '/' name '.mat']);
        if isfield(TVD.TVDdata,'cmPerPixY') %check whether the field exists and update scalar
            %ImageDepth = round(TVD.TVDdata.cmPerPixY*10,3); %round to 3 digits
            ImageDepth = round(TVD.TVDdata.Height * TVD.TVDdata.cmPerPixY,3)*10;
            handles.ID = ImageDepth;
            set(handles.ImDepthEdit,'String',num2str(ImageDepth));

            % Update handles structure
            guidata(hObject, handles);
        end
        if isfield(TVD.TVDdata,'Time') %check if also time exists as Telemed has uncostant framerate
            handles.Time = TVD.TVDdata.Time(1:end); %note the last timestamp is repeated in Echo Wave II recordings
        end
    else
        handles.Time = (double(1/handles.FrameRate):double(1/handles.FrameRate):double(handles.NumFrames/handles.FrameRate))';
    end

    % set the string in the frame_number box to the current frame value (1)
    set(handles.frame_number,'String',num2str(1))

    % set the limits on the slider - use number of frames to set maximum (min =
    % 1)
    set(handles.filename,'String',[handles.pname handles.fname])
    set(handles.frame_slider,'Min',1);
    set(handles.frame_slider,'Max',handles.NumFrames);
    set(handles.frame_slider,'Value',1);
    set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 10/handles.NumFrames]);
    set(handles.frame_rate,'String',handles.FrameRate(1))
    set(handles.vid_width,'String',handles.vidWidth(1))
    set(handles.vid_height,'String',handles.vidHeight(1))

    % update the image axes using show_image function (bottom)
    handles.ImStackOr = handles.ImStack;

    % crop
    handles = AutoCrop_Callback(hObject, eventdata, handles);

    %flip every time each video if this is what the user set
    if handles.flipimage.Value == 1%check based on flip tick box value
        handles = do_flip(hObject, eventdata, handles);
    end
    % show the new video
    handles = show_image(hObject,handles);

    %load fascicle automatically if exists
    if exist([path '/Fas_Data/Fas_' name '.mat'],'file')
        %go to last frame
        %set(handles.frame_slider,'Value',handles.NumFrames);
        %set(handles.frame_number,'String',handles.NumFrames);
        %set(handles.frame_)
        % load
        handles = menu_load_fascicle_Callback(hObject, eventdata, handles,[path '/Fas_Data/Fas_' name '.mat']);
    end

     
  
    %create a subfolder where to save results and tracked videos
    if ~exist([handles.pname 'Tracked/'], 'dir')
        mkdir([handles.pname 'Tracked/'])
    end
    handles.pname = [handles.pname 'Tracked/'];

    % process all based on what the ROI type is
    handles = process_all_Callback(hObject, eventdata, handles);

    % save
    save_video_Callback(hObject, eventdata, handles)
    Save_As_Mat_Callback(hObject, eventdata, handles)

end

function freq_lpf_Callback(hObject, eventdata, handles)
% hObject    handle to freq_lpf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_lpf as text
%        str2double(get(hObject,'String')) returns contents of freq_lpf as a double

handles.fc_lpf = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% If we have estimates, run state estimation
if isfield(handles, 'Region')
    do_state_estimation(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function freq_lpf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_lpf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.fc_lpf = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

function X_value_Callback(hObject, eventdata, handles)
% hObject    handle to X_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of X_value as text
%        str2double(get(hObject,'String')) returns contents of X_value as a double

handles.X = str2double(get(hObject,'String'));

% If we have estimates, run state estimation
if isfield(handles, 'Region')
    handles = do_state_estimation(hObject, eventdata, handles);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function X_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.X = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

function resize_fac_Callback(hObject, eventdata, handles)
% hObject    handle to resize_fac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resize_fac as text
%        str2double(get(hObject,'String')) returns contents of resize_fac as a double

handles.imresize_fac = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function resize_fac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resize_fac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.imresize_fac = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

function Qvalue_Callback(hObject, eventdata, handles)
% hObject    handle to Qvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Qvalue as text
%        str2double(get(hObject,'String')) returns contents of Qvalue as a double

handles.Q = abs(str2double(get(hObject,'String')));

% If we have estimates ALL frames, run state estimation otherwise just
% update Q value
if isfield(handles, 'Region')
    if sum(~cellfun(@isempty, handles.Region.Fascicle.fas_x, 'UniformOutput', true)) >= handles.NumFrames
        handles = do_state_estimation(hObject, eventdata, handles);
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Qvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Qvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.Q = abs(str2double(get(hObject,'String')));

% Update handles structure
guidata(hObject, handles);

function Nstat_Callback(hObject, eventdata, handles)
% hObject    handle to Nstat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nstat as text
%        str2double(get(hObject,'String')) returns contents of Nstat as a double

handles.NS = str2double(get(hObject,'String'));

if isfield(handles, 'Region')
    handles = do_state_estimation(hObject, eventdata, handles);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Nstat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nstat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.NS = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function manu_set_block_size_Callback(hObject, eventdata, handles)
% hObject    handle to manu_set_block_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = set_block_size(handles.BlockSize);

if sum(tmp ~= handles.BlockSize) ~= 0 %if the Block changed, then check and run Opticflow

    handles.BlockSize = tmp; % update blocksize according to the new values
    %re-run Optic flow automatically only if all frames were already
    %tracked
    if isfield(handles,'Region')

        for i = 1:length(handles.Region)

            %check whether all frames have been already tracked with
            %opticflow, if yes re-run it with the new block size
            if sum(~cellfun(@isempty, handles.Region.Fascicle.fas_x, 'UniformOutput', true)) >= handles.NumFrames
                %if size(handles.Region(i).Fascicle.analysed_frames,2) > 0 %double check this
                
                if handles.trackbck_chkBox.Value == 0
                    handles = process_all_UltraTrack(hObject, eventdata, handles);
                else
                    handles = process_all_UltraTrack_backwards(hObject, eventdata, handles);
                end
                %try estimation (depending on hough tracked or not
                try
                    % State estimation
                    handles = do_state_estimation(hObject, eventdata, handles);
                catch
                    fprintf('No TimTrack, so no estimator ran\n')
                end
                % update the image and data using functions
                show_data(hObject, handles);
                show_image(hObject, handles);
            end
        end

    end
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on selection change in ROI_type.
function ROI_type_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROI_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROI_type

ROI_options =  get(hObject,'String');
i = get(hObject, 'Value');

handles.ROItype = ROI_options{i};

% Run TimTrack in case some already ran UltraTrack for state estimator
% if isfield(handles, 'Region') && ~isfield(handles,'geofeatures') 
%     try     
%         handles = process_all_Callback(hObject, eventdata, handles);   
% 
%     catch        
         fprintf('ROI type chagned \n')
%     end
% end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ROI_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

ROI_options =  get(hObject,'String');
i = get(hObject, 'Value');

handles.ROItype = ROI_options{i};

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SA.
function SA_Callback(hObject, eventdata, handles)
% hObject    handle to SA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Qlog = 10.^(-7:3);
Qlin1 = (1:9) * 10^-3;
Qlin2 = (1:9) * 10^-4;

Qs = [sort([Qlog Qlin1 Qlin2]) inf];

color = parula(length(Qs));

for i = 1:length(Qs)
    handles.Q = Qs(i);

    handles = do_state_estimation(hObject, eventdata, handles);

    %         figure(2)
    FL(:,i) = handles.Region(1).fas_length(handles.start_frame:end);
    PEN(:,i) = handles.Region(1).fas_pen(handles.start_frame:end);

    % save
    Save_As_Mat_Callback(hObject, eventdata, handles)

end

N = size(FL,1);
dt = 1/ handles.FrameRate;
t = 0:dt:((N-1)*dt);

Wn = 1.5*handles.fc_lpf / (.5 * handles.FrameRate);
Wn(Wn>=1) = 1-1e-6;
Wn(Wn<=0) = 1e-6;
[b,a] = butter(2, Wn, 'high');

PEN_low = filtfilt(b,a,PEN);
FL_low = filtfilt(b,a,FL);

% noise = [std(FL_low, 0,1); std(PEN_low, 0,1)];
noise = [mean(abs(diff(diff(FL)))); mean(abs(diff(diff(PEN))))];
drift = [abs(trapz(t, FL-FL(:,end))); abs(trapz(t, PEN-PEN(:,end)))] / max(t);

%%
if ishandle(2), close(2); end; figure(2)
subplot(521)
set(gca,'colororder',parula(length(Qs))); hold on
plot(t, FL,'linewidth',2);
xlabel('Time (s)')
ylabel('Length (mm)');
title('Fascicle length')

subplot(522)
set(gca,'colororder',parula(length(Qs))); hold on
plot(t, PEN,'linewidth',2);
xlabel('Time (s)')
ylabel('Angle (deg)');
title('Fascicle angle')

units = {'(mm)', '(deg)'};
c = 0:.1:1;
colors = cool(length(c));

for j = 1:2
    subplot(5,2,2+j);
    semilogx(Qs, drift(j,:),'linewidth',2)
    xlabel('Q value')
    ylabel(['Drift ', units(j)]);

    subplot(5,2,4+j);
    semilogx(Qs, noise(j,:),'linewidth',2);
    xlabel('Q value')
    ylabel(['Noise ', units(j)]);

    drift_n = drift ./ std(drift,1,2);
    noise_n = noise ./ std(noise,1,2);

    cost = c(:) * drift_n(j,1:end-1) + (1-c(:))*noise_n(j,1:end-1);

    [mcost, id] = min(cost, [], 2);

    subplot(5,2,6+j);
    set(gca,'colororder', cool(length(c)),'xscale','log'); hold on
    semilogx(Qs(1:end-1), cost,'linewidth',2); hold on

    for i = 1:length(c)
        semilogx(Qs(id(i)), mcost(i),'o','color',colors(i,:).^2,'markerfacecolor',colors(i,:).^2)
    end

    xlabel('Q value')
    ylabel('Drift + Noise');

    subplot(5,2,8+j)
    semilogy(c, Qs(id),'linewidth',2);

    xlabel('Drift weighting')
    ylabel('Optimal Q');

end


for i = 1:10
    subplot(5,2,i)
    box off
    axis tight
end


% --------------------------------------------------------------------
function set_Tim_Track_Callback(hObject, eventdata, handles)
% hObject    handle to set_Tim_Track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if the cross is pressed nothing is updated anyway
parms = adjust_hough_parameters(handles.parms);
%overwrite updated TimTrack parms in the main UTT folder 
filename = [mfilename,'.m'];
fullpath = which(filename);
mainfoldername = erase(fullpath,filename);
handles.parms = parms;  
save([mainfoldername 'TimTrack_parms.mat'],'parms');
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in manual_estimate.
function manual_estimate_Callback(hObject, eventdata, handles)
% hObject    handle to manual_estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'h')
    if ~isvalid(handles.h)
        handles = rmfield(handles, 'h');
    end
end
if isfield(handles,'d')
    if ~isvalid(handles.d)
        handles = rmfield(handles, 'd');
    end
end
if isfield(handles,'s')
    if ~isvalid(handles.s)
        handles = rmfield(handles, 's');
    end
end

try
    i = 1;
    j = 1;
    
    frame_no = handles.start_frame + round(get(handles.frame_slider,'Value')) - 1;
    
    if ~isfield(handles, 'Region')
        handles = Auto_Detect_Callback(hObject, eventdata, handles);
    end
    
    if ~isfield(handles, 'h') || ~isfield(handles, 'd') || ~isfield(handles, 's')
        Supex = handles.Region(i).sup_x{frame_no};
        Supey = handles.Region(i).sup_y{frame_no};
        Deepx = handles.Region(i).deep_x{frame_no};
        Deepy = handles.Region(i).deep_y{frame_no};
        Fasx = handles.Region(i).Fascicle(j).fas_x{frame_no};
        Fasy = handles.Region(i).Fascicle(j).fas_y{frame_no};
    
        d = round(size(handles.ImStack,2)/2);
        axes(handles.axes1)
        handles.h = drawline('Position', [Fasx(1)+d Fasy(1); Fasx(2)+d Fasy(2)], 'color', 'red', 'linewidth',2,'StripeColor','w');
        handles.s = drawline('Position', [Supex(1)+d Supey(1); Supex(2)+d Supey(2)], 'color', 'blue', 'linewidth',2,'StripeColor','w');
        handles.d = drawline('Position', [Deepx(1)+d Deepy(1); Deepx(2)+d Deepy(2)], 'color', 'green', 'linewidth',2,'StripeColor','w');
    else
        
        handles = extract_estimates(hObject, eventdata, handles);
    
        % if the first frame, accept manual tracking
        if frame_no == handles.start_frame
            handles.Region(i).Fascicle(j).fas_x{frame_no} = handles.Region(i).Fascicle(j).fas_x_manual{frame_no};
            handles.Region(i).Fascicle(j).fas_y{frame_no} = handles.Region(i).Fascicle(j).fas_y_manual{frame_no};
        end
        
        handles = calc_fascicle_length_and_pennation(handles,frame_no);

        try
            handles = do_state_estimation(hObject, eventdata, handles);
        catch
            disp('Estimation avaiable only with ROI Type "Hough - Local" or "Hough - global"')
        end
    
        delete(handles.h)
        delete(handles.d)
        delete(handles.s)
    
        handles = rmfield(handles, 'h');
        handles = rmfield(handles, 'd');
        handles = rmfield(handles, 's');
    end
    
    show_data(hObject, handles);
    show_image(hObject, handles);
    
    guidata(hObject, handles);

catch
    warndlg('No video loaded!')
end
function manual_variance_Callback(hObject, eventdata, handles)
% hObject    handle to manual_variance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of manual_variance as text
%        str2double(get(hObject,'String')) returns contents of manual_variance as a double

handles.R_manual = abs(str2double(get(hObject,'String')));

try
    handles = do_state_estimation(hObject, eventdata, handles);
catch
    disp('No tracking yet')
end

show_data(hObject, handles);
show_image(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function manual_variance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to manual_variance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.R_manual = abs(str2double(get(hObject,'String')));

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in clear_manual.
function clear_manual_Callback(hObject, eventdata, handles)
% hObject    handle to clear_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

i = 1;
j = 1;

N = handles.NumFrames + handles.start_frame - 1;
handles.Region(i).fas_length_manual    = nan(N,1);
handles.Region(i).fas_pen_manual       = nan(N,1);

rmfields = {'sup_x_manual', 'sup_y_manual','deep_x_manual','deep_y_manual','ROIx_manual','ROIy_manual','fas_ang_manual'};

for j = 1:length(rmfields)
    handles.Region(i).(rmfields{j}) = [];
end

rmfields = {'fas_x_manual','fas_y_manual'};
for j = 1:length(rmfields)
    handles.Region(i).Fascicle.(rmfields{j}) = [];
end

show_data(hObject, handles);
show_image(hObject, handles);

try
    handles = do_state_estimation(hObject, eventdata, handles);
catch
    disp('No tracking yet')
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Spatial_calibration_Callback(hObject, eventdata, handles)
% hObject    handle to Spatial_calibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%just grub the shown frame even tho it's not strictly necessary for
%calibration, any frame would be fine (Ideally people don't change setting
%during a recording).

if isfield(handles,'ImStack')
        % find current frame number from slider
        frame_no = round(get(handles.frame_slider,'Value'));
    
        %msgbox('Please click twice to define length on the image')
        % select two point for calculating the calibration mmperpx
        set(handles.axes1);
        [~,distanceInPixels]=ginputYellow(2);
    
        %create a simple dlg to type the real value
        %Ask the user for the real-world distance.
        userPrompt = {'Enter real distance in mm'};
        dialogTitle = 'Calibration';
        def = {''};
        answer = inputdlg(userPrompt, dialogTitle, 1, def);
    
        while isnan(str2double(answer{1}))  %check if it's non numeric
            answer = inputdlg(userPrompt, dialogTitle, 1, def);
        end
        %get the answer
        dist_mm = str2double(answer{1});
    
        %calculate mmperpx factor
        calibration_value = dist_mm /  abs(round(diff(distanceInPixels)));
        calibration_value = round(calibration_value,3); %round otherwise the conversion crash in the calculation (don't ask why)
        
        ImDepth = handles.vidHeight .* calibration_value; %calculate img_depth
        %update handles and GUI
        set(handles.ImDepthEdit,"String",string(ImDepth)); %this should automatically update the plots with Fascicle data
        ImDepthEdit_Callback(handles.ImDepthEdit, [], handles); % Manually call the callback function because it doesn't do automatically, this also updates the field and the plot
else
       warndlg('No video loaded')
end



%%%% Ultratrack (KLT optic flow)
function [handles] = process_all_UltraTrack_backwards(hObject, eventdata, handles)

%need to flip here for optic flow?
% if handles.trackbck_chkBox.Value == 1
%     handles.geofeatures = flip(handles.geofeatures);
%     handles.ImStack = flip(handles.ImStack,3);
% end

im1 = handles.ImStack(:,:,handles.NumFrames-1);
h = waitbar(1,['Processing frame '  num2str(handles.NumFrames) '/', num2str(handles.NumFrames)],'Name','Running UltraTrack...'); 

tstart = tic;

frames = (handles.NumFrames+handles.start_frame)-1 : -1 : handles.start_frame; %<-- the last frame is treated as 1because start frame is only an index

% numIterations = length(frames);

for i = 1:length(handles.Region)
    for f = frames

        % extract image
        im = handles.ImStack(:,:,f);

        % make a copy
        I_fmasked = im;

        if isfield(handles,'geofeatures') && strcmp(handles.ROItype, 'Hough - local')
            M = zeros(size(im1,1), size(im1,2), handles.parms.fas.npeaks);

            for j = 1:handles.parms.fas.npeaks
                x1 = handles.geofeatures(f).x(j,1) * handles.imresize_fac;
                y1 = handles.geofeatures(f).y(j,1) * handles.imresize_fac;

                x2 = handles.geofeatures(f).x(j,2) * handles.imresize_fac;
                y2 = handles.geofeatures(f).y(j,2) * handles.imresize_fac;

                dy = 5; %what is 5?

                ROIx = [x1 x1 x2 x2 x1];
                ROIy = [y1-dy y1+dy y2+dy y2-dy y1-dy]';

                if sum(isfinite(ROIx)) == length(ROIx) && sum(isfinite(ROIy)) == length(ROIy)
                    M(:,:,j) = poly2mask(ROIx,ROIy, size(im1,1), size(im1,2));
                end
            end

            % mask
            fmask = sum(M,3);
            fmask(fmask>1) = 1;

            I_fmasked(fmask~=1) = 0;
        end


        % get current ROI
        if contains(handles.ROItype, 'Hough')
            I_amasked = im;

            n = handles.vidWidth;
            ROIx = [1 1 n n 1]';
            super_apo = handles.geofeatures(f).super_pos';
            deep_apo = handles.geofeatures(f).deep_pos';
            thickness = deep_apo - super_apo;

            r = .1;
            ROIy_fcor = round([super_apo(1)+thickness(1)*r; deep_apo-thickness*r; super_apo([2,1])+thickness([2,1])*r]);

            % mask
            fmask = poly2mask(ROIx, ROIy_fcor, size(im,1), size(im,2));
            I_fmasked(fmask~=1) = 0;

            m = handles.vidHeight;

            s = handles.parms.apo.super.cut;
            ROIys = [s(1) s(2) s(2) s(1) s(1)] * m;

            d = handles.parms.apo.deep.cut;
            ROIyd = [d(1) d(2) d(2) d(1) d(1)] * m;

            dmask =  poly2mask(ROIx, ROIyd, size(im,1), size(im,2));
            smask =  poly2mask(ROIx, ROIys, size(im,1), size(im,2));

            amask = dmask + smask;
            amask(amask > 1) = 1;
            I_amasked(amask~=1) = 0;
        end


        for j = 1:length(handles.Region(i).Fascicle)

            if f == (handles.NumFrames+handles.start_frame)-1
                % detect points
                fpoints = detectMinEigenFeatures(I_fmasked,'FilterSize',11, 'MinQuality', 0.005);

                if contains(handles.ROItype, 'Hough')
                    fpoints = fpoints.selectStrongest(300);
                end

                % get location
                fpoints = double(fpoints.Location);

                % must be in ROI
                ROIy = handles.Region(i).ROIy{f};
                ROIx = handles.Region(i).ROIx{f};
                inPoints = inpolygon(fpoints(:,1), fpoints(:,2), ROIx, ROIy);
                fpoints = fpoints(inPoints,:);

                % define fascicle tracker
                fpointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',handles.BlockSize);
                initialize(fpointTracker,fpoints,im);

                if contains(handles.ROItype, 'Hough')
                    % define aponeurosis tracker
                    apoints = detectMinEigenFeatures(I_amasked,'FilterSize',11, 'MinQuality', 0.005);
                    apoints =  double(apoints.Location);
                    apointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',handles.BlockSize);
                    initialize(apointTracker,apoints,im);
                end
            else

                if ~exist('fpointTracker','var')

                    fpoints = detectMinEigenFeatures(I_fmasked,'FilterSize',11, 'MinQuality', 0.005);

                    if contains(handles.ROItype, 'Hough')
                        fpoints = fpoints.selectStrongest(300);
                    end

                    % get location
                    fpoints = double(fpoints.Location);

                    % must be in ROI
                    ROIy = handles.Region(i).ROIy{f};
                    ROIx = handles.Region(i).ROIx{f};
                    inPoints = inpolygon(fpoints(:,1), fpoints(:,2), ROIx, ROIy);
                    fpoints = fpoints(inPoints,:);

                    % define fascicle tracker
                    fpointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',handles.BlockSize);
                    initialize(fpointTracker,fpoints,im);

                    if contains(handles.ROItype, 'Hough')
                        % define aponeurosis tracker
                        apoints = detectMinEigenFeatures(I_amasked,'FilterSize',11, 'MinQuality', 0.005);
                        apoints =  double(apoints.Location);
                        apointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',handles.BlockSize);
                        initialize(apointTracker,apoints,im);
                    end

                end
                % Compute the flow and new roi
                [fpointsNew, isFound] = step(fpointTracker, im);
                [wf,~] = estimateGeometricTransform2D(fpoints(isFound,:), fpointsNew(isFound,:), 'affine', 'MaxDistance',50);
                handles.Region(i).warp(:,:,f+1) = wf;

                if contains(handles.ROItype, 'Hough')
                    % Compute the flow and new roi
                    [apointsNew, isFound] = step(apointTracker, im);
                    [wa,~] = estimateGeometricTransform2D(apoints(isFound,:), apointsNew(isFound,:), 'affine', 'MaxDistance',50);
                    handles.Region(i).awarp(:,:,f+1) = wa;
                end

                % apply the warp to fascicles
                fas_prev = [handles.Region(i).Fascicle(j).fas_x{f+1} handles.Region(i).Fascicle(j).fas_y{f+1}];
                fas_new = transformPointsForward(wf, fas_prev);

                % save
                handles.Region(i).Fascicle(j).fas_x{f} = fas_new(:,1);
                handles.Region(i).Fascicle(j).fas_y{f} = fas_new(:,2);

                % make a copy
                handles.Region(i).Fascicle(j).fas_x_original{f} = handles.Region(i).Fascicle(j).fas_x{f};
                handles.Region(i).Fascicle(j).fas_y_original{f} = handles.Region(i).Fascicle(j).fas_y{f};

                % in the new version, aponeurosis has its own warp (in the old version it doesnt)
                if contains(handles.ROItype, 'Hough')

                    % apply warp to aponeurosis
                    super_prev = [handles.Region(i).sup_x{f+1} handles.Region(i).sup_y{f+1}];
                    super_new = transformPointsForward(wa, super_prev);

                    deep_prev = [handles.Region(i).deep_x{f+1} handles.Region(i).deep_y{f+1}];
                    deep_new = transformPointsForward(wa, deep_prev);

                    % save
                    handles.Region(i).sup_x{f} = [1 n]';
                    handles.Region(i).sup_y{f} = super_new(:,2);
                    handles.Region(i).deep_x{f} = [1 n]';
                    handles.Region(i).deep_y{f} = deep_new(:,2);

                    % get ROI
                    ROIx = [1 1 n n 1]';
                    super_apo = handles.geofeatures(f).super_pos';
                    deep_apo = handles.geofeatures(f).deep_pos';

                    ROIy = [super_apo(1); deep_apo; super_apo([2,1])];
                    ROIy_fcor = round([super_apo(1)+thickness(1)*r; deep_apo-thickness*r; super_apo([2,1])+thickness([2,1])*r]);

                else
                    % apply warp to ROI
                    ROIpos = transformPointsForward(wf, [handles.Region(i).ROIx{f+1} handles.Region(i).ROIy{f+1}]);

                    ROIx = ROIpos(:,1);
                    ROIy = ROIpos(:,2);

                    ROIx(ROIx > handles.vidWidth) = handles.vidWidth;
                    ROIy(ROIy > handles.vidHeight) = handles.vidHeight;
                    ROIx(ROIx < 1) = 1;
                    ROIy(ROIy < 1) = 1;

                    handles.Region(i).sup_x{f}  = ROIx([1,4]);
                    handles.Region(i).sup_y{f}  = ROIy([1,4]);
                    handles.Region(i).deep_x{f} = ROIx([2,3]);
                    handles.Region(i).deep_y{f} = ROIy([2,3]);
                end

                % save ROI
                handles.Region(i).ROIx{f} = ROIx;
                handles.Region(i).ROIy{f} = ROIy;

                % calculate the length and pennation for the current frame
                handles = calc_fascicle_length_and_pennation(handles,f);

                % update the points
                fpoints = fpointsNew;

                % if drops below 100, define new points
                if length(fpoints) < 100 || ~strcmp(handles.ROItype(1:5), 'Hough')

                    % detect points
                    fpoints = detectMinEigenFeatures(I_fmasked,'FilterSize',11, 'MinQuality', 0.005);

                    if strcmp(handles.ROItype(1:5), 'Hough')
                        fpoints = fpoints.selectStrongest(300);
                    end

                    fpoints = double(fpoints.Location);
                end

                % must be in ROI
                if strcmp(handles.ROItype(1:5), 'Hough')
                    inPoints = inpolygon(fpoints(:,1),fpoints(:,2), ROIx, ROIy_fcor);
                else
                    inPoints = inpolygon(fpoints(:,1),fpoints(:,2), ROIx, ROIy);
                end

                fpoints = fpoints(inPoints,:);

                % set tracker
                setPoints(fpointTracker, fpoints);

                % update the points
                if strcmp(handles.ROItype(1:5), 'Hough')
                    apoints = apointsNew;

                    if length(apoints) < 500
                        % detect points
                        apoints = detectMinEigenFeatures(I_amasked,'FilterSize',11, 'MinQuality', 0.005);
                        apoints = double(apoints.Location);
                    end

                    m = handles.vidHeight;

                    s = handles.parms.apo.super.cut;
                    ROIys = [s(1) s(2) s(2) s(1) s(1)]' * m;

                    d = handles.parms.apo.deep.cut;
                    ROIyd = [d(1) d(2) d(2) d(1) d(1)]' * m;

                    % must be in ROI
                    dinPoints = inpolygon(apoints(:,1),apoints(:,2), ROIx, ROIyd);
                    sinPoints = inpolygon(apoints(:,1),apoints(:,2), ROIx, ROIys);
                    apoints = apoints(dinPoints | sinPoints,:);

                    % set tracker
                    setPoints(apointTracker, apoints);
                end
            end

            % save the points
            handles.points{f} = fpoints;

            if strcmp(handles.ROItype(1:5), 'Hough')
                handles.apoints{f} = apoints;
            end

        end

        frac_progress = ((f-handles.start_frame)+(get(handles.frame_slider,'Max')*(i-1))) / (get(handles.frame_slider,'Max')*length(handles.Region));
        waitbar(frac_progress,h, ['Processing frame ', num2str((f-handles.start_frame)), '/', num2str(get(handles.frame_slider,'Max'))])
    end

end

close(h)
handles.ProcessingTime(2) = toc(tstart);
% --- Executes during object creation, after setting all properties.
function trackbck_chkBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackbck_chkBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in trackbck_chkBox.
function trackbck_chkBox_Callback(hObject, eventdata, handles)
% hObject    handle to trackbck_chkBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = ~handles.trackbck_chkBox.Value; %trick to force previous value (in case tracking forwards, and then some check and change variance) 
if isfield(handles,'Region')

        %check whether all frames have been already tracked with
        %opticflow, if yes re-run it with the new block size
        if sum(~cellfun(@isempty, handles.Region.Fascicle.fas_x, 'UniformOutput', true)) >= handles.NumFrames
            %if size(handles.Region(i).Fascicle.analysed_frames,2) > 0 %double check this
            %force to not be possible to track backwards to avoid messy
            %estimations
            set(handles.trackbck_chkBox,'Value',tmp)
        
        end

end
% Hint: get(hObject,'Value') returns toggle state of trackbck_chkBox
 
