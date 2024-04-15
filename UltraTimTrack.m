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

% check to make sure there is a settings mat-file present, if not then make
% one in the directory where the m-file sits.
[program_directory, ~, ~] = fileparts(mfilename('fullpath'));

if ismac || isunix
    if isempty(dir([program_directory '/ultrasound_tracking_settings.mat']))
        ImageDepth = 65;
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
        ImageDepth = 65;
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
    if isfield(TVD.TVDdata,'cmPerPixY') %check whether the field exists and update scalar
        %ImageDepth = round(TVD.TVDdata.cmPerPixY*10,3); %round to 3 digits
        ImageDepth = round(TVD.TVDdata.Height * TVD.TVDdata.cmPerPixY,3)*10;
        handles.ID = ImageDepth;
        set(handles.ImDepthEdit,'String',num2str(ImageDepth));

        % Update handles structure
        guidata(hObject, handles);
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

handles.Xmin = 0;
% handles.Xmax = handles.vidWidth;

cd(handles.pname)

% make a timeline which corresponds to the ultrasound frame data
if strcmp(Ext,'.mat')
    handles.Time = handles.TimeStamps;
elseif exist('TVD','var') && isfield(TVD.TVDdata,'Time') %check if also time exists as Telemed has uncostant framerate
    handles.Time = TVD.TVDdata.Time(1:end-1); %-1 because last timestamp is repeated or UT doesn't load all of them
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
    if isfield(handles,'Region')
        handles = rmfield(handles,'Region');
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

handles.ImStack = autocrop_Tim(handles.ImStack);
handles.vidHeight = size(handles.ImStack,1);
handles.vidWidth = size(handles.ImStack,2);

% Im = handles.ImStack(:,:,1);
% 
% handles.vidHeight = size(Im,1);
% 
% %corners = detectAutoCrop(handles.movObj);
% handles.crop_rect = autocrop(handles.ImStack,handles.movObj.FrameRate/2);
% 
%     %Crop all images before updating
%     for ii = 1 : handles.NumFrames
%         tmp(:,:,ii) = imcrop(handles.ImStack(:,:,ii),handles.crop_rect);
%     end
% 
%     handles.ImStack = tmp; %overwrite the one used thourghout the script
%     handles.vidHeight = handles.crop_rect(4);
%     handles.vidWidth = handles.crop_rect(3);
%     clearvars tmp
%     % Clean axis from original image and tight axis on the cropped image
%     cla
%     % update the image axes using show_image function (bottom)
%     show_image(hObject,handles);
%     axis tight

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
function menu_save_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'Region')

    for i = 1:length(handles.Region)

        if isfield(handles,'ImStack') && isfield(handles.Region(i),'fas_length')

            time = handles.Time;

            %determine any non-zero entries in fascicle length array
            nz = logical(handles.Region(i).fas_length(:,1) ~= 0);

            T = time(nz)';

            if isfield(handles.Region(i),'fas_length_corr')
                R(i).FL = handles.Region(i).fas_length_corr(nz,:)';
                R(i).PEN = handles.Region(i).fas_pen_corr(nz,:)';
            else R(i).FL = handles.Region(i).fas_length(nz,:)';
                R(i).PEN = handles.Region(i).fas_pen(nz,:)';
            end
            
            Save_As_Txt_Callback(hObject,eventdata,handles);

            

        end

    end
end

% if strcmp(handles.ext,'.c3d')
%     % create a new c3d file or overwrite the existing one
%     fileout_suggest = [handles.pname handles.dat_file];
%     [fileout, pathout] = uiputfile(fileout_suggest,'Save fascicle data as...');
% 
%     % write the data to the c3dfile specified
%     cd(pathout);
%     btkWriteAcquisition(c3d_handle,fileout);
% 
% end

% --------------------------------------------------------------------
function menu_load_fascicle_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load_fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'));
    [fname, pname] = uigetfile('*.mat','Load tracking MAT file');
    load([pname fname]);

    if exist('Fdat','var') && isfield(Fdat,'Region') && isfield(Fdat.Region,'Fascicle')

        for i = 1:length(Fdat.Region)

            for k = 1:length(Fdat.Region(i).Fascicle)

                handles.Region(i).Fascicle(k).fas_x{frame_no} = Fdat.Region(i).Fascicle(k).fas_x;
                handles.Region(i).Fascicle(k).fas_y{frame_no} = Fdat.Region(i).Fascicle(k).fas_y;

                handles.Region(i).Fascicle(k).current_xy(1,1) = handles.Region(i).Fascicle(k).fas_x{frame_no}(1);
                handles.Region(i).Fascicle(k).current_xy(1,2) = handles.Region(i).Fascicle(k).fas_y{frame_no}(1);
                handles.Region(i).Fascicle(k).current_xy(2,1) = handles.Region(i).Fascicle(k).fas_x{frame_no}(2);
                handles.Region(i).Fascicle(k).current_xy(2,2) = handles.Region(i).Fascicle(k).fas_y{frame_no}(2);

                handles.Region(i).fas_pen(frame_no,k) = atan2(abs(diff(handles.Region(i).Fascicle(k).fas_y{frame_no})),abs(diff(handles.Region(i).Fascicle(k).fas_x{frame_no})));
                scalar = handles.ID;%/handles.vidHeight;
                handles.Region(i).fas_length(frame_no,k) = scalar*sqrt(diff(handles.Region(i).Fascicle(k).fas_y{frame_no}).^2 + diff(handles.Region(i).Fascicle(k).fas_x{frame_no}).^2);




                if ~isfield(handles.Region(i).Fascicle(k),'analysed_frames')
                    handles.Region(i).Fascicle(k).analysed_frames = frame_no;
                else
                    handles.Region(i).Fascicle(k).analysed_frames = sort([handles.Region(i).Fascicle(k).analysed_frames frame_no]);
                end
            end

            Nfascicle(i) = length(handles.Region(i).Fascicle);

            handles.Region(i).ROIx{frame_no} = Fdat.Region(i).ROIx;
            handles.Region(i).ROIy{frame_no} = Fdat.Region(i).ROIy;
            %handles.Region(i).ROI{frame_no} = Fdat.Region(i).ROI;


        end

          %% Add multiple regions and fascicles in a later version
%         Nregions = length(handles.Region);
%         Nfas = max(Nfascicle);
% 
%         
%         % Correct the drop down lists for regions and fascicles to match
%         %imported data
%         ROItoCorrString{1} = 'all';
%         for i = 1:Nregions
%             ROIlistString{i} = num2str(i);
%             ROItoCorrString{i+1} = num2str(i);
%         end
%         set(handles.no_tracked_regions,'String',ROIlistString);
%         set(handles.RegionsToCorrect,'String',ROItoCorrString);
% 
%         FAStoCorrString{1} = 'all';
%         for i = 1:Nfas
%             FASlistString{i} = num2str(i);
%             FAStoCorrString{i+1} = num2str(i);
%         end
%         set(handles.no_tracked_fascicles,'String',FASlistString);
%         set(handles.FasciclesToCorrect,'String',FAStoCorrString);
    else

        warndlg('Fascicle data not found','ERROR')

    end


    show_image(hObject,handles);

end

% --------------------------------------------------------------------
function menu_save_fascicle_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'movObj')||isfield(handles,'BIm')||isfield(handles,'ImStack')

    % find current frame number from slilder
    frame_no = round(get(handles.frame_slider,'Value'));

    for i = 1:length(handles.Region)
        for k = 1:length(handles.Region(i).Fascicle)

            if isfield(handles.Region(i).Fascicle(k),'fas_x')

                Fdat.Region(i).Fascicle(k).fas_x = handles.Region(i).Fascicle(k).fas_x{frame_no};
                Fdat.Region(i).Fascicle(k).fas_y = handles.Region(i).Fascicle(k).fas_y{frame_no};

            end

            %% TimTrack_downsample
            %Fdat.Region(i).Fascicle(k).fas_x_Hough = handles.Region(i).Fascicle(k).fas_x_Hough{frame_no};
            %Fdat.Region(i).Fascicle(k).fas_y_Hough = handles.Region(i).Fascicle(k).fas_y_Hough{frame_no};


        end
        % also save the ROI data
        Fdat.Region(i).ROIx = handles.Region(i).ROIx{frame_no};
        Fdat.Region(i).ROIy = handles.Region(i).ROIy{frame_no};
        %Fdat.Region(i).ROI = handles.Region(i).ROIp{frame_no};
    end

    % save to file
    [fname, pname] = uiputfile('*.mat', 'Save fascicle to MAT file');
    save([pname fname],'Fdat');

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
%     spos = [handles.S.Position(1:2) ceil(size(handles.ImStack,2)/2) handles.S.Position(4)];
%     dpos = [handles.D.Position(1:2) ceil(size(handles.ImStack,2)/2) handles.D.Position(4)];
% 
%     handles.ImTrack(spos(2):(spos(2)+spos(4)),spos(1):spos(3),3,:) = 230;
%     handles.ImTrack(dpos(2):(dpos(2)+dpos(4)),dpos(1):dpos(3),2,:) = 230;

%     [fileout, pathout, FI] = uiputfile('*.mp4', 'Save video as');
%     filename = [pathout fileout];
%     
    filename = [handles.pname, handles.fname(1:end-4), '_tracked_Q=',strrep(num2str(handles.Q),'.','')];
    vidObj = VideoWriter(filename,'MPEG-4');
    vidObj.FrameRate = handles.FrameRate;
    open(vidObj);

%     if FI > 0

        h = waitbar(0,['Saving frame 1/', num2str(handles.NumFrames)],'Name','Saving to video file...');
        i = 1;
        j = 1;
        
        for f = handles.start_frame:1:get(handles.frame_slider,'Max')
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

            % add ROI
            currentImage = insertShape(currentImage,'Polygon',[handles.Region(i).ROIx{f}(1)+d, handles.Region(i).ROIy{f}(1), ...
            handles.Region(i).ROIx{f}(2)+d, handles.Region(i).ROIy{f}(2),handles.Region(i).ROIx{f}(3)+d, handles.Region(i).ROIy{f}(3),...
            handles.Region(i).ROIx{f}(4)+d, handles.Region(i).ROIy{f}(4),handles.Region(i).ROIx{f}(5)+d, handles.Region(i).ROIy{f}(5)],'LineWidth',1, 'Color','red');
        
            % save
            ImTrack = currentImage;
            
            % show region
            spos = floor([handles.S.Position(1:2) ceil(size(handles.ImStack,2)/2) handles.S.Position(4)]);
            dpos = floor([handles.D.Position(1:2) ceil(size(handles.ImStack,2)/2) handles.D.Position(4)]);

            ImTrack(spos(2):(spos(2)+spos(4)),spos(1):(spos(1)+spos(3)),3) = 230;
            ImTrack(dpos(2):(dpos(2)+dpos(4)),dpos(1):(dpos(1)+dpos(3)),2) = 230;

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
function Save_Tracking_As_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Tracking_As (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
            nz = logical(handles.Region(i).fas_length(:,1) ~= 0);

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
            Fdat.Region(i).Time = time(nz);
        end
    end
end

filename = [handles.pname, handles.fname(1:end-4), '_tracked_Q=',strrep(num2str(handles.Q),'.','')];
save(filename,'TrackingData','Fdat');
% [~,handles.file,handles.ext] = fileparts(handles.fname);
% fileout_suggest = [handles.pname handles.file '.mat'];
% [fileout, pathout] = uiputfile(fileout_suggest,'Save fascicle data as...');
% cd(pathout);
% 
% if exist([pathout fileout], 'file') == 2
% 
%     if exist('Data','var') && exist('Fdat','var')
%         save(fileout,'-append','Fdat','Data','TrackingData')
%     elseif exist('Fdat','var')
%         save(fileout,'-append','Fdat','TrackingData')
%     elseif exist('Data','var')
%         save(fileout,'-append','Data')
%     else
%         warndlg('There was no data to save', 'WARNING!')
%     end
% 
% else
% 
%     if exist('Data','var') && exist('Fdat','var')
%         save(fileout,'Fdat','Data','TrackingData')
%     elseif exist('Fdat','var')
%         save(fileout,'Fdat','TrackingData')
%     elseif exist('Data','var')
%         save(fileout,'Data')
%     else
%         warndlg('There was no data to save', 'WARNING!')
%     end

% end


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
show_data(hObject, handles);

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

menu_clear_tracking_Callback(hObject, eventdata, handles)% clear any current tracking

%select file
[fname, pname] = uigetfile('*.mat', 'Pick a .MAT file');
load([pname fname],'TrackingData','Fdat');

handles.Region = Fdat.Region;
handles.start_frame = TrackingData.start_frame;
handles.NumFrames = TrackingData.NumFrames;
handles.S = TrackingData.S;
handles.D = TrackingData.D;
set(handles.frame_slider,'Min',1);
set(handles.frame_slider,'Max',handles.NumFrames);
set(handles.frame_slider,'Value',1);
set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 5/handles.NumFrames]);
% set the string in the frame_number to 1
set(handles.frame_number,'String',1);

show_image(hObject,handles);

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

axes(handles.length_plot); hold off;

if isfield(handles, 'Region')
for i = 1:length(handles.Region)
    if isfield(handles.Region(i),'Fascicle')
        if isfield(handles.Region(i),'fas_length')

            if ~isempty(handles.Region(i).fas_length)

                dt = 1/handles.FrameRate;
                nz = logical(handles.Region(i).fas_length(:,1) ~= 0);
                FL = handles.Region(i).fas_length(nz,:);
                PEN = handles.Region(i).fas_pen(nz,:) * 180/pi;

                time = 0:dt:((handles.NumFrames-1)*dt);

                for j = 1:length(handles.Region(i).Fascicle)

                    plot(handles.length_plot,time(nz),FL,'r','linewidth',2);
                    set(handles.length_plot,'ylim',[min(FL)*0.85 max(FL)*1.15],'xlim', [0 max(time)],'box','off'); %set axis 15% difference of min and and value,easier to read
                    xlabel('Time (s)'); ylabel('Fascicle Length (mm)');
                    
                    plot(handles.mat_plot,time(nz),PEN,'r','linewidth',2);
                    set(handles.mat_plot,'ylim',[min(PEN)*0.85 max(PEN)*1.15], 'xlim', [0 max(time)],'box','off'); %set axis 15% difference of min and and value,easier to read
                    handles.mat_plot.YLabel.String = 'Fascicle Angle (deg)';
                    handles.mat_plot.XLabel.String = 'Time (s)';
                    
                end
            end
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
    frame_no = round(get(handles.frame_slider,'Value'));
    
    i = 1;
    j = 1;

    if isfield(handles, 'Region')
        if length(handles.Region(i).Fascicle(j).fas_x) >= frame_no
            % create analyzed frame
            d = round(size(handles.ImStack,2)/2);

            ZeroPadL = 200*ones(size(handles.ImStack,1), ceil(size(handles.ImStack,2)/2),'uint8');
            ZeroPadR = 200*ones(size(handles.ImStack,1), floor(size(handles.ImStack,2)/2),'uint8');

            f = frame_no;
            % add padding
            currentImage = [ZeroPadL, handles.ImStack(:,:,frame_no+handles.start_frame-1), ZeroPadR];

            % add fascicle
            currentImage = insertShape(currentImage,'line',[handles.Region(i).Fascicle(j).fas_x{f}(1)+d, handles.Region(i).Fascicle(j).fas_y{f}(1), ...
            handles.Region(i).Fascicle(j).fas_x{f}(2)+d,handles.Region(i).Fascicle(j).fas_y{f}(2)], 'LineWidth',5, 'Color','red');

            currentImage = insertMarker(currentImage,[handles.Region(i).Fascicle(j).fas_x{f}(1)+d, handles.Region(i).Fascicle(j).fas_y{f}(1);...
                handles.Region(i).Fascicle(j).fas_x{f}(2)+d, handles.Region(i).Fascicle(j).fas_y{f}(2)], 'o', 'Color','red','size',5);

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
            Im = handles.ImStack(:,:,frame_no+handles.start_frame-1);
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
        handles.S = images.roi.Rectangle(gca,'position', [1 handles.parms.apo.super.cut(1)*handles.vidHeight handles.vidWidth diff(handles.parms.apo.super.cut)*handles.vidHeight],'color','blue');
        handles.D = images.roi.Rectangle(gca,'position', [1 handles.parms.apo.deep.cut(1)*handles.vidHeight handles.vidWidth diff(handles.parms.apo.deep.cut)*handles.vidHeight],'color','green');
    end
    
    % Update Im and NIm
    handles.Im = Im;
    handles.NIm = handles.Im;
    
    % remove previous vertical lines
    children = get(handles.length_plot, 'children');
    if length(children) > 1
        delete(children(1));
    end

    children = get(handles.mat_plot, 'children');
    if length(children) > 1
        delete(children(1));
    end

    if isfield(handles,'Region')
        FL = handles.Region(1).fas_length;
        PEN = handles.Region(1).fas_pen * 180/pi;

        % add new vertical lines
        line(handles.length_plot, 'xdata', handles.Time(frame_no) * ones(1,2), 'ydata', [.85*min(FL) 1.15*max(FL)],'color',[0 0 0]);
        line(handles.mat_plot, 'xdata', handles.Time(frame_no) * ones(1,2), 'ydata', [.85*min(PEN) 1.15*max(PEN)],'color', [0 0 0]);
      
    end

    handles.prev_frame_no = frame_no;

    % Update handles structure
    guidata(hObject, handles);

end

% --- Executes on button press in process_all.
function[handles] = process_all_Callback(hObject, eventdata, handles)
% hObject    handle to process_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% detect first frame
handles = Auto_Detect_Callback(hObject, eventdata, handles);

% Run TimTrack
if ~isnan(handles.Q)
    handles = process_all_TimTrack(hObject, eventdata, handles);
    
    % low-pass filter ROI
    handles = lowpass_ROI(hObject, eventdata, handles);
end

% Run UltraTrack (note: includes state estimation on ROI)
handles = process_all_UltraTrack(hObject, eventdata, handles);

% State estimation
% handles = do_state_estimation(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% update the image axes using show_image function (bottom)
show_data(hObject, handles);

function[handles] = process_all_UltraTrack(hObject, eventdata, handles)

    %% Optical flow and state estimation
    % setup current and new image
    frame_no = 1;
    
    im1 = handles.ImStack(:,:,handles.start_frame);
    h = waitbar(0,['Processing frame 1/', num2str(handles.NumFrames)],'Name','Running UltraTrack...');
   

    for i = 1:length(handles.Region)
            
        points = detectMinEigenFeatures(im1,'FilterSize',11, 'MinQuality', 0.005);
        points = double(points.Location);
        [inPoints] = inpolygon(points(:,1),points(:,2), handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no});
        points = points(inPoints,:);

        %% Initialize point tracker and run opticflow
        tstart = tic;
        %calculate block size according to ROI 
        width       = floor(max(abs(diff(handles.Region.ROIx{frame_no}))) * 0.0725); %0.20
        height      = floor(max(abs(diff(handles.Region.ROIy{frame_no}))) * 0.08); %0.40; thickness changes?

        handles.BlockSize = [width height]; %save as width and height for later comparison

        % Ensure vidWidth and vidHeight are both odd numbers
        if mod(width, 2) == 0
            width = width + 1; % Increment by 1 to make it odd
        end

        if mod(height, 2) == 0
            height = height + 1; % Increment by 1 to make it odd
        end

        pointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',[height width]);
%         pointTracker = vision.PointTracker('NumPyramidLevels',1,'MaxIterations',10,'MaxBidirectionalError',inf,'BlockSize',[21 71]);
        initialize(pointTracker,points,im1);

        for f = frame_no+1:handles.NumFrames
            % Get the current image
            im2 = handles.ImStack(:,:,f+handles.start_frame-1);
            handles.NIm = im2;

            % Compute the flow and new roi
            [pointsNew, isFound] = step(pointTracker, im2);
            [w,~] = estimateGeometricTransform2D(points(isFound,:), pointsNew(isFound,:), 'affine', 'MaxDistance',50);
            handles.Region(i).warp(:,:,f-1) = w;

            for j = 1:length(handles.Region(i).Fascicle)

                % apply the warp to fascicles
                fas_prev = [handles.Region(i).Fascicle(j).fas_x{f-1}' handles.Region(i).Fascicle(j).fas_y{f-1}'];
                fas_new = transformPointsForward(w, fas_prev);
                handles.Region(i).Fascicle(j).fas_x{f} = fas_new(:,1)';
                handles.Region(i).Fascicle(j).fas_y{f} = fas_new(:,2)';
                
                % make a copy
                handles.Region(i).Fascicle(j).fas_x_original{f} =  handles.Region(i).Fascicle(j).fas_x{f};
                handles.Region(i).Fascicle(j).fas_y_original{f} =  handles.Region(i).Fascicle(j).fas_y{f};
                
                % calculate the length and pennation for the current frame
                handles.Region(i).fas_pen(f,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y{f})),...
                    abs(diff(handles.Region(i).Fascicle(j).fas_x{f})));

                scalar = handles.ID/handles.vidHeight;

                handles.Region(i).fas_length(f,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{f}).^2 +...
                    diff(handles.Region(i).Fascicle(j).fas_x{f}).^2);

                % optionally: calc ROI from optical flow
                if isnan(handles.Q)
                    ROI_prev = [handles.Region(i).ROIx{f-1} handles.Region(i).ROIy{f-1}];
                    ROI_new = transformPointsForward(w, ROI_prev);
                    handles.Region(i).ROIx{f} = ROI_new(:,1);
                    handles.Region(i).ROIy{f} = ROI_new(:,2);
                end
                    
                % set the points
                inPoints = inpolygon(points(:,1),points(:,2), handles.Region(i).ROIx{f}, handles.Region(i).ROIy{f});
                points = points(inPoints,:);
                setPoints(pointTracker, points);
                
            end

%             profile viewer
            frac_progress = (f+(get(handles.frame_slider,'Max')*(i-1))) / (get(handles.frame_slider,'Max')*length(handles.Region));
            waitbar(frac_progress,h, ['Processing frame ', num2str(f), '/', num2str(get(handles.frame_slider,'Max'))])
        end
        
    end
    
close(h)
handles.ProcessingTime(2) = toc(tstart);

function[handles] = process_all_TimTrack(hObject, eventdata, handles)

% run TimTrack on all frames
frames = 1:handles.NumFrames;
numIterations = length(frames);

parms = handles.parms;
parms.extrapolation = 0;

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

            tstart = tic;
            hwb = waitbar(0,'','Name','Running TimTrack...');
            for f = frames

                geofeatures(f) = auto_ultrasound(im2(:,:,f), parms);
                waitbar(f / numIterations, hwb, sprintf('Processing frame %d/%d', f, numIterations));

            end
            close(hwb)
            handles.ProcessingTime(1) = toc(tstart);
        end

       % Adjust the parameter of geofeatures
        for kk = 1:length(geofeatures)
            geofeatures(kk).super_coef(2) = geofeatures(kk).super_coef(2) * handles.imresize_fac;
            geofeatures(kk).deep_coef(2) = geofeatures(kk).deep_coef(2) * handles.imresize_fac;
            geofeatures(kk).thickness = geofeatures(kk).thickness * handles.imresize_fac;
        end

        handles.geofeatures = geofeatures;


    end
end

function[handles] = lowpass_ROI(hObject, eventdata, handles)
    i = 1;
    n = handles.vidWidth;

    frames = 1:handles.NumFrames;
    
    APO_TT = nan(5,length(frames));
    for f = frames
        handles.Region(i).ROIy{f} = round([polyval(handles.geofeatures(f).super_coef, 1) polyval(handles.geofeatures(f).deep_coef, [1 n]) polyval(handles.geofeatures(f).super_coef, [n 1])])';
        handles.Region(i).ROIx{f} = [1 1 n n 1]';
            
        APO_TT(:,f) = handles.Region(i).ROIy{f};
    end
        
    y = nan(size(APO_TT));

    if ~isnan(handles.fc_lpf)
        Wn = handles.fc_lpf / (.5 * handles.FrameRate);
        Wn(Wn>=1) = 1-1e-6;
        Wn(Wn<=0) = 1e-6;
        [b,a] = butter(2, Wn);
    
        for kk = 1:size(APO_TT,1)
            y(kk,:) = filtfilt(b,a,APO_TT(kk,:));
        end
    else
        y = APO_TT;
    end

    for f = frames
        handles.Region.ROIy{f} = round(y(:,f));
    end

    
function[handles] = estimate_variance(hObject, eventdata, handles)

geofeatures = handles.geofeatures;


x = nan(handles.NumFrames,1);
for f = 1:handles.NumFrames
    x(f,1) = geofeatures(f).alpha;
end

Wn = 1.5*handles.fc_lpf / (.5 * handles.FrameRate);
Wn(Wn>=1) = 1-1e-6;
Wn(Wn<=0) = 1e-6;
[b,a] = butter(2, Wn);

y = nan(size(x));
for i = 1:size(x,2)
    y(:,i) = filtfilt(b,a,x(:,i));
end

handles.R = var(x-y);

function[handles] = do_state_estimation(hObject, eventdata, handles)
% hObject    handle to do_state_estimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get Qmax
if ~isnan(handles.Q)
    handles = estimate_variance(hObject, eventdata, handles);

% forward state estimation
for f = 2:get(handles.frame_slider,'Max')
    for i = 1:length(handles.Region)
        for j = 1:length(handles.Region(i).Fascicle)
            % state estimation
            handles = state_estimator(handles,f,f-1,'forward');
        end
    end
end

% backward state estimation
for f = (get(handles.frame_slider,'Max')-1):-1:2
    for i = 1:length(handles.Region)
        for j = 1:length(handles.Region(i).Fascicle)
            % state estimation
            handles = state_estimator(handles,f,f+1,'backward');
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

% a posteriori variance estimate
K.P_minus = k.P_prev + k.Q;

% kalman xshiftcor
K.K = K.P_minus / (K.P_minus + k.R);

if (K.K < 0) || (K.K > 1)
    disp('warning: kalman gain outside 0-1 interval');

    K.K(K.K<0) = 0;
    K.K(K.K>1) = 1;
end

% estimated state
K.x_plus = k.x_minus + K.K * (k.y - k.x_minus);

% estimated variance
K.P_plus = (1-K.K) * K.P_minus;


function[handles] = state_estimator(handles,frame_no,prev_frame_no, direction)

i = 1; j = 1;

% define the initial variance
handles.Region(1).ROIp{1} = zeros(1,5);
handles.Region(i).Fascicle(j).fas_p{1} = zeros(1,2);

%% Apply warp
fas_prev = [handles.Region(i).Fascicle(j).fas_x{prev_frame_no}' handles.Region(i).Fascicle(j).fas_y{prev_frame_no}'];
alpha_prev = handles.Region(i).Fascicle(j).alpha{prev_frame_no};

if strcmp(direction,'forward')
    w = handles.Region(i).warp(:,:,prev_frame_no);
    fas_new = transformPointsForward(w, fas_prev);
else
    w = handles.Region(i).warp(:,:,frame_no);
    fas_new = transformPointsInverse(w, fas_prev);
end

% Estimate the change in fascicle angle from the change in points
dalpha = abs(atan2d(diff(fas_new(:,2)), diff(fas_new(:,1)))) - abs(atan2d(diff(fas_prev(:,2)), diff(fas_prev(:,1))));
alpha_new = alpha_prev + dalpha;

% A priori state estimate
x_minus = [fas_new(2,1) alpha_new];

% fit the current aponeurosis
ROI         = [handles.Region(i).ROIx{frame_no} handles.Region(i).ROIy{frame_no}];
super_apo   = ROI([1,4],:);
deep_apo    = ROI([2,3],:);
super_coef  = polyfit(super_apo(:,1), super_apo(:,2), 1);
deep_coef   = polyfit(deep_apo(:,1), deep_apo(:,2), 1);

%% State estimation superficial aponeurosis attachment
% previous estimate covariance
P_prev = handles.Region(i).Fascicle(j).fas_p{prev_frame_no};

% get the process noise and measurement noise covariance
dx = sqrt((fas_new(2,1)-fas_prev(2,1)).^2 + (fas_new(2,2)-fas_prev(2,2)).^2);
s.Q = getQ(handles, dx);
s.R = handles.X;

% a priori estimate from optical flow
s.x_minus = x_minus(1);

% 'measurement', here is the first value
s.y = handles.Region(i).Fascicle(j).fas_x{1}(2);

% previous state covariance
s.P_prev = P_prev(1);

% run kalman filter
S = run_kalman_filter(s);

% ascribe
fasx2 = S.x_plus;
supP = S.P_plus;

if isinf(handles.Q)
    fasx2 = s.y;
    supP = s.R;
elseif isinf(handles.X)
    fasx2 = s.x_minus;
    supP = s.Q;
end

% get the vertical point from the estimated aponeurosis
fasy2 = super_coef(2) + fasx2*super_coef(1);

%% Fascicle angle estimate
% get the process noise and measurement noise covariance
dx = abs(dalpha);
dx(dx<0.005) = 0;
f.Q = getQ(handles, dx);
f.R = handles.R(1);

% apriori estimate from optical flow
f.x_minus = x_minus(2);

% measurement from Hough transform
f.y = handles.geofeatures(frame_no).alpha;

% previous state covariance
f.P_prev = P_prev(2);

% run kalman filter
F = run_kalman_filter(f);

% ascribe
alpha = F.x_plus;
Kgain = F.K;
fasP = F.P_plus;

if isinf(handles.Q)
    alpha = f.y;
    Kgain = 1;
    fasP = f.R;
end

% get the deep attachment point from the superficial point and the angle
fas_coef(1) = -tand(alpha);
fas_coef(2) =  fasy2 - fas_coef(1) * fasx2;
fasx1 = (fas_coef(2) - deep_coef(2)) / (deep_coef(1) - fas_coef(1));
fasy1 = deep_coef(2) + fasx1*deep_coef(1);

%% update
% state and dependent variables
handles.Region(i).Fascicle(j).alpha{frame_no}   = alpha;
handles.Region(i).Fascicle(j).fas_x{frame_no}   = [fasx1 fasx2];
handles.Region(i).Fascicle(j).fas_y{frame_no}   = [fasy1 fasy2];

% state covariance
handles.Region(i).Fascicle(j).fas_p{frame_no} = [supP fasP];

% kalman xshiftcor for fascicle
handles.Region(i).Fascicle(j).K(frame_no) = Kgain;

% calculate the length and pennation for the current frame
handles.Region(i).fas_pen(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y{frame_no})),...
    abs(diff(handles.Region(i).Fascicle(j).fas_x{frame_no})));

scalar = handles.ID;%/handles.vidHeight;
scalar = handles.ID/handles.vidHeight;

handles.Region(i).fas_pen(frame_no,j) = handles.Region(i).Fascicle(j).alpha{frame_no}/180*pi;

handles.Region(i).fas_length(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{frame_no}).^2 +...
    diff(handles.Region(i).Fascicle(j).fas_x{frame_no}).^2);


% --- Executes on button press in Auto_Detect.
function [handles] = Auto_Detect_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_Detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% find current frame number from slider
frame_no = round(get(handles.frame_slider,'Value'));

w = handles.vidWidth;

if ~isvalid(handles.S)
    handles.S = images.roi.Rectangle(gca,'position', [1 handles.parms.apo.super.cut(1)*handles.vidHeight w diff(handles.parms.apo.super.cut)*handles.vidHeight],'color','blue');
end

if ~isvalid(handles.D)
    handles.D = images.roi.Rectangle(gca,'position', [1 handles.parms.apo.deep.cut(1)*handles.vidHeight w diff(handles.parms.apo.deep.cut)*handles.vidHeight],'color','green');
end

%% Aponeurosis detection
axes(handles.axes1); 

handles.parms.apo.super.cut = [handles.S.Position(2) handles.S.Position(2)+handles.S.Position(4)] / handles.vidHeight;
handles.parms.apo.deep.cut = [handles.D.Position(2) handles.D.Position(2)+handles.D.Position(4)] / handles.vidHeight;

set(handles.S, 'EdgeAlpha',0,'FaceAlpha',0.1,'InteractionsAllowed','none')
set(handles.D, 'EdgeAlpha',0,'FaceAlpha',0.1,'InteractionsAllowed','none')

% max range
% parms.fas.range = 90 - [-90 89];

% % detect orientation
data = imresize(handles.ImStack(:,:,frame_no+handles.start_frame-1), 1/handles.imresize_fac);
% [geofeatures, ~, parms] = auto_ultrasound(data, parms);
% alphas = geofeatures.alphas;
% alphas(alphas==90 | alphas == 180) = [];
% % 
% if median(alphas) < 90 || isempty(alphas)
%     parms.fas.range = [5 80];
% else
%     parms.fas.range = [5 80] + 90;
% end

% run TimTrack
geofeatures = auto_ultrasound(data, handles.parms);

% scale
geofeatures.thickness = geofeatures.thickness * handles.imresize_fac;
geofeatures.super_coef(2) = geofeatures.super_coef(2)     .* [handles.imresize_fac];
geofeatures.deep_coef(2) = geofeatures.deep_coef(2)       .* [handles.imresize_fac];
geofeatures.fas_coef(2) = geofeatures.fas_coef(2)         .* [handles.imresize_fac];

n = handles.vidWidth;
i = 1; j = 1;

Deep_intersect_x = round((geofeatures.deep_coef(2) - geofeatures.fas_coef(2))   ./ (geofeatures.fas_coef(1) - geofeatures.deep_coef(1)));
Super_intersect_x = round((geofeatures.super_coef(2) - geofeatures.fas_coef(2)) ./ (geofeatures.fas_coef(1) - geofeatures.super_coef(1)));
Super_intersect_y = polyval(geofeatures.super_coef, Super_intersect_x);
Deep_intersect_y = polyval(geofeatures.deep_coef, Deep_intersect_x);

% draw a fascicle
% h = drawline('Position', [Deep_intersect_x Deep_intersect_y; Super_intersect_x Super_intersect_y], 'color', 'red', 'linewidth',2);

handles.Region.Fascicle.fas_x{frame_no} = [Deep_intersect_x Super_intersect_x];
handles.Region.Fascicle.fas_y{frame_no} = [Deep_intersect_y Super_intersect_y];

handles.Region.Fascicle.fas_x_original{frame_no} = handles.Region.Fascicle.fas_x{frame_no};
handles.Region.Fascicle.fas_y_original{frame_no} = handles.Region.Fascicle.fas_y{frame_no};

handles.Region(i).ROIx{frame_no} = [1 1 n n 1]';
handles.Region(i).ROIy{frame_no} = round([polyval(geofeatures.super_coef, 1) polyval(geofeatures.deep_coef, [1 n]) polyval(geofeatures.super_coef, [n 1])])';

for i = 1:size(handles.Region,2)
    handles.Region(i).Fascicle(j).current_xy = [handles.Region(i).Fascicle(j).fas_x{frame_no};handles.Region(i).Fascicle(j).fas_y{frame_no}]';
    handles.Region(i).fas_pen(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y{frame_no})),abs(diff(handles.Region(i).Fascicle(j).fas_x{frame_no})));

    if i == 1
        scalar = handles.ID/handles.vidHeight;
        handles.CIm = handles.Im;
    end

    handles.Region(i).fas_length(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{frame_no}).^2 + diff(handles.Region(i).Fascicle(j).fas_x{frame_no}).^2);

    if ~isfield(handles.Region(i).Fascicle(j),'analysed_frames')
        handles.Region(i).Fascicle(j).analysed_frames = frame_no;
    else handles.Region(i).Fascicle(j).analysed_frames = sort([handles.Region(i).Fascicle(j).analysed_frames frame_no]);
    end

    handles.Region(i).Fascicle(j).alpha{frame_no} = handles.Region(i).fas_pen(frame_no,j) * 180/pi;
end

% Update handles structure
guidata(hObject, handles);

show_image(hObject,handles);
show_data(hObject, handles)

function xshiftcor_Callback(hObject, eventdata, handles)
% hObject    handle to xshiftcor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xshiftcor as text
%        str2double(get(hObject,'String')) returns contents of xshiftcor as a double

handles.fcor(3) = str2double(get(hObject,'String')); % [Hz]

% Update handles structure
guidata(hObject, handles);

% If we already processed, run state estimation
if isfield(handles.Region,'Fascicle')
    if length(handles.Region.Fascicle.fas_x) == handles.NumFrames
    do_state_estimation(hObject, eventdata, handles)
    end
end

% --- Executes during object creation, after setting all properties.
function xshiftcor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xshiftcor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.fcor(3) = str2double(get(hObject,'String')); % [Hz]

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
handles.flipimage = get(hObject,'Value');
if isfield(handles,"Region")
    updateX = @(fas_x) flip(handles.vidWidth - fas_x); %only here we need correction as axis starts from 1 (plotting)
   
    for i = 1:numel(handles.Region)
        %adjust ROI coordinates and logic mask image
        handles.Region(i).ROIx = cellfun(updateX, handles.Region(i).ROIx, 'UniformOutput', false);
        handles.Region(i).ROIy = cellfun(@flip, handles.Region(i).ROIy, 'UniformOutput', false);
        
        if isfield(handles.Region(i),"ROI")
            handles.Region(i).ROI = cellfun(@fliplr, handles.Region(i).ROI, 'UniformOutput', false);
        end
        %adjust each fascicle's pts
        for j = 1:numel(handles.Region(i).Fascicle)
            handles.Region(i).Fascicle(j).fas_x = cellfun(updateX, handles.Region(i).Fascicle(j).fas_x, 'UniformOutput', false);
            handles.Region(i).Fascicle(j).fas_y = cellfun(@flip, handles.Region(i).Fascicle(j).fas_y, 'UniformOutput', false);
        end
    end

end

% If statement not necessary, if tick flip else flip back, so everytime flipimage
% changes which depends on the callback, flip the image
%if handles.flipimage 

if isfield(handles, 'ImStack')
    handles.ImStack = flip(handles.ImStack, 2);
end

%end
guidata(hObject, handles);
show_image(hObject,handles);


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


% --- Function for detecting borders of the images and returning the rect matrix to crop.
function ImCropped = autocrop_Tim(Im)

dIm = sum(abs(diff(Im,1, 3)),3);

ndIm = dIm / max(dIm(:));

th = .05;
ndIm(ndIm>th) = 1;
ndIm(ndIm<=th)= 0;

BW2 = bwareaopen(ndIm,50);
% 
% figure(10)
% imshow(BW2, [])

% Find connected components in the filtered matrix.
stats = regionprops(BW2, 'BoundingBox');

% Extract the bounding box information.
boundingBoxes = cat(1, stats.BoundingBox);

% Calculate the overall bounding box that encompasses all smaller bounding boxes.
B = round([min(boundingBoxes(:, 1)), min(boundingBoxes(:, 2)), ...
    max(boundingBoxes(:, 1) + boundingBoxes(:, 3)) - min(boundingBoxes(:, 1)), ...
    max(boundingBoxes(:, 2) + boundingBoxes(:, 4)) - min(boundingBoxes(:, 2))]);

ImCropped = Im(B(2):(B(2)+B(4)-1), B(1):(B(1)+B(3)-1),:);


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
files= dir(dir_data);

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
clearvars -except files eventdata hObject handles

% for j = 1:length(subfolders)
%     foldername = subfolders{j};
%
%     cd([mainfolder foldername]);
%
%     files = dir('*.mp4');
%

for k = 1:numel(files) %foreach file

    handles.fname = files(k).name;
    handles.pname = [files(k).folder,'/'];
    cd(handles.pname)
    handles.movObj = VideoReader([handles.pname handles.fname]);
    mb = waitbar(0,'Loading Video....');

    % get info
    handles.vidHeight = handles.movObj.Height;
    handles.vidWidth = handles.movObj.Width;
    handles.NumFrames = handles.movObj.NumFrames;
    handles.FrameRate = handles.movObj.FrameRate;

    i=1;

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
            handles.Time = TVD.TVDdata.Time(1:end-1); %-1 because last timestamp is repeated or UT doesn't load all of them
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

    %
    if handles.flipimage % get(hObject,'Value');
        handles.ImStack = flip(handles.ImStack, 2);
    end
    % show
    handles = show_image(hObject,handles);
    
    %create a subfolder where to save results
    if ~exist([handles.pname 'Tracked/'], 'dir')
        mkdir([handles.pname 'Tracked/'])
    end
    handles.pname = [handles.pname 'Tracked/'];

    for kk = 1:2 %
        if kk == 1
            Qs = nan;
            handles.Q = Qs;
        else
            Qs = [0, 10.^(-4:0), 1000, inf];
            handles.Q = 1;
        end

        % process
        handles = process_all_Callback(hObject, eventdata, handles);

        for i = 1:length(Qs)
            handles.Q = Qs(i);

            handles = do_state_estimation(hObject, eventdata, handles);

            % save
            save_video_Callback(hObject, eventdata, handles)
            Save_As_Mat_Callback(hObject, eventdata, handles)

        end
    end
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
    handles = lowpass_ROI(hObject, eventdata, handles);
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

% Update handles structure
guidata(hObject, handles);

% If we have estimates, run state estimation
if isfield(handles, 'Region')
    do_state_estimation(hObject, eventdata, handles)
end

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

handles.Q = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% If we have estimates, run state estimation
if isfield(handles, 'Region')
    do_state_estimation(hObject, eventdata, handles)
end

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

handles.Q = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);
