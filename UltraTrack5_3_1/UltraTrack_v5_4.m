function varargout = UltraTrack_v5_4(varargin)
% ULTRATRACK_V5_4 M-file for UltraTrack_v5_4.fig
%      ULTRATRACK_V5_4, by itself, creates a new ULTRATRACK_V5_4 or raises the existing
%      singleton*.
%      H = ULTRATRACK_V5_4 returns the handle to a new ULTRATRACK_V5_4 or the handle to
%      the existing singleton*.
%
%      ULTRATRACK_V5_4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ULTRATRACK_V5_4.M with the given input arguments.
%
%      ULTRATRACK_V5_4('Property','Value',...) creates a new ULTRATRACK_V5_4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UltraTrack_v5_4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UltraTrack_v5_4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UltraTrack_v5_4

% Last Modified by GUIDE v2.5 22-Mar-2024 12:57:53
% Last Modified by Paolo Tecchio 17/08/2022
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UltraTrack_v5_4_OpeningFcn, ...
    'gui_OutputFcn',  @UltraTrack_v5_4_OutputFcn, ...
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

% --- Executes just before UltraTrack_v5_4 is made visible.
function UltraTrack_v5_4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UltraTrack_v5_4 (see VARARGIN)

% Check if Parallel Computing Toolbox is installed
% chkParallelToolBox(); %if exists it runs infinitely till Ultratimtrack is closed

% Choose default command line output for UltraTrack_v5_4
handles.output = hObject;

%add automatically all files and subfolders dynamically
filename = [mfilename,'.m'];
fullpath = which(filename);
foldername = erase(fullpath,filename);
mainfoldername = erase(foldername,'UltraTrack5_3_1/'); % need to move one up to include the submodule in path
addpath(genpath(mainfoldername));

% check to make sure there is a settings mat-file present, if not then make
% one in the directory where the m-file sits.
[program_directory, ~, ~] = fileparts(mfilename('fullpath'));

if ismac || isunix
    if isempty(dir([program_directory '/ultrasound_tracking_settings.mat']))
        ImageDepth = 65;
        Sigma = 3;
        S_Step = 3;
        Position = get(gcf,'Position');
        Default_Directory = cd;
        save([program_directory '/ultrasound_tracking_settings.mat'], 'ImageDepth', 'Sigma',...
            'S_Step', 'Position', 'Default_Directory');
    else
        load ultrasound_tracking_settings.mat
    end
end
%
if ispc
    if isempty(dir([program_directory '\ultrasound_tracking_settings.mat']))
        ImageDepth = 65;
        Sigma = 3;
        S_Step = 3;
        Position = get(gcf,'Position');
        Default_Directory = cd;
        save([program_directory '\ultrasound_tracking_settings.mat'], 'ImageDepth', 'Sigma',...
            'S_Step', 'Position', 'Default_Directory');
    else
        load ultrasound_tracking_settings.mat
    end

end

% load the settings mat-file and set default settings
handles.ID = ImageDepth;
set(handles.ImDepthEdit,'String',num2str(ImageDepth));
handles.SIGMA = Sigma;
handles.S_STEP = S_Step;
set(0,'RecursionLimit',3000)
set(gcf,'DoubleBuffer','on','Position',Position);
%cd(Default_Directory)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UltraTrack_v5_4 wait for user response (see UIRESUME)
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

            handles.vidHeight = handles.movObj.Height;
            handles.vidWidth = handles.movObj.Width;
            i=1;
            while hasFrame(handles.movObj)
                waitbar(handles.movObj.CurrentTime/handles.movObj.Duration,mb)
                if regexp(handles.movObj.VideoFormat,'RGB')
                    %handles.ImStack(:,:,i) = fliplr(im2gray(readFrame(handles.movObj)));
                    handles.ImStack(:,:,i) = im2gray(readFrame(handles.movObj));
                else
                    %handles.ImStack(:,:,i) = fliplr(readFrame(handles.movObj));
                    handles.ImStack(:,:,i) = readFrame(handles.movObj);
                end
                i=i+1;
            end
            handles.NumFrames = size(handles.ImStack,3);
            handles.FrameRate = handles.NumFrames/handles.movObj.Duration;

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
    load([path '/' name '.mat']);
    if isfield(TVDdata,'cmPerPixY') %check whether the field exists and update scalar
        ImageDepth = round(TVDdata.cmPerPixY*10,3); %round to 3 digits
        handles.ID = ImageDepth;
        set(handles.ImDepthEdit,'String',num2str(ImageDepth));

        % Update handles structure
        guidata(hObject, handles);
    end
    clearvars TVDdata path name
end

% display the path and name of the file in the filename text box
set(handles.filename,'String',[handles.pname handles.fname])

% set the string in the frame_number box to the current frame value (1)
set(handles.frame_number,'String',num2str(1))

handles.start_frame = 1;

% set the limits on the slider - use number of frames to set maximum (min =
% 1)
set(handles.frame_slider,'Min',1);
set(handles.frame_slider,'Max',handles.NumFrames);
set(handles.frame_slider,'Value',1);
set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 10/handles.NumFrames]);


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
elseif exist("TVDdata",'var') %if the mat file associated with mp4 exists, load that because echowave has unconstant framerate
    handles.Time = TVDdata.Time(1:end-1); %-1 because last timestamp is repeated or UT doesn't load all of them?
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

show_image(hObject,handles)
waitbar(1,mb)
close(mb)

% --------------------------------------------------------------------
function Load_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[dfile,dpath] = uigetfile('*.*','Select Associated Data File');
handles.dat_file = [dpath dfile];

% determine the file prefix and filetype
[~,handles.file,handles.ext] = fileparts(handles.dat_file);

% open the file using relevant method depending on filetype
switch handles.ext

    case '.mat'
        % if it is a matfile then load
        data = load(handles.dat_file);

        % find the variables in the matfile
        s = whos('-file', handles.dat_file);

        for i = 1:length(s)
            if strcmp(s(i).class,'struct')

                % if data is analog data then make time and value column
                if isfield(data.(s(i).name),'values')
                    D(i).type = 'Analog';
                    D(i).name = data.(s(i).name).title;
                    D(i).time = data.(s(i).name).start:data.(s(i).name).interval:data.(s(i).name).start+(data.(s(i).name).interval*(data.(s(i).name).length-1));
                    D(i).value = data.(s(i).name).values;

                    % if data is analog data then make time and value column
                elseif isfield(data.(s(i).name),'times')
                    D(i).type = 'Digital';
                    D(i).name = data.(s(i).name).title;
                    D(i).time =  [(data.(s(i).name).times)'; (data.(s(i).name).times)'];
                    D(i).value =  [(zeros(size(data.(s(i).name).times)))';(ones(size(data.(s(i).name).times))*5)'];

                else warning('The file contains structures which cannot be resolved');
                    D(i).type = 'Unknown';
                    D(i).name = 'Unknown';
                    D(i).time = [];
                    D(i).value = [];
                end

            end
        end

    case '.c3d'

        % at this stage we will limit this data to marker, analog and
        % forceplate (calibrated analog) data. We can also add angles etc.
        D = [];
        % if file is c3d then load using the BTK matlab wrapper
        acq = btkReadAcquisition(handles.dat_file);

        % get marker data
        [markers, markersInfo] = btkGetMarkers(acq);
        % find the names of the markers from the field names of the
        % structure
        marker_names = fieldnames(markers);
        % if there are marker names then store to structure
        if ~isempty(marker_names)
            for i = 1:length(marker_names)
                M(i).type = 'Marker';
                M(i).name = ['Marker ' marker_names{i}];
                M(i).time = 1/markersInfo.frequency:1/markersInfo.frequency:length(markers.(marker_names{i}))/markersInfo.frequency;
                M(i).value = markers.(marker_names{i});
            end
            if isempty(D)
                D = M;
            else D(end+1:end+length(M)) = M;
            end
        end

        % get analog data
        [analogs, analogsInfo] = btkGetAnalogs(acq);
        % find the names of the analog channels from the field names of the
        % structure
        analog_names = fieldnames(analogs);
        % if there are analog channel names then store to structure
        if ~isempty(analog_names)
            for i = 1:length(analog_names)
                A(i).type = 'Analog';
                A(i).name = ['Analog ' analog_names{i}];
                A(i).time = 1/analogsInfo.frequency:1/analogsInfo.frequency:length(analogs.(analog_names{i}))/analogsInfo.frequency;
                A(i).value = analogs.(analog_names{i});
            end
            if isempty(D)
                D = A;
            else D(end+1:end+length(A)) = A;
            end
        end

    case '.dat'
        % load the biodex dat file
        biodex = biodex_load_file(handles.dat_file);

        timeDiff = biodex.data(end,1) - handles.TimeStamps(end)

        for i = 2:size(biodex.data,2)
            D(i-1).type = 'Biodex_Analog';
            switch i
                case 2
                    D(i-1).name = 'Torque';
                case 3
                    D(i-1).name = 'Velocity';
                case 4
                    D(i-1).name = 'Position';
                otherwise
                    D(i-1).name = ['Channel_' num2str(i-1)];
            end
            D(i-1).time =  biodex.data(:,1)-biodex.data(1,1);
            D(i-1).value =  biodex.data(:,i);
        end

    case '.txt'
        % load an ascii file --> must ensure that the first column is a
        % time column
        [~, data] = hdrload(handles.dat_file);

        for i = 2:size(data,2)
            D(i-1).type = 'Ascii_Analog';
            D(i-1).name = ['Channel_' num2str(i-1)];
            D(i-1).time =  data(:,1);
            D(i-1).value =  data(:,i);
        end

    otherwise
        warndlg('Unrecognised file extension - cannot load')
        return
end
% save name of variable to variable list
handles.Avail_Variables = {D.name};

% fill in the variable list for plotting
set(handles.variable_list,'String',handles.Avail_Variables);

handles.D = D;
show_image(hObject,handles)
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = UltraTrack_v5_4_OutputFcn(~, eventdata, handles)
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

    set(handles.frame_slider,'Min',1);
    set(handles.frame_slider,'Max',handles.NumFrames);
    set(handles.frame_slider,'Value',1);
    set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 5/handles.NumFrames]);

    % set the string in the frame_number to 1
    set(handles.frame_number,'String',1);

    % update the image axes using show_image function (bottom)
    show_image(hObject,handles)
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
    show_image(hObject,handles)
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
show_image(hObject,handles)

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
show_image(hObject,handles)

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
    show_image(hObject,handles)

end


% --- Executes on button press in define_roi.
function define_roi_Callback(hObject, eventdata, handles)
% hObject    handle to define_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    axes(handles.axes1)

    frame_no = round(get(handles.frame_slider,'Value'));

    i = get(handles.no_tracked_regions,'Value');
    [handles.Region(i).ROI{frame_no},handles.Region(i).ROIx{frame_no},...
        handles.Region(i).ROIy{frame_no}] = roipoly;
    handles.CIm = handles.Im;

    show_image(hObject,handles)

end

% --- Executes on button press in define_fascicle.
function define_fascicle_Callback(hObject, eventdata, handles)
% hObject    handle to define_fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'));

    i = get(handles.no_tracked_regions,'Value');
    j = get(handles.no_tracked_fascicles,'Value');


    % select the two points that define the fascicle line of action
    %set(handles.axes1,'XLimMode','auto')
    [p1,p2]=rbline(handles.axes1);
    handles.Region(i).Fascicle(j).fas_x{frame_no}(1) = p1(1);
    handles.Region(i).Fascicle(j).fas_x{frame_no}(2) = p2(1);
    handles.Region(i).Fascicle(j).fas_y{frame_no}(1) = p1(2);
    handles.Region(i).Fascicle(j).fas_y{frame_no}(2) = p2(2);

    % %%adjust this ???
    if handles.Region(i).Fascicle(j).fas_y{frame_no}(1) < handles.Region(i).Fascicle(j).fas_y{frame_no}(2)
        handles.Region(i).Fascicle(j).fas_y{frame_no} = fliplr(handles.Region(i).Fascicle(j).fas_y{frame_no});
        handles.Region(i).Fascicle(j).fas_x{frame_no} = fliplr(handles.Region(i).Fascicle(j).fas_x{frame_no});
    end

    handles.Region(i).Fascicle(j).current_xy = [handles.Region(i).Fascicle(j).fas_x{frame_no};handles.Region(i).Fascicle(j).fas_y{frame_no}]';

    % calculate the length and pennation for the current frame
    handles.Region(i).fas_pen(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y{frame_no})),abs(diff(handles.Region(i).Fascicle(j).fas_x{frame_no})));
    scalar = handles.ID;%/handles.vidHeight;
    handles.Region(i).fas_length(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{frame_no}).^2 + diff(handles.Region(i).Fascicle(j).fas_x{frame_no}).^2);

    handles.CIm = handles.Im;

    if ~isfield(handles.Region(i).Fascicle(j),'analysed_frames')
        handles.Region(i).Fascicle(j).analysed_frames = frame_no;
    else 
        handles.Region(i).Fascicle(j).analysed_frames = sort([handles.Region(i).Fascicle(j).analysed_frames frame_no]);
    end

    show_image(hObject,handles)

end

% --------------------------------------------------------------------
function menu_define_fascicle_Callback(hObject, eventdata, handles)
% hObject    handle to menu_define_fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

define_fascicle_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_define_roi_Callback(hObject, eventdata, handles)
% hObject    handle to menu_define_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

define_roi_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_edit_current_Callback(hObject, eventdata, handles)
% hObject    handle to menu_process_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clear_fascicle.
function clear_fascicle_Callback(hObject, eventdata, handles)
% hObject    handle to clear_fascicle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')
    if isfield(handles,'Region')
        handles = rmfield(handles,'Region');
    end

    %     if ~isempty(get(handles.keyframe_list,'String'))
    %         set(handles.keyframe_list,'String',[])
    %     end

    set(handles.no_tracked_regions,'Value',1,'String','1')
    set(handles.no_tracked_fascicles,'Value',1,'String','1')
    %     set(handles.RegionsToCorrect,'String',{'all','1'},'Value',1);
    %     set(handles.FasciclesToCorrect,'String',{'all','1'},'Value',1);


    % set current frame to 1
    set(handles.frame_slider,'Value',1);
    set(handles.frame_number,'String',1);

    cla(handles.length_plot) %clean fascicle length data
    cla(handles.mat_plot)%clean fascicle angle data
    
    cla(handles.axes1) %clean image data    
    set(handles.axes1)
    axis tight %first center the image
    axis auto 
    show_image(hObject,handles)
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
Sigma = handles.SIGMA;
S_Step = handles.S_STEP;
Position = get(gcf,'Position');
Default_Directory = cd;

[directory, ~, ~] = fileparts(mfilename('fullpath'));

if ismac || isunix
    save([directory '/ultrasound_tracking_settings.mat'], 'ImageDepth', 'Sigma',...
        'S_Step', 'Position', 'Default_Directory');
end

if ispc
    save([directory '\ultrasound_tracking_settings.mat'], 'ImageDepth', 'Sigma',...
        'S_Step', 'Position', 'Default_Directory');
end

msgbox('Settings Saved')

% --------------------------------------------------------------------
function menu_affine_settings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_affine_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'SIGMA')
    handles.SIGMA = 3;
end

if ~isfield(handles, 'S_STEP')
    handles.SIGMA = 3;
end

A = inputdlg({'Sigma','Sample Step'}, 'Affine Flow Settings', 1,...
    {num2str(handles.SIGMA),num2str(handles.S_STEP)});
handles.SIGMA = str2double(A{1});
handles.S_STEP = str2double(A{2});

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_clear_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to menu_clear_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    if isfield(handles,'Region')
        handles = rmfield(handles,'Region');
    end


    if ~isempty(get(handles.keyframe_list,'String'))
        set(handles.keyframe_list,'String',[])
    end

    % set current frame to 1
    set(handles.frame_slider,'Value',1);
    set(handles.frame_number,'String',1);

    cla(handles.length_plot)
    cla(handles.mat_plot)
    cla(handles.axes1) %clean image data    

    show_image(hObject,handles)

    % Correct the drop down lists for regions and fascicles to match
    % imported data

    ROIlistString = {'1'};
    ROItoCorrString = {'all','1'};
    set(handles.no_tracked_regions,'String',ROIlistString);
    set(handles.RegionsToCorrect,'String',ROItoCorrString);

    FASlistString = {'1'};
    FAStoCorrString = {'all','1'};
    set(handles.no_tracked_fascicles,'String',FASlistString);
    set(handles.FasciclesToCorrect,'String',FAStoCorrString);
end


% --------------------------------------------------------------------
function AutoCrop_Callback(hObject, eventdata, handles)
% hObject    handle to AutoCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Im = handles.ImStack(:,:,1);

handles.vidHeight = size(Im,1);

%corners = detectAutoCrop(handles.movObj);
handles.crop_rect = detectAutoCrop2(handles.ImStack,handles.movObj.FrameRate/2);
handles.ImStackOr = handles.ImStack; %back up original
    %tmp = zeros(size(handles.ImStack));

    %Crop all images before updating
    for ii = 1 : handles.NumFrames
        tmp(:,:,ii) = imcrop(handles.ImStack(:,:,ii),handles.crop_rect);
    end

    handles.ImStack = tmp; %overwrite the one used thourghout the script
    handles.vidHeight = handles.crop_rect(4);
    handles.vidWidth = handles.crop_rect(3);
    clearvars tmp
    % Clean axis from original image and tight axis on the cropped image
    cla
    % update the image axes using show_image function (bottom)
    show_image(hObject,handles)
    axis tight

%{ 
added by Tim
% gIm = imbinarize(imgradient(Im));
gIm = imbinarize(mean(diff(handles.ImStack,3),3));

% sum across depth
P = nan(2,2);
for j = 1

    if j == 2
        gIm = gIm(:,P(1,1):P(1,2));
    end

    sgIm = sum(gIm,j);

    % scale between -1 and 1
    y = 2*sgIm / max(sgIm) - 1;


    %     y = sgIm - mean(sgIm);

    % fit
    x = 1:length(y);

    if j == 1
        costfun = @(p,x,y) -sum(((x(:) * -1) + 2 * (x(:) > p(1) & x(:) < p(2))) .* y(:));
        P(j,:) = fminsearch(@(p) costfun(p, x(:), y(:)), [length(y)/2 length(y)/2]);
    else
        costfun = @(p,x,y) -sum(((x(:) * -1) + 2 * (x(:) > p(1))) .* y(:));
        P(j,1) = fminsearch(@(p) costfun(p, x(:), y(:)), 50);
    end

end


handles.crop_rect(1) = round(P(1,1)) + 2; % add two to be safe
handles.crop_rect(3) = round(P(1,2)-P(1,1)) - 4; % add two to be safe

if j == 1
    handles.crop_rect(2) = 1;
    handles.crop_rect(4) = handles.vidHeight;
else

    handles.crop_rect(2) = round(P(2,1)) + 2;
    handles.crop_rect(4) = handles.vidHeight - handles.crop_rect(2);
end

handles.vidHeight = handles.crop_rect(4);
handles.vidWidth = handles.crop_rect(3);

%}

%%

% % work out how to trim the edges of the ultrasound (edges are intensity value of
% % 56]
% m = find(Im(:,1) == 56); % find where the image starts vertically - at top of image the gray changes from 66 to 56.
% n = find(mean((Im(end-50:end,:) == 56))<0.4); % use the averaage value of the bottom 50 rows to determine whether this is a high value (i.e. gray) or low (i.e. real image)
%
% if isempty(m)
%     return;
% end
% if isempty(n)
%     return;
% end
%
% handles.crop_rect = [n(find(diff(n)>1)+1) m(1) n(end)-n(find(diff(n)>1)+1) length(m)-1]; % define the new rectangle for autocroping
% handles.vidHeight = handles.crop_rect(4);
% handles.vidWidth = handles.crop_rect(3);
% Update handles structure
guidata(hObject, handles);
% update the image axes using show_image function (bottom)
show_image(hObject,handles)

% --------------------------------------------------------------------
function menu_crop_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_crop_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    % define current axes
    set(handles.axes1)
    %axes(handles.axes1)

    % use imcrop tool to determine croppable area
    [~,handles.crop_rect] = imcrop(handles.axes1);

    handles.crop_rect = round(handles.crop_rect);
    handles.vidHeight = handles.crop_rect(4);
    handles.vidWidth = handles.crop_rect(3);
    
    handles.ImStackOr = handles.ImStack; %back up original
    %tmp = zeros(size(handles.ImStack));

    %Crop all images before updating
    for ii = 1 : handles.NumFrames
        tmp(:,:,ii) = imcrop(handles.ImStack(:,:,ii),handles.crop_rect);
    end    
   
    handles.ImStack = tmp; %overwrite the one used thourghout the script
    clearvars tmp
    % Clean axis from original image and tight axis on the cropped image
    cla
    % update the image axes using show_image function (bottom)
    show_image(hObject,handles)
    axis tight

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
    show_image(hObject,handles)

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

% --- Executes on button press in export_check.
function export_check_Callback(hObject, eventdata, handles)
% hObject    handle to export_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of export_check



% --------------------------------------------------------------------
function menu_save_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'Region')

    for i = 1:length(handles.Region)

        if isfield(handles,'ImStack') && isfield(handles.Region(i),'fas_length')

            dt = 1/handles.FrameRate;
            %time = (dt:dt:length(handles.Region(i).fas_length)*dt)+((handles.start_frame-1)*dt);
            time = handles.Time(handles.start_frame:handles.start_frame + length(handles.Region(i).fas_length)-1);


            %determine any non-zero entries in fascicle length array
            nz = logical(handles.Region(i).fas_length(:,1) ~= 0);

            T = time(nz)';

            if isfield(handles.Region(i),'fas_length_corr')
                R(i).FL = handles.Region(i).fas_length_corr(nz,:)';
                R(i).PEN = handles.Region(i).fas_pen_corr(nz,:)';
            else R(i).FL = handles.Region(i).fas_length(nz,:)';
                R(i).PEN = handles.Region(i).fas_pen(nz,:)';
            end

            [~,handles.file,handles.ext] = fileparts(handles.dat_file);

            if strcmp(handles.ext,'.c3d')

                % Create a data acquisition object from the c3d file and store
                % a handle to it
                if i==1
                    c3d_handle = btkReadAcquisition(handles.dat_file);

                    Fq = btkGetPointFrequency(c3d_handle);
                    Fn = btkGetPointFrameNumber(c3d_handle);% Get the number of frames of point data
                    T_c3d = 0:1/Fq:(Fn-1)*(1/Fq);% create a time array at the rate of point data

                    % find the nearest times in the c3d file to the start and end
                    % of digitising
                    [~,SD] = min(abs(T_c3d-T(1)));
                    [~,ED] = min(abs(T_c3d-T(end)));
                end


                for j = 1:size(R(i).FL,1)

                    % create a vector which is the ultrasound fascicle length interpolated at
                    % the same time points that the marker data is determined
                    L_data = interp1(T,R(i).FL(j,:)',T_c3d(SD:ED)','linear','extrap');
                    L_data(isnan(L_data)) = 0;
                    P_data = interp1(T,R(i).PEN(j,:)',T_c3d(SD:ED)','linear','extrap');
                    P_data(isnan(P_data)) = 0;

                    FasDataOut = zeros(btkGetPointFrameNumber(c3d_handle),3);% create an appropriately sized array
                    FasDataOut(SD:ED,1) = L_data; % Insert the fascicle length data
                    PenDataOut = zeros(btkGetPointFrameNumber(c3d_handle),3);
                    PenDataOut(SD:ED,1) = P_data;

                    % Append the length data as point (type = scalar but, could be
                    % changed)

                    Points = btkGetPoints(c3d_handle);
                    Fdata_label = ['R' num2str(i) '_F' num2str(j) '_FasLength'];
                    Pdata_label = ['R' num2str(i) '_F' num2str(j) '_Pennation'];

                    if isfield(Points,Fdata_label)
                        btkSetPoint(c3d_handle,Fdata_label,FasDataOut);
                        btkSetPoint(c3d_handle,Pdata_label,PenDataOut);
                    else
                        btkAppendPoint(c3d_handle,'scalar',Fdata_label,FasDataOut);
                        btkAppendPoint(c3d_handle,'scalar',Pdata_label,PenDataOut);
                    end

                    clear L_data P_data FasDataOut PenDataOut Fdata_label Pdata_label
                end

            else

                Save_As_Txt_Callback(hObject,eventdata,handles);

            end

        end

    end
end

if strcmp(handles.ext,'.c3d')
    % create a new c3d file or overwrite the existing one
    fileout_suggest = [handles.pname handles.dat_file];
    [fileout, pathout] = uiputfile(fileout_suggest,'Save fascicle data as...');

    % write the data to the c3dfile specified
    cd(pathout);
    btkWriteAcquisition(c3d_handle,fileout);

end

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
            handles.Region(i).ROI{frame_no} = Fdat.Region(i).ROI;


        end
        Nregions = length(handles.Region);
        Nfas = max(Nfascicle);

        % Correct the drop down lists for regions and fascicles to match
        % imported data
        ROItoCorrString{1} = 'all';
        for i = 1:Nregions
            ROIlistString{i} = num2str(i);
            ROItoCorrString{i+1} = num2str(i);
        end
        set(handles.no_tracked_regions,'String',ROIlistString);
        set(handles.RegionsToCorrect,'String',ROItoCorrString);

        FAStoCorrString{1} = 'all';
        for i = 1:Nfas
            FASlistString{i} = num2str(i);
            FAStoCorrString{i+1} = num2str(i);
        end
        set(handles.no_tracked_fascicles,'String',FASlistString);
        set(handles.FasciclesToCorrect,'String',FAStoCorrString);
    else

        warndlg('Fascicle data not found','ERROR')

    end


    show_image(hObject,handles)

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

            % if fascicles have been corrected for drift then save this value
            % rather than original tracking
            if isfield(handles.Region(i).Fascicle(k),'fas_x_corr')

                Fdat.Region(i).Fascicle(k).fas_x = handles.Region(i).Fascicle(k).fas_x_corr{frame_no};
                Fdat.Region(i).Fascicle(k).fas_y = handles.Region(i).Fascicle(k).fas_y_corr{frame_no};

            elseif isfield(handles.Region(i).Fascicle(k),'fas_x')

                Fdat.Region(i).Fascicle(k).fas_x = handles.Region(i).Fascicle(k).fas_x{frame_no};
                Fdat.Region(i).Fascicle(k).fas_y = handles.Region(i).Fascicle(k).fas_y{frame_no};

            end

            % TimTrack_downsample
            Fdat.Region(i).Fascicle(k).fas_x_Hough = handles.Region(i).Fascicle(k).fas_x_Hough{frame_no};
            Fdat.Region(i).Fascicle(k).fas_y_Hough = handles.Region(i).Fascicle(k).fas_y_Hough{frame_no};


        end
        % also save the ROI data
        Fdat.Region(i).ROIx = handles.Region(i).ROIx{frame_no};
        Fdat.Region(i).ROIy = handles.Region(i).ROIy{frame_no};
        Fdat.Region(i).ROI = handles.Region(i).ROI{frame_no};
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

% determine maximal horizontal range
for i =  1:1:get(handles.frame_slider,'Max')
    x(i,:) = handles.Region.Fascicle.fas_x{i};
end

xrange = [min(x(:)) max(x(:))];

if isfield(handles,'ImStack')
    [fileout, pathout, FI] = uiputfile('*.mp4', 'Save video as');

    vidObj = VideoWriter([pathout fileout],'MPEG-4');
    vidObj.FrameRate = handles.FrameRate;
    open(vidObj);

    if FI > 0

        %         waitbar(0, 'Saving video ...');

        axes(handles.axes1)

        for i = 1:1:get(handles.frame_slider,'Max')

            set(handles.frame_slider,'Value',i);
            set(handles.frame_number,'String',num2str(i));

            show_image(hObject,handles)
            img = handles.axes1;
            set(img,'color','w')
            F = getframe(img);
            writeVideo(vidObj,F)
        end
    end
    close(vidObj)

end

% --- Executes on selection change in variable_list.
function variable_list_Callback(hObject, eventdata, handles)
% hObject    handle to variable_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns variable_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from variable_list
if isfield(handles,'D')
    show_image(hObject,handles)
end

% --- Executes during object creation, after setting all properties.
function variable_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to variable_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in legend_check.
function legend_check_Callback(hObject, eventdata, handles)
% hObject    handle to legend_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of legend_check
if isfield(handles,'D')
    show_image(hObject,handles)
end


% --------------------------------------------------------------------
function menu_process_all_Callback(hObject, eventdata, handles)
% hObject    handle to menu_process_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

process_all_Callback(hObject, eventdata, handles)


% --- Executes on selection change in no_tracked_fascicles.
function no_tracked_fascicles_Callback(hObject, eventdata, handles)
% hObject    handle to no_tracked_fascicles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns no_tracked_fascicles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from no_tracked_fascicles


% --- Executes during object creation, after setting all properties.
function no_tracked_fascicles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_tracked_fascicles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in no_tracked_regions.
function no_tracked_regions_Callback(hObject, eventdata, handles)
% hObject    handle to no_tracked_regions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns no_tracked_regions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from no_tracked_regions


% --- Executes during object creation, after setting all properties.
function no_tracked_regions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_tracked_regions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddRegionButton.
function AddRegionButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddRegionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Curr_ROI_list = get(handles.no_tracked_regions,'String');
if length(Curr_ROI_list)==1
    no_ROI = str2num(Curr_ROI_list);
else
    no_ROI = str2num(Curr_ROI_list{end});
end

add_ROI = [Curr_ROI_list; {num2str(no_ROI + 1)}];
set(handles.no_tracked_regions,'String',add_ROI);

Curr_ROIcorrect_list = get(handles.RegionsToCorrect,'String');
add_ROIcorr = [Curr_ROIcorrect_list; {num2str(no_ROI + 1)}];
set(handles.RegionsToCorrect,'String',add_ROIcorr)


% --- Executes on button press in AddFasButton.
function AddFasButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddFasButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Curr_FAS_list = get(handles.no_tracked_fascicles,'String');
if length(Curr_FAS_list) == 1
    no_FAS = str2num(Curr_FAS_list);
else
    no_FAS = str2num(Curr_FAS_list{end});
end
add_FAS = [Curr_FAS_list; {num2str(no_FAS + 1)}];
set(handles.no_tracked_fascicles,'String',add_FAS);

Curr_FAScorrect_list = get(handles.FasciclesToCorrect,'String');
add_FAScorr = [Curr_FAScorrect_list; {num2str(no_FAS + 1)}];
set(handles.FasciclesToCorrect,'String',add_FAScorr)


% --- Executes on selection change in RegionsToCorrect.
function RegionsToCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to RegionsToCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RegionsToCorrect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RegionsToCorrect


% --- Executes during object creation, after setting all properties.
function RegionsToCorrect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RegionsToCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FasciclesToCorrect.
function FasciclesToCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to FasciclesToCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FasciclesToCorrect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FasciclesToCorrect


% --- Executes during object creation, after setting all properties.
function FasciclesToCorrect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FasciclesToCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Save_Tracking_As_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Tracking_As (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Save_As_C3D_Callback(hObject, eventdata, handles)
% hObject    handle to Save_As_C3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if strcmp(handles.ext,'.c3d')
    menu_save_tracking_Callback(hObject,eventdata,handles)
else
    warndlg('DATA NOT SAVED - Saving to C3D format requires the input data file to be a C3D file','WARNING!')
end

% --------------------------------------------------------------------
function Save_As_Txt_Callback(hObject, eventdata, handles)
% hObject    handle to Save_As_Txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'Region')

    for i = 1:length(handles.Region)


        if isfield(handles,'movObj')||isfield(handles,'BIm')||isfield(handles,'ImStack') && isfield(handles.Region(i),'fas_length')

            dt = 1/handles.FrameRate;
            %time = (dt:dt:length(handles.Region(i).fas_length)*dt)+((handles.start_frame-1)*dt);
            time = handles.Time(handles.start_frame:handles.start_frame + length(handles.Region(i).fas_length)-1);

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

            col_head{i} = ['Time\tFascicle Length R' num2str(i) '_F1\tPennation Angle R' num2str(i) '_F1\t'];
            col_type{i} = '%3.4f\t%3.6f\t%3.6f';

            data_out{i} = [T R(i).FL(1,:)' R(i).PEN(1,:)'];

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

    if isfield(handles, 'D') && get(handles.export_check,'Value')
        % determine the variables highlighted in list
        cur_var = get(handles.variable_list,'Value');

        D_col_head = [];
        D_col_type = [];
        D_data_out = [];

        for i = 1:length(cur_var)

            new_col_head = [];
            new_col_type = [];


            if strcmp(handles.D(cur_var(i)).type,'Digital')
                new_col = zeros(size(T));
                for j = 1:size(handles.D(cur_var(i)).time,2)
                    new_col(find(abs(T-handles.D(cur_var(i)).time(1,j)) == ...
                        min(abs(T-handles.D(cur_var(i)).time(1,j)))),1) = ...
                        handles.D(cur_var(i)).value(2,j);
                end

                new_col_head = [new_col_head handles.D(cur_var(i)).name '\t'];
                new_col_type = [new_col_type '%3.6f\t'];
            else new_col = interp1(handles.D(cur_var(i)).time, ...
                    handles.D(cur_var(i)).value, T, 'linear', 'extrap');
                for j = 1:size(handles.D(cur_var(i)).value,2)
                    new_col_head = [new_col_head handles.D(cur_var(i)).name '\t'];
                    new_col_type = [new_col_type '%3.6f\t'];
                end
            end

            D_col_head = [D_col_head new_col_head];
            D_col_type = [D_col_type new_col_type];
            D_data_out = [D_data_out new_col];
        end
        header = [header D_col_head];
        type = [type D_col_type];
        dout = [dout D_data_out];
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

            dt = 1/handles.FrameRate;
            %time = (dt:dt:length(handles.Region(i).fas_length)*dt)+((handles.start_frame-1)*dt);
            time = handles.Time(handles.start_frame:handles.start_frame + length(handles.Region(i).fas_length)-1);
            %determine any non-zero entries in fascicle length array
            nz = logical(handles.Region(i).fas_length(:,1) ~= 0);

            T = time(nz)';

            TrackingData.Region = rmfield(handles.Region,'ROI');% remove the ROI to reduce file size
            TrackingData.start_frame = handles.start_frame;
            TrackingData.res = handles.ID;
            TrackingData.NumFrames = handles.NumFrames;
            %info about tracking for replication purposes
            TrackingData.ProcessingTime = handles.ProcessingTime; %two
            TrackingData.BlockSize = handles.BlockSize;
            TrackingData.Gains = handles.fcor / handles.FrameRate; %[handles.apogain handles.posgain handles.fasgain];
            TrackingData.Parallel = handles.do_parfor.Value;
            TrackingData.info = "Processing [TimTrack; Opticflow], %%\nBlockSize [width; height], %%\nGains [Apo, Position, Angle]";

            if isfield(handles.Region(i),'fas_length_corr') && ~isempty(handles.Region(i).fas_length_corr)
                Fdat.Region(i).FL = handles.Region(i).fas_length_corr(nz,:)';
                Fdat.Region(i).PEN = handles.Region(i).fas_pen_corr(nz,:)';
                Fdat.Region(i).Time = T';

                for j = 1:size(handles.Region(i).fas_length_corr,2)
                    if sum(handles.Region(i).fas_length_corr(nz,j)) == 0

                        Fdat.Region(i).FL(j,:) = handles.Region(i).fas_length(nz,j)';
                        Fdat.Region(i).PEN(j,:) = handles.Region(i).fas_pen(nz,j)';

                    end
                end

            else
                Fdat.Region(i).FL = handles.Region(i).fas_length(nz,:)';
                Fdat.Region(i).PEN = handles.Region(i).fas_pen(nz,:)';
                Fdat.Region(i).Time = time(nz);
            end

        end
    end
end

if isfield(handles, 'D') && get(handles.export_check,'Value')

    %ExpAtRate = inputdlg('Do you want to maintain the original data sampling rate? (y/n)');

    ExpAtRate = 'y';
    % determine the variables highlighted in list
    %cur_var = get(handles.variable_list,'Value');
    var_list = get(handles.variable_list,'String');

    if ~iscell(var_list)
        var_list = {var_list};
    end

    Data =[];
    for i = 1:length(var_list)

        if strcmp(ExpAtRate,'y')

            varname{i} = handles.D(i).name;
            Data.(regexprep(varname{i},' ','_')) = handles.D(i).value(:,:);
            if strcmp(handles.D(i).type,'Digital') && ~exist('Tdigital','var')
                Tdigital  = handles.D(cur_var(i)).time;
            elseif strcmp(handles.D(i).type,'Analog') && ~exist('Tanalog','var')
                Tanalog  = handles.D(i).time;
            end
            if strcmp(handles.D(i).type,'Biodex_Analog') && ~exist('Tanalog','var')
                Tanalog  = handles.D(i).time;
            end

        else
            if strcmp(handles.D(i).type,'Digital')
                new_col = zeros(size(T));
                for j = 1:size(handles.D(i).time,2)
                    new_col(find(abs(T-handles.D(i).time(1,j)) == ...
                        min(abs(T-handles.D(i).time(1,j)))),1) = ...
                        handles.D(i).value(2,j);
                end
            else
                new_col = interp1(handles.D(i).time, ...
                    handles.D(i).value, T, 'linear', 'extrap');

            end

            varname{i} = handles.D(i).name;
            Data.(regexprep(varname{i},' ','_')) = new_col;
        end
    end
    if strcmp(ExpAtRate,'y')
        if exist('Tanalog','var')
            Data.AnalogTime = Tanalog;
        end
        if exist('Tdigital','var')
            Data.DigitalTime = Tdigital;
        end
    else
        Data.Time = T;
    end

end



[~,handles.file,handles.ext] = fileparts(handles.fname);
fileout_suggest = [handles.pname handles.file '.mat'];
[fileout, pathout] = uiputfile(fileout_suggest,'Save fascicle data as...');
cd(pathout);

if exist([pathout fileout], 'file') == 2

    if exist('Data','var') && exist('Fdat','var')
        save(fileout,'-append','Fdat','Data','TrackingData')
    elseif exist('Fdat','var')
        save(fileout,'-append','Fdat','TrackingData')
    elseif exist('Data','var')
        save(fileout,'-append','Data')
    else
        warndlg('There was no data to save', 'WARNING!')
    end

else

    if exist('Data','var') && exist('Fdat','var')
        save(fileout,'Fdat','Data','TrackingData')
    elseif exist('Fdat','var')
        save(fileout,'Fdat','TrackingData')
    elseif exist('Data','var')
        save(fileout,'Data')
    else
        warndlg('There was no data to save', 'WARNING!')
    end

end


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

        show_image(hObject,handles)
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

        show_image(hObject,handles)
    end
end


% --- Executes on button press in data_informed_checkbox.
function data_informed_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to data_informed_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of data_informed_checkbox



function MinDiff_Callback(hObject, eventdata, handles)
% hObject    handle to MinDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinDiff as text
%        str2double(get(hObject,'String')) returns contents of MinDiff as a double


% --- Executes during object creation, after setting all properties.
function MinDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlayButton.
function PlayButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'Interruptible','on');
show_data(hObject, handles)

frame_no = round(get(handles.frame_slider,'Value'));
end_frame= get(handles.frame_slider,'Max');
handles.stop = 0;
guidata(hObject,handles);

%set(gcf,'DoubleBuffer','on')

% cla(handles.axes1) 
%not sure this is necessary, if the fascicle is smaller the image (i.e. no
%extrapolation), it is not nice to see

for i = 1:end_frame
    x(i,:) = handles.Region.Fascicle.fas_x{i};
    y(i,:) = handles.Region.Fascicle.fas_y{i};
end
%initiliaze them
lim_x = [0, handles.vidWidth];
lim_y = [0, handles.vidHeight];
%x axis adjustment
if min(x(:)) < 0
    lim_x(1) = min(x(:));
end
if max(x(:)) > handles.vidWidth
    lim_x(2) = max(x(:));
end
%y axis adjustment
if min(y(:)) < 0
    lim_y(1) = min(y(:));
end
if max(y(:)) > handles.vidHeight
    lim_y(2) = max(y(:));
end
%set(handles.axes1,'XLim', [0 -abs(min(x(:))) handles.vidWidth + max(x(:))],'YLim', [min(y(:)) max(y(:))])
set(handles.axes1,'XLim', lim_x,'YLim', lim_y);


for f = frame_no:end_frame
    %    profile on
    stop = get(handles.StopButton,'Value');

    if stop
        set(handles.StopButton,'Value',0.0);
        return
    else
        set(handles.frame_slider,'Value',f);

        %         tic
        % set the string in the frame_number box to the current frame value
        set(handles.frame_number,'String',num2str(f));
        %         toc

        handles = show_image(hObject,handles);
        drawnow
    end
    %    profile viewer

end



% --- Executes on button press in StopButton.
function StopButton_Callback(hObject, eventdata, handles)
% hObject    handle to StopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.StopButton,'Value',1.0)



% --- Executes on button press in Undo_Correction.
function Undo_Correction_Callback(hObject, eventdata, handles)
% hObject    handle to Undo_Correction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i = 1:length(handles.Region)

    if isfield(handles.Region(i).Fascicle,'fas_x_corr')

        handles.Region(i).Fascicle = rmfield(handles.Region(i).Fascicle,'fas_x_corr');
        handles.Region(i).Fascicle = rmfield(handles.Region(i).Fascicle,'fas_y_corr');

    end

    if isfield(handles.Region(i).Fascicle,'fas_x_keyframe')
        handles.Region(i).Fascicle = rmfield(handles.Region(i).Fascicle,'fas_x_keyframe');
        handles.Region(i).Fascicle = rmfield(handles.Region(i).Fascicle,'fas_y_keyframe');
        handles.Region(i).Fascicle = rmfield(handles.Region(i).Fascicle,'current_xy_keyframe');
    end
end

if isfield(handles.Region, 'fas_length_corr')
    handles.Region = rmfield(handles.Region,'fas_length_corr');
    handles.Region = rmfield(handles.Region,'fas_pen_corr');
    if isfield(handles.Region,'fas_pen_keyframe')
        handles.Region = rmfield(handles.Region,'fas_pen_keyframe');
        handles.Region = rmfield(handles.Region,'fas_length_keyframe');
    end

end

handles.KeyframeInd = logical(0.0);
show_image(hObject,handles)



% --------------------------------------------------------------------
function Load_All_Tracked_Frames_Callback(hObject, eventdata, handles)
% hObject    handle to Load_All_Tracked_Frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%load the avi file

menu_clear_tracking_Callback(hObject, eventdata, handles)% clear any current tracking

%select file
[fname, pname] = uigetfile('*.mat', 'Pick a .MAT file');

load([pname fname],'TrackingData');

for i = 1:length(TrackingData.Region)

    % Recompute the ROI using roipoly and the vertices that were saved
    sf = TrackingData.start_frame;% find the start frame and end frame of previous tracking
    ef = (length(TrackingData.Region(i).ROIx)-1)+TrackingData.start_frame;
    s = size(handles.ImStack(:,:,sf:ef));% reorganise the image stack into a temporary cell array
    I = shiftdim(mat2cell(handles.ImStack(:,:,sf:ef),...
        s(1),s(2),ones(1,s(3))),1);
    x = TrackingData.Region(i).ROIx; % get the coordinates of the vertices
    y = TrackingData.Region(i).ROIy;
    TrackingData.Region(i).ROI = cellfun(@roipoly,I,x,y,'UniformOutput',false);%compute ROI

    for j = 1:length(TrackingData.Region(i).Fascicle)

        if isfield(TrackingData.Region(i).Fascicle(j),'fas_x_corr')

            TrackingData.Region(i).Fascicle(j).fas_x = TrackingData.Region(i).Fascicle(j).fas_x_corr;
            TrackingData.Region(i).Fascicle(j).fas_y = TrackingData.Region(i).Fascicle(j).fas_y_corr;
            TrackingData.Region(i).fas_length =  TrackingData.Region(i).fas_length_corr;
            TrackingData.Region(i).fas_pen =  TrackingData.Region(i).fas_pen_corr;

        end

    end
    if isfield(TrackingData.Region(i).Fascicle(j),'fas_x_corr')
        TrackingData.Region(i).Fascicle = rmfield(TrackingData.Region(i).Fascicle,{'fas_x_corr','fas_y_corr'});
    end
end

if isfield(TrackingData.Region(i).Fascicle(j),'fas_length_corr')
    TrackingData.Region = rmfield(TrackingData.Region,{'fas_length_corr','fas_pen_corr'});
end

handles.Region=TrackingData.Region;

handles.start_frame = TrackingData.start_frame;
handles.NumFrames = TrackingData.NumFrames;
set(handles.frame_slider,'Min',1);
set(handles.frame_slider,'Max',handles.NumFrames);
set(handles.frame_slider,'Value',1);
set(handles.frame_slider,'SliderStep',[1/handles.NumFrames 5/handles.NumFrames]);
% set the string in the frame_number to 1
set(handles.frame_number,'String',1);

show_image(hObject,handles);



% --- Executes on button press in FixedROI_on.
function FixedROI_on_Callback(hObject, eventdata, handles)
% hObject    handle to FixedROI_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FixedROI_on


% --- Executes on selection change in EventList.
function EventList_Callback(hObject, eventdata, handles)
% hObject    handle to EventList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EventList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EventList
show_image(hObject,handles);

% --- Executes during object creation, after setting all properties.
function EventList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EventList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadEvents.
function LoadEvents_Callback(hObject, eventdata, handles)
% hObject    handle to LoadEvents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.ext,'.c3d')

    acq = btkReadAcquisition(handles.dat_file);
    events = btkGetEvents(acq);
    eventnames = fieldnames(events);
    set(handles.EventList,'String',eventnames)
    handles.c3d_events = events;

    show_image(hObject,handles);
else

    warndlg('Can only load events from c3d format');
    return

end

% --- Executes on button press in DeleteEvent.
function DeleteEvent_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Event2keyframe.
function Event2keyframe_Callback(hObject, eventdata, handles)
% hObject    handle to Event2keyframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

elist = get(handles.EventList,'String');

if isempty(elist)
    warndlg('no events loaded')
    return
else
    i = get(handles.EventList,'Value');
    ename = elist{i};
    acq = btkReadAcquisition(handles.dat_file);
    events = btkGetEvents(acq);
    eventT = events.(ename);
    for k = 1:length(eventT)
        [~,eventF(k)] = min(abs(handles.Time-eventT(k)));
    end
    eventF(eventF<handles.start_frame)=[];
    eventF(eventF>(handles.start_frame+get(handles.frame_slider,'Max')))=[];
    eventF=unique((eventF+1)-handles.start_frame);
    for k = 1:length(eventF)
        eventFc{k} = num2str(eventF(k));
    end
    set(handles.keyframe_list,'String',eventFc);
    set(handles.keyframe_list,'Value',1);
end


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

for i = 1:length(handles.Region)
    if isfield(handles.Region(i),'Fascicle')
        if isfield(handles.Region(i),'fas_length')

            if ~isempty(handles.Region(i).fas_length)

                dt = 1/handles.FrameRate;
                nz = logical(handles.Region(i).fas_length(:,1) ~= 0);
                FL = handles.Region(i).fas_length(nz,:);
                PEN = handles.Region(i).fas_pen(nz,:) * 180/pi;

                time = (dt:dt:length(handles.Region(i).fas_length(:,1))*dt)+((handles.start_frame-1)*dt);

                for j = 1:length(handles.Region(i).Fascicle)

                    axes(handles.length_plot); hold off;
                    plot(handles.length_plot,time(nz),FL,'r','linewidth',2);
                    set(handles.length_plot,'ylim',[min(FL)*0.85 max(FL)*1.15]); %set axis 15% difference of min and and value,easier to read
                    box off;

                    axes(handles.mat_plot); hold off;
                    plot(handles.mat_plot,time(nz),PEN,'r','linewidth',2);
                    set(handles.mat_plot,'ylim',[min(PEN)*0.85 max(PEN)*1.15]); %set axis 15% difference of min and and value,easier to read
                    box off;
                end
            end
        end
    end
end

set(handles.length_plot,'XLim',[0 handles.Time(end)]);
%     L_range = get(handles.length_plot,'YLim');

xlabel(handles.length_plot,'Time (s)')
ylabel(handles.length_plot,'Fascicle Length (mm)')

set(handles.mat_plot,'XLim',[0 handles.Time(end)]);
xlabel(handles.mat_plot,'Time (s)')
ylabel(handles.mat_plot,'Pennation (deg)')

%---------------------------------------------------------
% Function to show image with appropriate image processing
%---------------------------------------------------------
function handles = show_image(hObject, handles)

if isfield(handles,'ImStack')
    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'));

    % show the image
    if isfield(handles,'ImStack')
        Im = handles.ImStack(:,:,frame_no+handles.start_frame-1);

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

    % Update Im and NIm
    handles.Im = Im;
    handles.NIm = handles.Im;

    % Show the fascicle and ROI
    if isfield(handles,'Region')
        for i = 1:length(handles.Region)
            for j = 1:length(handles.Region(i).Fascicle)

                if ~isfield(handles, 'fascicle') || ~isvalid(handles.fascicle)
                    handles.fascicle = line(handles.axes1,'xdata',handles.Region(i).Fascicle(j).fas_x{frame_no}, 'ydata', handles.Region(i).Fascicle(j).fas_y{frame_no}, ...
                        'color','Red', 'linewidth',2,'marker','o','markerfacecolor','Red');
                else
                    set(handles.fascicle, 'xdata',handles.Region(i).Fascicle(j).fas_x{frame_no}, 'ydata', handles.Region(i).Fascicle(j).fas_y{frame_no})
                end

                if ~isfield(handles, 'ROI') || ~isvalid(handles.ROI)
                    handles.ROI = line(handles.axes1,'xdata',handles.Region(i).ROIx{frame_no}, 'ydata', handles.Region(i).ROIy{frame_no}, ...
                        'linestyle',':','color','Red', 'linewidth',2);
                else
                    set(handles.ROI, 'xdata',handles.Region(i).ROIx{frame_no}, 'ydata', handles.Region(i).ROIy{frame_no})
                end

            end
        end
    end
    
    
    % remove previous vertical lines
    children = get(handles.length_plot, 'children');
    if length(children) > 1
        delete(children(1));
    end

    children = get(handles.mat_plot, 'children');
    if length(children) > 1
        delete(children(1));
    end

    FL = handles.Region(1).fas_length;
    PEN = handles.Region(1).fas_pen * 180/pi;
    
    % add new vertical lines
    line(handles.length_plot, 'xdata', handles.Time(frame_no) * ones(1,2), 'ydata', [.85*min(FL) 1.15*max(FL)],'color',[0 0 0]);
    line(handles.mat_plot, 'xdata', handles.Time(frame_no) * ones(1,2), 'ydata', [.85*min(PEN) 1.15*max(PEN)],'color', [0 0 0]);
       
    handles.prev_frame_no = frame_no;

    % Update handles structure
    guidata(hObject, handles);

%     profile viewer
    
%     keyboard


end

% --- Executes on button press in process_all.
function process_all_Callback(hObject, eventdata, handles)
% hObject    handle to process_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Run TimTrack
handles = process_all_TimTrack(hObject, eventdata, handles);

% Run UltraTrack (note: includes state estimation on ROI)
handles = process_all_UltraTrack(hObject, eventdata, handles);
   
% Update handles structure
guidata(hObject, handles);

% update the image axes using show_image function (bottom)
show_image(hObject,handles)

show_data(hObject, handles)

function[handles] = process_all_UltraTrack(hObject, eventdata, handles)

    %% Optical flow and state estimation
    % setup current and new image
    frame_no = 1;
    
    im1 = imcrop(handles.ImStack(:,:,1),handles.crop_rect);
    h = waitbar(0,['Processing frame 1/', num2str(handles.NumFrames)],'Name','Running UltraTrack...');

    for i = 1:length(handles.Region)
            
        points = detectMinEigenFeatures(im1,'FilterSize',11, 'MinQuality', 0.005);
        points = double(points.Location);
        [inPoints] = inpolygon(points(:,1),points(:,2), handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no});
        points = points(inPoints,:);

        %% Define a smaller area with features pts around the fascicle line
        %get points around fascicle line
        %ptsFas = [ handles.Region.Fascicle.fas_x{1,1}, handles.Region.Fascicle.fas_y{1,1}] ;

        % Calculate the distance between each point of interest and the line defined by ptsFas
        %distances = pointToLineDistance(points, ptsFas);

        % Set the maximum distance within which a point is considered around the line
        %maxDistance = round(2 / handles.ID); % 2mm, adjust according to 10% of pixels to mm

        % Keep only the points that are within the maximum distance from the line
        % selectedPoints = points(distances <= maxDistance, :);
        %
        % points = selectedPoints;

        %% Initialize point tracker and run opticflow
        tstart = tic;
        %calculate block size according to ROI 
        width = floor(max(abs(diff(handles.ROI.XData))) * 0.20);
        height = floor(max(abs(diff(handles.ROI.YData))) * 0.40); %thickness changes?

        handles.BlockSize = [width height]; %save as width and height for later comparison

        % Ensure width and height are both odd numbers
        if mod(width, 2) == 0
            width = width + 1; % Increment by 1 to make it odd
        end

        if mod(height, 2) == 0
            height = height + 1; % Increment by 1 to make it odd
        end

        pointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',[height width]);
        %             pointTracker = vision.PointTracker('NumPyramidLevels',1,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',[21 71]);
        initialize(pointTracker,points,im1);


        for f = frame_no+1:get(handles.frame_slider,'Max')
            % Get the current image

%           profile on
            im2 = imcrop(handles.ImStack(:,:,handles.start_frame+f-1),handles.crop_rect);
            handles.NIm = im2;

            % Compute the flow and new roi
            [pointsNew, isFound] = step(pointTracker, im2);
            [w,~] = estimateGeometricTransform2D(points(isFound,:), pointsNew(isFound,:), 'affine', 'MaxDistance',50);
            handles.Region(i).warp(:,:,f) = w;

            handles.CIm = handles.NIm;

            for j = 1:length(handles.Region(i).Fascicle)

                % optical flow
                handles = apply_transform(handles,f,f-1,i,j);

                % state estimation
                handles = ROI_state_estimator(handles,f,i,j);

                % option 1: detect new points in each frame
%                 points = detectMinEigenFeatures(im2,'FilterSize',11, 'MinQuality', 0.005);
%                 points = double(points.Location);
                
                % option 2: use point from previous
                points = pointsNew;
                
                % set the points
                inPoints = inpolygon(points(:,1),points(:,2), handles.Region(i).ROIx{f},handles.Region(i).ROIy{f});
                points = points(inPoints,:);
                setPoints(pointTracker, points);
                
                % calculate the length and pennation for the current frame
                handles.Region(i).fas_pen(f,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y{f})),...
                abs(diff(handles.Region(i).Fascicle(j).fas_x{f})));

                scalar = handles.ID;%/handles.vidHeight;

                handles.Region(i).fas_length(f,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{f}).^2 +...
                diff(handles.Region(i).Fascicle(j).fas_x{f}).^2);
        
            end

%             profile viewer
            frac_progress = (f+(get(handles.frame_slider,'Max')*(i-1))) / (get(handles.frame_slider,'Max')*length(handles.Region));
            waitbar(frac_progress,h, ['Processing frame ', num2str(f), '/', num2str(get(handles.frame_slider,'Max'))])
        
        end


        % correct


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
    [geofeatures, ~, parms] = auto_ultrasound(im2(:,:,1), parms);
    
    % call again a bunch of times to get estimate the total duration
    for i = 1:min([size(im2,3), 5])
        tstart = tic;
        geofeatures = auto_ultrasound(im2(:,:,1), parms);
        dt = toc(tstart);
    end
    
    est_duration = dt * numIterations;
    
	% prompt to optionally change processing based on computational time
    answer = 'Undefined';
    if (est_duration < 60) && handles.do_parfor.Value == 1
        answer = questdlg(['Estimated TimTrack duration: ', num2str(est_duration), ' s, would you like to use parallel pool?'], 'Type of computation', 'Yes','No','Cancel','No');
    elseif (est_duration > 60) && handles.do_parfor.Value == 0
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

       % Adjust the parameter of geofeat
        for kk = 1:length(geofeatures)
            geofeatures(kk).super_coef(2) = geofeatures(kk).super_coef(2) * handles.imresize_fac;
            geofeatures(kk).deep_coef(2) = geofeatures(kk).deep_coef(2) * handles.imresize_fac;
            geofeatures(kk).faslen = geofeatures(kk).faslen * handles.imresize_fac;
        end

        handles.geofeatures = geofeatures;

    end
end


function[handles] = ROI_state_estimator(handles,frame_no,i,j)

geofeatures = handles.geofeatures;

n = handles.vidWidth;

% Gains
g = handles.fcor / handles.FrameRate;
g(g>1) = 1;
g(g<0) = 0;

% Hough estimate of vertical aponeurosis position
ROIy_est = round([polyval(geofeatures(frame_no).super_coef, 1) polyval(geofeatures(frame_no).deep_coef, [1 n]) polyval(geofeatures(frame_no).super_coef, [n 1])])';

% State estimation on ROI
ROIy_cor = handles.Region(i).ROIy{frame_no} + g(1) * (ROIy_est - handles.Region(i).ROIy{frame_no});
handles.Region(i).ROIy{frame_no} = ROIy_cor;

function[handles] = state_estimator(handles,frame_no,i,j, direction)

geofeatures = handles.geofeatures;

% Gains
g = handles.fcor / handles.FrameRate;
g(g>1) = 1;
g(g<0) = 0;

%% Get the current fascicle
% Fascicle - Aponeurosis intersection points from optical flow
x = handles.Region(i).Fascicle(j).fas_x{frame_no};
y = handles.Region(i).Fascicle(j).fas_y{frame_no};

% Length and angle from optical flow
alpha = atan2d(abs(diff(y)), abs(diff(x)));
L = sqrt(diff(x).^2 + diff(y).^2);

%% Superficial intersection drift estimate
% With respect to the initial position
x20 = handles.Region(i).Fascicle(j).fas_x{1}(2);
dx2_drift = x(2) - x20;

% For vertical position, assume Hough is correct
y20 = polyval(geofeatures(frame_no).super_coef, x20);
dy2_drift = y(2) - y20;

%% Fascicle angle drift estimate
dalpha_drift = alpha - geofeatures(frame_no).alpha;

%% Deep intersection drift estimate
% option 1: take length from TimTrack
% L_TT = geofeatures(frame_no).faslen;

% option 2: use the deep aponeurosis from TimTrack
fas_coef(1) = -tand(alpha);
fas_coef(2) =  y(2) - fas_coef(1) * x(2);
x1_TT = (fas_coef(2) - geofeatures(frame_no).deep_coef(2)) / (geofeatures(frame_no).deep_coef(1) - fas_coef(1));
y1_TT = polyval(geofeatures(frame_no).deep_coef, x1_TT);
L_TT = sqrt(diff([x1_TT x(2)]).^2 + diff([y1_TT y(2)]).^2); % length before Hough

% correct
dL_drift = L - L_TT;

%% Correct the drift
% correction opposing the drift
dx2_cor = -g(2) * dx2_drift;
dy2_cor = -g(2) * dy2_drift;

% length and angle
alpha_cor   = -g(3) * dalpha_drift;
L_cor       = -g(2) * dL_drift;

% determine new lengh and angle to calculate deep correction
alpha_new = alpha + alpha_cor; 
L_new = L + L_cor;
x2_new = x(2) + dx2_cor;
y2_new = y(2) + dy2_cor;

x1_new = x2_new - cosd(alpha_new) * L_new;
y1_new = y2_new + sind(alpha_new) * L_new;

% infer the correction
dx1_cor = x1_new - x(1);
dy1_cor = y1_new - y(1);

if strcmp(direction,'forward') % apply to all frames in the future
    for k = frame_no:handles.NumFrames
        handles.Region(i).Fascicle(j).fas_x{k} = handles.Region(i).Fascicle(j).fas_x{k} + [dx1_cor dx2_cor];
        handles.Region(i).Fascicle(j).fas_y{k} = handles.Region(i).Fascicle(j).fas_y{k} + [dy1_cor dy2_cor];
    end
    
elseif strcmp(direction,'backward') % apply to all frames in the past
    for k = frame_no:-1:1
        handles.Region(i).Fascicle(j).fas_x{k} = handles.Region(i).Fascicle(j).fas_x{k} + [dx1_cor dx2_cor];
        handles.Region(i).Fascicle(j).fas_y{k} = handles.Region(i).Fascicle(j).fas_y{k} + [dy1_cor dy2_cor];
    end
end

% calculate the length and pennation for the current frame
handles.Region(i).fas_pen(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y{frame_no})),...
    abs(diff(handles.Region(i).Fascicle(j).fas_x{frame_no})));

scalar = handles.ID;%/handles.vidHeight;

handles.Region(i).fas_length(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{frame_no}).^2 +...
    diff(handles.Region(i).Fascicle(j).fas_x{frame_no}).^2);


% function to perform transformation (warp) of the points
function handles = apply_transform(handles,frame_no,prev_frame_no,i,j)

% extract the warp matrix
w = handles.Region(i).warp(:,:,frame_no);

% ROI
ROIpos = transformPointsForward(w, [handles.Region(i).ROIx{prev_frame_no} handles.Region(i).ROIy{prev_frame_no}]);
ROIpos(:,1) = handles.Region(i).ROIx{prev_frame_no};

[handles.Region(i).ROI{frame_no},handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no}] = roipoly(handles.NIm, ROIpos(:,1), ROIpos(:,2));

% Adjust the ROI if it goes outside the image
handles.Region(i).ROIy{frame_no}(handles.Region(i).ROIy{frame_no} > handles.vidHeight) = handles.vidHeight;
handles.Region(i).ROIy{frame_no}(handles.Region(i).ROIy{frame_no} < 1) = 1;

%% Fascicle
% loop through all fasicles defined for the region and apply
% the warp and calculate new fascicle length

% book keeping
handles.Region(i).Fascicle(j).current_xy(1,1) = handles.Region(i).Fascicle(j).fas_x{prev_frame_no}(1);
handles.Region(i).Fascicle(j).current_xy(1,2) = handles.Region(i).Fascicle(j).fas_y{prev_frame_no}(1);
handles.Region(i).Fascicle(j).current_xy(2,1) = handles.Region(i).Fascicle(j).fas_x{prev_frame_no}(2);
handles.Region(i).Fascicle(j).current_xy(2,2) = handles.Region(i).Fascicle(j).fas_y{prev_frame_no}(2);

% actual transform
handles.Region(i).Fascicle(j).new_xy = transformPointsForward(w, handles.Region(i).Fascicle(j).current_xy);

% update
handles.Region(i).Fascicle(j).current_xy = handles.Region(i).Fascicle(j).new_xy;

% define Hough variables
handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no} = [];

% more book keeping
handles.Region(i).Fascicle(j).fas_x{frame_no}(1) = handles.Region(i).Fascicle(j).current_xy(1,1);
handles.Region(i).Fascicle(j).fas_y{frame_no}(1) = handles.Region(i).Fascicle(j).current_xy(1,2);
handles.Region(i).Fascicle(j).fas_x{frame_no}(2) = handles.Region(i).Fascicle(j).current_xy(2,1);
handles.Region(i).Fascicle(j).fas_y{frame_no}(2) = handles.Region(i).Fascicle(j).current_xy(2,2);

% save the original (never changes)
handles.Region(i).Fascicle(j).fas_x_original{frame_no} = handles.Region(i).Fascicle(j).fas_x{frame_no};
handles.Region(i).Fascicle(j).fas_y_original{frame_no} = handles.Region(i).Fascicle(j).fas_y{frame_no};


% --- Executes on button press in Auto_Detect.
function Auto_Detect_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_Detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% add 3 to path and load TimTrack parameters
%     ccd = cd;
%     cd(handles.TimTrackfolder.String)
%     addpath(genpath(cd));
%     cd('Parameters')
load('parms.mat','parms')
%     cd(ccd)

% find current frame number from slider
frame_no = round(get(handles.frame_slider,'Value'));

% show the image
handles.IM = handles.ImStack(:,:,frame_no+handles.start_frame-1);

img = im2double(handles.IM);
axes(handles.axes1); hold off

im2 = imresize(img, 1/handles.imresize_fac);

%% Aponeurosis detection
for jj = 1:2
    [~,apCentre] = ginputYellow(1);

    apCentre = round(apCentre/handles.vidHeight,2)*100;
    apRound(jj) = round(apCentre,-1);
end

%data = img(handles.crop_rect(2):(handles.crop_rect(2)+handles.crop_rect(4)-1), handles.crop_rect(1):(handles.crop_rect(1)+handles.crop_rect(3)-1));
data = im2;

% don't use TimTrack's figure display, because we already have this GUI
parms.show = 0;
parms.fas.show = 0;

% need to be more lenient for broad range of muscles
parms.apo.deep.maxangle = 10;

% some default parameters
range = 15;

% make range dependent on user-picked locations
parms.apo.super.cut = [max(apRound(1)-range, 0), apRound(1)+range] / 100;
parms.apo.deep.cut = [apRound(2)-range, min(apRound(2)+range, 100)] / 100;

% run TimTrack
[geofeatures, ~, parms] = auto_ultrasound(data, parms);

% scale back
geofeatures.super_coef(2) = geofeatures.super_coef(2) * handles.imresize_fac;
geofeatures.deep_coef(2) = geofeatures.deep_coef(2) * handles.imresize_fac;
geofeatures.fas_coef(2) = geofeatures.fas_coef(2) * handles.imresize_fac;
geofeatures.faslen = geofeatures.faslen * handles.imresize_fac;

handles.parms = parms;

%frame_no = 1;

n = handles.vidWidth;
i = 1; j = 1;

Deep_intersect_x = round((geofeatures.deep_coef(2) - geofeatures.fas_coef(2))   ./ (geofeatures.fas_coef(1) - geofeatures.deep_coef(1)));
Super_intersect_x = round((geofeatures.super_coef(2) - geofeatures.fas_coef(2)) ./ (geofeatures.fas_coef(1) - geofeatures.super_coef(1)));
Super_intersect_y = polyval(geofeatures.super_coef, Super_intersect_x);
Deep_intersect_y = polyval(geofeatures.deep_coef, Deep_intersect_x);

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
        scalar = handles.ID;
        handles.CIm = handles.Im;
    end

    handles.Region(i).fas_length(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{frame_no}).^2 + diff(handles.Region(i).Fascicle(j).fas_x{frame_no}).^2);

    if ~isfield(handles.Region(i).Fascicle(j),'analysed_frames')
        handles.Region(i).Fascicle(j).analysed_frames = frame_no;
    else handles.Region(i).Fascicle(j).analysed_frames = sort([handles.Region(i).Fascicle(j).analysed_frames frame_no]);
    end

end


% Update handles structure
guidata(hObject, handles);

show_image(hObject,handles)


function gain_Callback(hObject, eventdata, handles)
% hObject    handle to gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gain as text
%        str2double(get(hObject,'String')) returns contents of gain as a double

handles.fcor(3) = str2double(get(hObject,'String')); % [Hz]

% Update handles structure
guidata(hObject, handles);

% If we already processed, run state estimation
if length(handles.Region.Fascicle.fas_x) == handles.NumFrames
do_state_estimation(hObject, eventdata, handles)
end



% --- Executes during object creation, after setting all properties.
function gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gain (see GCBO)
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


function do_state_estimation(hObject, eventdata, handles)
% hObject    handle to do_state_estimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% start with the original
for f = 1:get(handles.frame_slider,'Max')
    for i = 1:length(handles.Region)
        for j = 1:length(handles.Region(i).Fascicle)
            handles.Region(i).Fascicle(j).fas_x{f} = handles.Region(i).Fascicle(j).fas_x_original{f};
            handles.Region(i).Fascicle(j).fas_y{f} = handles.Region(i).Fascicle(j).fas_y_original{f};
        end
    end
end

% forward state estimation
for f = 1:get(handles.frame_slider,'Max')
    for i = 1:length(handles.Region)
        for j = 1:length(handles.Region(i).Fascicle)
            % state estimation
            handles = state_estimator(handles,f,i,j,'forward');
        end
    end
end

% backward state estimation
for f = get(handles.frame_slider,'Max'):-1:1
    for i = 1:length(handles.Region)
        for j = 1:length(handles.Region(i).Fascicle)
            % state estimation
            handles = state_estimator(handles,f,i,j,'backward');
        end
    end
end

show_data(hObject, handles)
guidata(hObject, handles);

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
handles.ImStack = flip(handles.ImStack, 2);
%end
guidata(hObject, handles);
show_image(hObject,handles)



function apo_gain_Callback(hObject, eventdata, handles)
% hObject    handle to apo_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of apo_gain as text
%        str2double(get(hObject,'String')) returns contents of apo_gain as a double

handles.fcor(1) = str2double(get(hObject,'String'));

% Run optical flow
handles = process_all_UltraTrack(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% Show image
show_image(hObject,handles)


% --- Executes during object creation, after setting all properties.
function apo_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apo_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.fcor(1) = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);



function pos_gain_Callback(hObject, eventdata, handles)
% hObject    handle to pos_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pos_gain as text
%        str2double(get(hObject,'String')) returns contents of pos_gain as a double

handles.fcor(2) = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% If we have estimates, run state estimation
if length(handles.Region.Fascicle.fas_x) == handles.NumFrames 
    do_state_estimation(hObject, eventdata, handles)
end


% --- Executes during object creation, after setting all properties.
function pos_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pos_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.fcor(2) = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);


% --- Function to calculate euclidean distance form the fascicle line and
% the featuers points to keep around that fascicle (NOT USED NOW)
function distances = pointToLineDistance(points, line)
    % Calculate distance between each point and the line
    x1 = line(1, 1);
    y1 = line(1, 3);
    x2 = line(1, 2);
    y2 = line(1, 4);

    distances = abs((y2 - y1) * points(:, 1) - (x2 - x1) * points(:, 2) + x2 * y1 - y2 * x1) / sqrt((y2 - y1)^2 + (x2 - x1)^2);

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
disp('Parallel pool shutted down!');
% Hint: delete(hObject) closes the figure
delete(hObject);



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

