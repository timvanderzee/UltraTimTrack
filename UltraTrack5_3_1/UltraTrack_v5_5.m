function varargout = UltraTrack_v5_5(varargin)
% ULTRATRACK_V5_4 M-file for UltraTrack_v5_5.fig
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
%      applied to the GUI before UltraTrack_v5_5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UltraTrack_v5_5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UltraTrack_v5_5

% Last Modified by GUIDE v2.5 01-Feb-2024 10:02:31
% Last Modified by Paolo Tecchio 17/08/2022
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UltraTrack_v5_5_OpeningFcn, ...
    'gui_OutputFcn',  @UltraTrack_v5_5_OutputFcn, ...
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


% --- Executes just before UltraTrack_v5_5 is made visible.
function UltraTrack_v5_5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UltraTrack_v5_5 (see VARARGIN)

% Choose default command line output for UltraTrack_v5_5
handles.output = hObject;

%add automatically all files and subfolders dynamically
tmp = mfilename;
tmp = erase(which(tmp),[tmp '.m']);
addpath(genpath(tmp));

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

% UIWAIT makes UltraTrack_v5_5 wait for user response (see UIRESUME)
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
function varargout = UltraTrack_v5_5_OutputFcn(~, eventdata, handles)
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

    if handles.Region(i).Fascicle(j).fas_y{frame_no}(1) > handles.Region(i).Fascicle(j).fas_y{frame_no}(2)
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
    else handles.Region(i).Fascicle(j).analysed_frames = sort([handles.Region(i).Fascicle(j).analysed_frames frame_no]);
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

    axes(handles.length_plot)
    cla

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

    axes(handles.length_plot)
    cla

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

% added by Tim
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
    axes(handles.axes1)

    % use imcrop tool to determine croppable area
    [~,handles.crop_rect] = imcrop;

    handles.crop_rect = round(handles.crop_rect);
    handles.vidHeight = handles.crop_rect(4);
    handles.vidWidth = handles.crop_rect(3);

    % update the image axes using show_image function (bottom)
    show_image(hObject,handles)

end

% --------------------------------------------------------------------
function menu_reset_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_reset_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack')

    % set the image croppable area to the maximum area
    handles.crop_rect = [1 1 handles.vidWidth handles.vidHeight];

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
            
%             waitbar(i/get(handles.frame_slider,'Max'))
            
            Im = handles.ImStack(:,:,i);
            % make the image axes the current axes
            
            %cla(handles.axes1)
            % crop image with current cropping value
            Im = imcrop(Im,handles.crop_rect); 
            
            % show the image
            hold off; 
            imshow(Im,'Border','loose'); hold on
            colormap(gray(256)); axis off;
            
             plot(handles.axes1, handles.Region.Fascicle.fas_x{i},handles.Region.Fascicle.fas_y{i},'ro-','LineWidth',2)
             plot(handles.axes1, [xrange(1) xrange(2)], polyval(handles.super_coef{i}, [xrange(1) xrange(2)]),'r--');
             plot(handles.axes1, [xrange(1) xrange(2)], polyval(handles.deep_coef{i}, [xrange(1) xrange(2)]),'r--');
             
             set(handles.axes1, 'Xlim', [xrange(1) xrange(2)]);
             
            
%             set(handles.frame_slider,'Value',i)
            
%             profile on
%             show_image(hObject, handles)
%             profile viewer
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

frame_no = round(get(handles.frame_slider,'Value'));
end_frame= get(handles.frame_slider,'Max');
handles.stop = 0;
guidata(hObject,handles);

%set(gcf,'DoubleBuffer','on')

for f = frame_no:end_frame;
    stop = get(handles.StopButton,'Value');
    if stop
        set(handles.StopButton,'Value',0.0);
        return
    else
        set(handles.frame_slider,'Value',f);
        % set the string in the frame_number box to the current frame value
        set(handles.frame_number,'String',num2str(f));
        show_image(hObject,handles)
        drawnow
    end

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


function ImDepthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ImDepthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImDepthEdit as text
%        str2double(get(hObject,'String')) returns contents of ImDepthEdit as a double
handles.ID = str2double(get(handles.ImDepthEdit,'String'));
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


% --- Executes on button press in U_fix.
function U_fix_Callback(hObject, eventdata, handles)
% hObject    handle to U_fix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k = get(handles.no_tracked_regions,'Value');
j = get(handles.no_tracked_fascicles,'Value');

if isfield(handles.Region(k).Fascicle(j),'fas_y')

    if isfield(handles.Region(k).Fascicle(j),'fas_y_corr') && ~isempty(handles.Region(k).Fascicle(j).fas_y_corr)
        handles.Region(k).Fascicle(j).fas_x = handles.Region(k).Fascicle(j).fas_x_corr;
        handles.Region(k).Fascicle(j).fas_y = handles.Region(k).Fascicle(j).fas_y_corr;
    else
        handles.Region(k).Fascicle(j).fas_x_corr = handles.Region(k).Fascicle(j).fas_x;
        handles.Region(k).Fascicle(j).fas_y_corr = handles.Region(k).Fascicle(j).fas_y;
        handles.Region(k).fas_length_corr(:,j) = handles.Region(k).fas_length(:,j);
        handles.Region(k).fas_pen_corr(:,j) = handles.Region(k).fas_pen(:,j);
    end

    % find current frame number from slilder
    frame_no = round(get(handles.frame_slider,'Value'));

    gcf;
    axes(handles.axes1);
    ax=gca;

    p = [handles.Region(k).Fascicle(j).fas_x{frame_no}(2),...
        handles.Region(k).Fascicle(j).fas_y{frame_no}(2)];

    Cx = handles.Region(k).Fascicle(j).fas_x{frame_no};
    Cy = handles.Region(k).Fascicle(j).fas_y{frame_no};

    [p1,p2]=rbline2(ax,p);

    handles.Region(k).Fascicle(j).fas_x_corr{frame_no}(1)=p2(1);
    handles.Region(k).Fascicle(j).fas_y_corr{frame_no}(1)=p2(2);
    handles.Region(k).fas_pen_corr(frame_no,j) = atan2(abs(diff(handles.Region(k).Fascicle(j).fas_y_corr{frame_no})),abs(diff(handles.Region(k).Fascicle(j).fas_x_corr{frame_no})));
    scalar = handles.ID;%/handles.vidHeight;
    handles.Region(k).fas_length_corr(frame_no,j) = scalar*sqrt(diff(handles.Region(k).Fascicle(j).fas_y_corr{frame_no}).^2 + diff(handles.Region(k).Fascicle(j).fas_x_corr{frame_no}).^2);

    % check if the ROI needs adjusting too
    if ~get(handles.FixedROI_on,'Value')
        answr = questdlg('Redefine Region of Interest?');
        switch answr
            case 'Yes'
                axes(handles.axes1)
                [handles.Region(k).ROI{frame_no},handles.Region(k).ROIx{frame_no},...
                    handles.Region(k).ROIy{frame_no}] = roipoly;

                show_image(hObject,handles)

            case 'No'

            case 'Cancel'
        end
    else
        answr = 'No';
    end

    wb = waitbar(0,'Recalculating ...');

    EF = get(handles.frame_slider,'Max');


    for f = frame_no+1:EF

        if ~isempty(find(handles.Region(k).Fascicle(j).analysed_frames == f,1))

            handles = apply_transform(handles,f,f-1,k,j);

        else

            waitbar(1,wb);
            break

        end
        waitbar((f-frame_no)/(EF-frame_no),wb);
    end


    close(wb)
    show_image(hObject,handles)

else
    warndlg('Selected Fascicle Data Not Found','WARNING')
end


% --- Executes on button press in L_fix.
function L_fix_Callback(hObject, eventdata, handles)
% hObject    handle to L_fix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k = get(handles.no_tracked_regions,'Value');
j = get(handles.no_tracked_fascicles,'Value');

if isfield(handles.Region(k).Fascicle(j),'fas_y')

    if isfield(handles.Region(k).Fascicle(j),'fas_y_corr') && ~isempty(handles.Region(k).Fascicle(j).fas_y_corr)
        handles.Region(k).Fascicle(j).fas_x = handles.Region(k).Fascicle(j).fas_x_corr;
        handles.Region(k).Fascicle(j).fas_y = handles.Region(k).Fascicle(j).fas_y_corr;
    else
        handles.Region(k).Fascicle(j).fas_x_corr = handles.Region(k).Fascicle(j).fas_x;
        handles.Region(k).Fascicle(j).fas_y_corr = handles.Region(k).Fascicle(j).fas_y;
        handles.Region(k).fas_length_corr(:,j) = handles.Region(k).fas_length(:,j);
        handles.Region(k).fas_pen_corr(:,j) = handles.Region(k).fas_pen(:,j);
    end

    % find current frame number from slilder
    frame_no = round(get(handles.frame_slider,'Value'));

    gcf;
    axes(handles.axes1);
    ax=gca;

    p = [handles.Region(k).Fascicle(j).fas_x{frame_no}(1),...
        handles.Region(k).Fascicle(j).fas_y{frame_no}(1)];

    Cx = handles.Region(k).Fascicle(j).fas_x{frame_no};
    Cy = handles.Region(k).Fascicle(j).fas_y{frame_no};

    [p1,p2]=rbline2(ax,p);

    handles.Region(k).Fascicle(j).fas_x_corr{frame_no}(2)=p2(1);
    handles.Region(k).Fascicle(j).fas_y_corr{frame_no}(2)=p2(2);
    handles.Region(k).fas_pen_corr(frame_no,j) = atan2(abs(diff(handles.Region(k).Fascicle(j).fas_y_corr{frame_no})),abs(diff(handles.Region(k).Fascicle(j).fas_x_corr{frame_no})));
    scalar = handles.ID;%/handles.vidHeight;
    handles.Region(k).fas_length_corr(frame_no,j) = scalar*sqrt(diff(handles.Region(k).Fascicle(j).fas_y_corr{frame_no}).^2 + diff(handles.Region(k).Fascicle(j).fas_x_corr{frame_no}).^2);


    % check if the ROI needs adjusting too
    if ~get(handles.FixedROI_on,'Value')
        answr = questdlg('Redefine Region of Interest?');
        switch answr
            case 'Yes'
                axes(handles.axes1)
                [handles.Region(k).ROI{frame_no},handles.Region(k).ROIx{frame_no},...
                    handles.Region(k).ROIy{frame_no}] = roipoly;

                show_image(hObject,handles)

            case 'No'

            case 'Cancel'
        end
    else
        answr ='No';
    end

    wb = waitbar(0,'Recalculating ...');

    EF = get(handles.frame_slider,'Max');

    for f = frame_no+1:EF

        if ~isempty(find(handles.Region(k).Fascicle(j).analysed_frames == f,1))

            handles = apply_transform(handles,f,f-1,k,j);

        else

            waitbar(1,wb);
            break

        end
        waitbar((f-frame_no)/(EF-frame_no),wb);
    end


    close(wb)
    show_image(hObject,handles)

else
    warndlg('Selected Fascicle Data Not Found','WARNING')
end


% --- Executes on button press in manual_correct_keyframe.
function manual_correct_keyframe_Callback(hObject, eventdata, handles)
% hObject    handle to manual_correct_keyframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles.Region(1).Fascicle(1),'fas_x_keyframe')
    L = get(handles.keyframe_list,'String');
    f = get(handles.keyframe_list,'Value');
    frame_no = str2double(L{f});
    keyframes = cellfun(@str2double,L);

    i = get(handles.no_tracked_regions,'Value');
    j = get(handles.no_tracked_fascicles,'Value');

    % select the two points that define the fascicle line of action
    %set(handles.axes1,'XLimMode','auto')
    [p1,p2]=rbline(handles.axes1);
    handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no}(1) = p1(1);
    handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no}(2) = p2(1);
    handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no}(1) = p1(2);
    handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no}(2) = p2(2);

    if handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no}(1) > handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no}(2)
        handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no} = fliplr(handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no});
        handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no} = fliplr(handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no});
    end

    handles.Region(i).Fascicle(j).current_xy_keyframe = [handles.Region(i).Fascicle(j).fas_x{frame_no};handles.Region(i).Fascicle(j).fas_y{frame_no}]';
    handles.Region(i).fas_pen_keyframe(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no})),...
        abs(diff(handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no})));
    scalar = handles.ID;%/handles.vidHeight;
    handles.Region(i).fas_length_keyframe(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no})...
        .^2 + diff(handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no}).^2);

    DT = keyframes/handles.FrameRate;
    Dfas_x = [handles.Region(i).Fascicle(j).fas_x_keyframe{keyframes}] - [handles.Region(i).Fascicle(j).fas_x{keyframes}];
    Dfas_y = [handles.Region(i).Fascicle(j).fas_y_keyframe{keyframes}] - [handles.Region(i).Fascicle(j).fas_y{keyframes}];
    Dfas_x = reshape(Dfas_x,2,length(keyframes))';
    Dfas_y = reshape(Dfas_y,2,length(keyframes))';

    Time = (1:length(handles.Region(i).Fascicle(j).fas_x))/handles.FrameRate;

    for f = 1:length(handles.Region(i).Fascicle(j).fas_x)

        handles.Region(i).Fascicle(j).fas_x_corr{f} = handles.Region(i).Fascicle(j).fas_x{f} +...
            interp1(DT', Dfas_x, Time(f),'linear');
        handles.Region(i).Fascicle(j).fas_y_corr{f} = handles.Region(i).Fascicle(j).fas_y{f} +...
            interp1(DT', Dfas_y, Time(f),'linear');

        % calculate the length and pennation for the current frame
        handles.Region(i).fas_pen_corr(f,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y_corr{f})),...
            abs(diff(handles.Region(i).Fascicle(j).fas_x_corr{f})));
        scalar = handles.ID;%/handles.vidHeight;
        handles.Region(i).fas_length_corr(f,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y_corr{f}).^2 +...
            diff(handles.Region(i).Fascicle(j).fas_x_corr{f}).^2);

    end
    guidata(hObject, handles);
    show_image(hObject, handles)
end


% --- Executes on button press in correct_keyframes.
function correct_keyframes_Callback(hObject, eventdata, handles)
% hObject    handle to correct_keyframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ImStack') && length(get(handles.keyframe_list,'String'))>1;

    %First check which regions the user wishes to apply the correction to
    %and create a variable to manage the iterations of the for loop later
    SelectedRegions = get(handles.RegionsToCorrect,'Value');
    RegionList = get(handles.RegionsToCorrect,'String');
    ApplyToRegions = RegionList{SelectedRegions};
    if strcmp(ApplyToRegions,'all');
        allregions = get(handles.no_tracked_regions,'String');
        Rloop = 1:length(allregions);
        Rloop(isnan(Rloop))=[];
    else
        Rloop = str2double(ApplyToRegions);
        Rloop(isnan(Rloop))=[];
    end
    %Then check which fascicles the user wishes to apply the correction to
    %and create a variable to manage the iterations of the for loop later
    SelectedFascicles = get(handles.FasciclesToCorrect,'Value');
    FasList = get(handles.FasciclesToCorrect,'String');
    ApplyToFascicles = FasList{SelectedFascicles};

    if strcmp(ApplyToFascicles,'all');
        for k = 1:length(RegionList)-1
            Floop{k} = 1:length(handles.Region(k).Fascicle);
        end
    else
        for k = 1:length(RegionList)-1
            Floop{k} = str2double(ApplyToFascicles);
        end
    end

    if length(find(handles.Region(Rloop(1)).fas_length(:,Floop{1}(1)) ~= 0)) == handles.NumFrames

        contents = get(handles.keyframe_list,'String');% Get the keyframes
        image_depth = handles.ID;% get the image depth

        % store all of the original values
        Region = handles.Region;

        for f = 1:length(get(handles.keyframe_list,'String')) %loop through each keyframe

            set(handles.keyframe_list,'Value',f)% set the selected keyframe

            frame_no = str2double(contents{f});% set frame # to current keyframe

            % set the frame slider to the new frame value
            set(handles.frame_slider,'Value',frame_no);

            % set the string in the frame_number box to the current frame value
            set(handles.frame_number,'String',num2str(frame_no));

            % Load the image for the current keyframe
            Im = imcrop(handles.ImStack(:,:,frame_no+handles.start_frame-1),handles.crop_rect);
            % make the image axes the current axes
            axes(handles.axes1)
            cla
            % show the image
            imshow(Im)
            hold on
            handles.Im = Im;
            % this image is the new image
            handles.NIm = handles.Im;

            %Loop through each selected region and fascicle, checking if
            %that fascicle has already had a correction done on it. If it
            %has, make the corrected value the 'uncorrected value'
            for i = Rloop(1):Rloop(end)
                for j = Floop{i}(1):Floop{i}(end)
                    if isfield(handles.Region(i).Fascicle(j),'fas_x_corr') && ~isempty(handles.Region(i).Fascicle(j).fas_x_corr)
                        handles.Region(i).Fascicle(j).fas_x = handles.Region(i).Fascicle(j).fas_x_corr;
                        handles.Region(i).Fascicle(j).fas_y = handles.Region(i).Fascicle(j).fas_y_corr;
                        handles.Region(i).Fascicle(j).fas_pen = handles.Region(i).fas_pen_corr;
                        handles.Region(i).Fascicle(j).fas_length = handles.Region(i).fas_length_corr;
                    end
                end
            end

            % If this is the first keyframe, this is our reference value
            % and is not required to be adjusted for drift. Thus, we do not
            % adjust it and move on to the next iteration
            if f == 1

                for i = Rloop(1):Rloop(end)
                    for j = Floop{i}(1):Floop{i}(end)
                        handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no} = handles.Region(i).Fascicle(j).fas_x{frame_no};
                        handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no} = handles.Region(i).Fascicle(j).fas_y{frame_no};
                        handles.Region(i).Fascicle(j).current_xy_keyframe(1,1) = handles.Region(i).Fascicle(j).fas_x{frame_no}(1);
                        handles.Region(i).Fascicle(j).current_xy_keyframe(1,2) = handles.Region(i).Fascicle(j).fas_y{frame_no}(1);
                        handles.Region(i).Fascicle(j).current_xy_keyframe(2,1) = handles.Region(i).Fascicle(j).fas_x{frame_no}(2);
                        handles.Region(i).Fascicle(j).current_xy_keyframe(2,2) = handles.Region(i).Fascicle(j).fas_y{frame_no}(2);
                        handles.Region(i).fas_pen_keyframe(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no})),...
                            abs(diff(handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no})));
                        scalar = image_depth;%/handles.vidHeight;
                        handles.Region(i).fas_length_keyframe(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no})...
                            .^2 + diff(handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no}).^2);
                    end
                end

                % make new image the current image for next iteration
                handles.CIm = handles.NIm;
                handles.prev_frame_no = frame_no;


            else %if it's not the first keyframe then we must compute the warp between the current keyframe and the previous one
                %setup the images and the affine flow parameters
                im1 = handles.CIm;
                im2 = handles.NIm;
                points = detectMinEigenFeatures(im1,'FilterSize',11, 'MinQuality', 0.005);
                points = double(points.Location);

                for i = Rloop(1):Rloop(end)

                    [inPoints] = inpolygon(points(:,1),points(:,2), handles.Region(i).ROIx{handles.prev_frame_no}, handles.Region(i).ROIy{handles.prev_frame_no});
                    points = points(inPoints,:);
                    pointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',[21 71]);
                    initialize(pointTracker,points,im1);
                    [pointsNew, isFound] = step(pointTracker, im2);
                    [w,~] = estimateGeometricTransform2D(points(isFound,:), pointsNew(isFound,:), 'affine', 'MaxDistance',50);
                    handles.Region(i).warp(:,:,frame_no) = w;

                    for j = Floop{i}(1):Floop{i}(end)
                        %determine original values before transformation
                        handles.Region(i).Fascicle(j).current_xy = [handles.Region(i).Fascicle(j).fas_x{frame_no}' handles.Region(i).Fascicle(j).fas_y{frame_no}'];

                        % determine transformation
                        handles = apply_transform(handles,frame_no,handles.prev_frame_no,i,j);

                        % Assign variables to keyframes
                        handles.Region(i).Fascicle(j).current_xy_keyframe = [handles.Region(i).Fascicle(j).fas_x{frame_no}' handles.Region(i).Fascicle(j).fas_y{frame_no}'];
                        handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no} = handles.Region(i).Fascicle(j).fas_x{frame_no};
                        handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no} = handles.Region(i).Fascicle(j).fas_y{frame_no};
                        handles.Region(i).fas_pen_keyframe(frame_no,j) = handles.Region(i).fas_pen(frame_no,j);
                        handles.Region(i).fas_length_keyframe(frame_no,j) = handles.Region(i).fas_length(frame_no,j);

                        plot(handles.axes1,handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no},handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no},'y','LineWidth',2);

                    end

                end
            end

            handles.CIm = handles.NIm;
            handles.prev_frame_no = frame_no;

        end


        %reset the old values so we can compare to keyframes
        for f = 1:length(get(handles.keyframe_list,'String'))
            frame_no = str2double(contents{f});
            for i = Rloop(1):Rloop(end)
                handles.Region(i).ROIx{frame_no} = Region(i).ROIx{frame_no};
                handles.Region(i).ROIy{frame_no} = Region(i).ROIy{frame_no};
                for j = Floop{i}(1):Floop{i}(end)
                    handles.Region(i).Fascicle(j).fas_x{frame_no} = Region(i).Fascicle(j).fas_x{frame_no};
                    handles.Region(i).Fascicle(j).fas_y{frame_no} = Region(i).Fascicle(j).fas_y{frame_no};
                    % calculate the difference between the keyframe x & y coordinates
                    % and the original tracked x & y coordinates
                    Region(i).Fascicle(j).DT(f) = frame_no/handles.FrameRate;
                    Region(i).Fascicle(j).Dfas_x(f,:) = handles.Region(i).Fascicle(j).fas_x_keyframe{frame_no} - handles.Region(i).Fascicle(j).fas_x{frame_no};
                    Region(i).Fascicle(j).Dfas_y(f,:) = handles.Region(i).Fascicle(j).fas_y_keyframe{frame_no} - handles.Region(i).Fascicle(j).fas_y{frame_no};

                end
            end

        end

        % Create variables representing drift over time
        for i = Rloop(1):Rloop(end)
            for j = Floop{i}(1):Floop{i}(end)
                Region(i).Fascicle(j).DT = [0 Region(i).Fascicle(j).DT get(handles.frame_slider,'Max')/handles.FrameRate];
                Region(i).Fascicle(j).Dfas_x = [0 0; Region(i).Fascicle(j).Dfas_x; Region(i).Fascicle(j).Dfas_x(end,:)];
                Region(i).Fascicle(j).Dfas_y = [0 0; Region(i).Fascicle(j).Dfas_y; Region(i).Fascicle(j).Dfas_y(end,:)];

                % if the last keyframe is the last frame, then this frame
                % will have been duplicated and we need to remove the last
                % value from the above variables or it will cause an error
                % in the interpolation.
                if Region(i).Fascicle(j).DT(end)<=Region(i).Fascicle(j).DT(end-1)
                    Region(i).Fascicle(j).DT(end)=[];
                    Region(i).Fascicle(j).Dfas_x(size(Region(i).Fascicle(j).Dfas_x,1),:)=[];
                    Region(i).Fascicle(j).Dfas_y(size(Region(i).Fascicle(j).Dfas_y,1),:)=[];
                end
            end
        end

        % Now we can make the correction. Drift is assumed to be linear
        % over time and is subtracted from the x,y coordinates of the
        % fascicle in each frame
        for i = Rloop(1):Rloop(end)
            for j = Floop{i}(1):Floop{i}(end)

                Time = (1:length(handles.Region(i).Fascicle(j).fas_x))/handles.FrameRate;

                for f = 1:length(handles.Region(i).Fascicle(j).fas_x)

                    handles.Region(i).Fascicle(j).fas_x_corr{f} = handles.Region(i).Fascicle(j).fas_x{f} +...
                        interp1(Region(i).Fascicle(j).DT', Region(i).Fascicle(j).Dfas_x, Time(f),'linear');
                    handles.Region(i).Fascicle(j).fas_y_corr{f} = handles.Region(i).Fascicle(j).fas_y{f} +...
                        interp1(Region(i).Fascicle(j).DT', Region(i).Fascicle(j).Dfas_y, Time(f),'linear');

                    % calculate the length and pennation for the current frame
                    handles.Region(i).fas_pen_corr(f,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y_corr{f})),...
                        abs(diff(handles.Region(i).Fascicle(j).fas_x_corr{f})));
                    scalar = image_depth;%/handles.vidHeight;
                    handles.Region(i).fas_length_corr(f,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y_corr{f}).^2 +...
                        diff(handles.Region(i).Fascicle(j).fas_x_corr{f}).^2);

                end
                clear Time
            end
        end


    end

    handles.KeyframeInd = logical(1.0);
    show_image(hObject, handles)

end

%---------------------------------------------------------
% Function to show image with appropriate image processing
%---------------------------------------------------------
function show_image(hObject, handles)

set(gcf,'DoubleBuffer','on');

if isfield(handles,'ImStack')

    axes(handles.length_plot)
    cla
    hold on
    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'));
    % get the image depth
    image_depth = handles.ID;

    % show the image
    if isfield(handles,'ImStack')
        Im = handles.ImStack(:,:,frame_no+handles.start_frame-1);
        % make the image axes the current axes
        axes(handles.axes1)
        %cla(handles.axes1)
        % crop image with current cropping value
        Im = imcrop(Im,handles.crop_rect);
        hold off
        % show the image
        imshow(Im,'Border','loose')
        colormap(gray(256)); axis off;
        if isfield(handles,'Vidax_X')
            set(handles.axes1,'XLimMode','manual','XLim',handles.Vidax_X,'YLim',handles.Vidax_Y');
        else
            set(handles.axes1,'XLimMode','manual','XLim',[0-(0.1*handles.crop_rect(3)) handles.crop_rect(3)+(0.1*handles.crop_rect(3))],'YLimMode','manual','YLim',[0 handles.vidHeight])
        end

        hold on
    end

    % Update Im and NIm
    handles.Im = Im;
    handles.NIm = handles.Im;

    if isfield(handles,'Region')
        % go through each defined muscle region, check if the frame has been processed
        % and, if not, perform the affine flow calculation and apply to any defined fascicles
        if isfield(handles.Region,'ROI')

            % define the parameters of the LKT / affine algorithm
            % ADD HERE

            for i = 1:length(handles.Region)

                if isfield(handles.Region(i),'Fascicle')

                    if ~isempty(handles.Region(i).Fascicle)

                        for j = 1:length(handles.Region(i).Fascicle)

                            if isempty(find(handles.Region(i).Fascicle(j).analysed_frames==frame_no, 1))

                                %setup current and new image
                                im1 = handles.CIm;
                                im2 = handles.NIm;

                                points = detectMinEigenFeatures(im1,'FilterSize',11, 'MinQuality', 0.005);
                                points = double(points.Location);
                                [inPoints] = inpolygon(points(:,1),points(:,2), handles.Region(i).ROIx{handles.prev_frame_no}, handles.Region(i).ROIy{handles.prev_frame_no});
                                points = points(inPoints,:);
                                pointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',[21 71]);
                                initialize(pointTracker,points,im1);
                                [pointsNew, isFound] = step(pointTracker, im2);
                                [w,~] = estimateGeometricTransform2D(points(isFound,:), pointsNew(isFound,:), 'affine', 'MaxDistance',50);
                                handles.Region(i).warp(:,:,frame_no) = w;

                                handles.CIm = handles.NIm;

                                handles = apply_transform(handles,frame_no,handles.prev_frame_no,i,j);
                            end
                        end

                    else
                        for j = 1:length(handles.Region(i).Fascicle)
                            if isempty(find(handles.Region(i).Fascicle(j).analysed_frames==frame_no, 1))
                                handles.Region(i).Fascicle(j).current_xy(1,1) = handles.Region(i).Fascicle(j).fas_x{frame_no}(1);
                                handles.Region(i).Fascicle(j).current_xy(1,2) = handles.Region(i).Fascicle(j).fas_y{frame_no}(1);
                                handles.Region(i).Fascicle(j).current_xy(2,1) = handles.Region(i).Fascicle(j).fas_x{frame_no}(2);
                                handles.Region(i).Fascicle(j).current_xy(2,2) = handles.Region(i).Fascicle(j).fas_y{frame_no}(2);
                                % calculate the length and pennation for the current frame
                                handles.Region(i).fas_pen(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y{frame_no})),abs(diff(handles.Region(i).Fascicle(j).fas_x{frame_no})));
                                scalar = image_depth;%/handles.vidHeight;
                                handles.Region(i).fas_length(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{frame_no}).^2 + diff(handles.Region(i).Fascicle(j).fas_x{frame_no}).^2);
                            end
                        end

                    end



                    if isfield(handles.Region(i),'fas_length')
                        if ~isempty(handles.Region(i).fas_length)
                            %plot the fascicle length on the graph relative to time
                            % define the line type
                            switch i
                                case 1
                                    c = '-';
                                case 2
                                    c = ':';
                                case 3
                                    c = '-.';
                                otherwise
                                    c = '--';
                            end

                            dt = 1/handles.FrameRate;
                            nz = logical(handles.Region(i).fas_length(:,1) ~= 0);
                            FL = handles.Region(i).fas_length(nz,:);
                           
                            time = (dt:dt:length(handles.Region(i).fas_length(:,1))*dt)+((handles.start_frame-1)*dt);

                            for j = 1:length(handles.Region(i).Fascicle)

                                if isfield(handles.Region(i).Fascicle(j),'fas_x_corr') && ~isempty(handles.Region(i).Fascicle(j).fas_x_corr)
                                    nz = logical(handles.Region(i).fas_length_corr(:,j) ~= 0);
                                    time = (dt:dt:length(handles.Region(i).fas_length_corr(:,1))*dt)+((handles.start_frame-1)*dt);
                                    %plot the fascicle length on the graph relative to time
                                    if isfield(handles.Region(i),'fas_length_keyframe')
                                        dt = 1/handles.FrameRate;
                                        nnz = logical(handles.Region(i).fas_length_keyframe(:,j) ~= 0);
                                        FL_keyframe = handles.Region(i).fas_length_keyframe(nnz,j);
                                        time_keyframe = (dt:dt:length(handles.Region(i).fas_length_keyframe(:,j))*dt)+((handles.start_frame-1)*dt);
                                        %plot(handles.length_plot,time(nz),handles.Region(i).fas_length_corr(nz,j),c,time_keyframe(nnz),FL_keyframe,'ro');
                                        plot(handles.length_plot,time(nz),handles.Region(i).fas_length_corr(nz,j),time_keyframe(nnz),FL_keyframe,'ro');
                                    else
                                        plot(handles.length_plot,time(nz),handles.Region(i).fas_length_corr(nz,j),'r'); hold on
                                        plot(handles.length_plot,time(nz),handles.Region(i).fas_length_Hough(nz,j),'y')
                                    end

                                    axes(handles.length_plot)

                                    if size(handles.Region(i).fas_length_corr,1)>=frame_no

                                        text(handles.Time(frame_no+handles.start_frame-1)+0.2,handles.Region(i).fas_length_corr(frame_no,j)...
                                            ,num2str(round(handles.Region(i).fas_length_corr(frame_no,j))),...
                                            'FontSize',10,'FontWeight','bold')
                                    end

                                    %for j = 1:length(handles.Region(i).Fascicle)
                                    %plot the fascicle
                                    plot(handles.axes1,handles.Region(i).Fascicle(j).fas_x_corr{frame_no},handles.Region(i).Fascicle(j).fas_y_corr{frame_no},'yo-','LineWidth',2)
                                    %end
                                else

                                    
                                    
                                    if isfield(handles.Region(i), 'fas_length_Hough')
                                        Hough_time = (dt:dt:length(handles.Region(i).fas_length_Hough(:,1))*dt)+((handles.start_frame-1)*dt);
                                        plot(handles.length_plot,Hough_time(:), handles.Region(i).fas_length_Hough(:,j),'y')
                                    end

                                    plot(handles.length_plot,time(nz),FL,'r','linewidth',2); hold on
                                    set(handles.length_plot,'ylim',[min(FL)*0.85 max(FL)*1.15]); %set axis 15% difference of min and and value,easier to read

                                    if isfield(handles.Region(i),'fas_length_original')
                                    plot(handles.length_plot,time(nz), handles.Region(i).fas_length_original(nz),'r:','linewidth',1); hold on
                                    end
                                    
                                    if isfield(handles.Region(i),'adjusted_flength')&& size(handles.Region(i).adjusted_flength,1) > 1
                                        Adj_FL = handles.Region(i).adjusted_flength(nz,:);
                                        plot(handles.length_plot,time(nz),Adj_FL)
                                    end

                                    axes(handles.length_plot)


                                    if size(FL,1)>=frame_no
                                        text(handles.Time(frame_no+handles.start_frame-1)+0.2,FL(frame_no,j)...
                                            ,num2str(round(FL(frame_no,j))),...
                                            'FontSize',10,'FontWeight','bold')
                                    end


                                    %for j = 1:length(handles.Region(i).Fascicle)

                                    %plot the fascicle
                                    if isfield(handles.Region(i),'fas_length_original')
                                        plot(handles.axes1,handles.Region(i).Fascicle(j).fas_x_original{frame_no},handles.Region(i).Fascicle(j).fas_y_original{frame_no},'ro:','LineWidth',1)
                                    end
                                    
                                    plot(handles.axes1,handles.Region(i).Fascicle(j).fas_x{frame_no},handles.Region(i).Fascicle(j).fas_y{frame_no},'ro-','LineWidth',2)
                                    
                                    if isfield(handles.Region(i).Fascicle(j),'fas_xy_Hough')
                                        if ~isempty(handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no})


                                            handles.Region(i).Fascicle(j).Hough_fas_x{frame_no}(1) = handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(1,1);
                                            handles.Region(i).Fascicle(j).Hough_fas_y{frame_no}(1) = handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(1,2);
                                            handles.Region(i).Fascicle(j).Hough_fas_x{frame_no}(2) = handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(2,1);
                                            handles.Region(i).Fascicle(j).Hough_fas_y{frame_no}(2) = handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(2,2);

                                            plot(handles.axes1,handles.Region(i).Fascicle(j).Hough_fas_x{frame_no},handles.Region(i).Fascicle(j).Hough_fas_y{frame_no},'yo-','LineWidth',2)
                                        end
                                    end
                                    
                                    if isfield(handles.Region(i),'APOROIy')
                                        if ~isempty(handles.Region(i).APOROIy{frame_no})
                                            ROIcor = [handles.Region(i).APOROIx{frame_no}, handles.Region(i).APOROIy{frame_no}(:)];
                                            plot(handles.axes1,ROIcor(:,1), ROIcor(:,2),'y--','LineWidth',2)
                                        end
                                    end
                                end
                            end

                            set(handles.length_plot,'XLim',[0 handles.Time(end)]);
                        end
                    end

                end

                % plot the ROI
                plot(handles.axes1,handles.Region(i).ROIx{frame_no},handles.Region(i).ROIy{frame_no},'r:','LineWidth',2)

            end

        end
    end

    axes(handles.length_plot)
    set(handles.length_plot,'XLim',[0 handles.Time(end)]);
    L_range = get(handles.length_plot,'YLim');
    plot([handles.Time(frame_no+handles.start_frame-1) handles.Time(frame_no+handles.start_frame-1)],[L_range(1) L_range(2)],'k');
    hold off
    

    xlabel(handles.length_plot,'Time')
    ylabel(handles.length_plot,'Fascicle Length (mm)')

    % plot the data from the variable list of the accompanying dat file
    if isfield(handles,'D')
        plot_data(hObject,handles)
        axes(handles.mat_plot)
        y_range = get(handles.mat_plot,'YLim');
        hold on
        plot([handles.Time(frame_no+handles.start_frame-1) handles.Time(frame_no+handles.start_frame-1)],...
            [y_range(1) y_range(2)],'k')

        hold off
    end
    set(handles.mat_plot,'XLim',[0 handles.Time(end)]);

    xlabel(handles.mat_plot,'Time')
    ylabel(handles.mat_plot,'Data')
    handles.prev_frame_no = frame_no;

    % Update handles structure
    guidata(hObject, handles);



end

% --- Executes on button press in process_all.
function process_all_Callback(hObject, eventdata, handles)
% hObject    handle to process_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% all frames
frames = 1:handles.NumFrames;

% if there exists a ApoDetect_downsample estimate
if strcmp(handles.ApoDetect_downsample.String,'0')
    handles.ApoDetect_frames = [];
else 
    handles.ApoDetect_frames = 1:str2double(handles.ApoDetect_downsample.String):handles.NumFrames;
end

% if there exists a ApoDetect_downsample estimate
if strcmp(handles.TimTrack_downsample.String,'0')
    handles.Hough_frames = [];
else 
    handles.Hough_frames = 1:str2double(handles.TimTrack_downsample.String):handles.NumFrames;
end

parms = get_TimTrack_parms();

if isfield(handles,'ImStack')
    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'));

        h = waitbar(0,'Please wait while processing forward...');

        for i = 1:length(handles.Region)

            % TimTrack
            for f = frames

                    im2 = imcrop(handles.ImStack(:,:,f),handles.crop_rect);
                    handles.NIm = im2;

                    [geofeatures(f), apovecs(f)] = auto_ultrasound(im2, parms);
                    
                    % optional: fix to edges
                    APOROI(:,1) = [1; 1; handles.vidWidth+1; handles.vidWidth+1; 1];

                    % make sure ROIs are on aponeuroses
                    APOROI([2,3],2) = polyval(geofeatures(f).deep_coef, APOROI([2,3],1));
                    APOROI([1,4,5],2) = polyval(geofeatures(f).super_coef, APOROI([1,4,5],1));

                    [handles.Region(i).APOROI{frame_no},handles.Region(i).APOROIx{frame_no}, handles.Region(i).APOROIy{frame_no}] = roipoly(handles.NIm, APOROI(:,1),APOROI(:,2));
            end
        end
                    
            

            % setup current and new image
            points = detectMinEigenFeatures(im1,'FilterSize',11, 'MinQuality', 0.005);
            points = double(points.Location);
            [inPoints] = inpolygon(points(:,1),points(:,2), handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no});
            points = points(inPoints,:);
            pointTracker = vision.PointTracker('NumPyramidLevels',4,'MaxIterations',50,'MaxBidirectionalError',inf,'BlockSize',[21 71]);
            initialize(pointTracker,points,im1);


            for f = frame_no+1:get(handles.frame_slider,'Max')
                % Get the current image

                im2 = imcrop(handles.ImStack(:,:,handles.start_frame+f-1),handles.crop_rect);
                handles.NIm = im2;

                %Compute the flow and new roi
                [pointsNew, isFound] = step(pointTracker, im2);
                [w,~] = estimateGeometricTransform2D(points(isFound,:), pointsNew(isFound,:), 'affine', 'MaxDistance',50);
                handles.Region(i).warp(:,:,f) = w;

                handles.CIm = handles.NIm;

                for j = 1:length(handles.Region(i).Fascicle)

                    handles = apply_transform(handles,f,f-1,i,j);
                    
                    % correct ROI
                    if isfield(handles.Region, 'APOROIy')
                        % correct
                        alpha = str2double(handles.ApoGain.String);
                        if isnan(alpha)
                            alpha = 0.5;
                        end
                        
                        ROIpos = [handles.Region(i).ROIx{f}, handles.Region(i).ROIy{f}];
                        ROIcor = [handles.Region(i).APOROIx{f}, handles.Region(i).APOROIy{f}(:)];
                            
                        ROIpos = ROIpos + alpha * (ROIcor - ROIpos);

                        [handles.Region(i).ROI{f},handles.Region(i).ROIx{f}, handles.Region(i).ROIy{f}] = roipoly(handles.NIm, ROIpos(:,1),ROIpos(:,2));
                    end


                    points = detectMinEigenFeatures(im2,'FilterSize',11, 'MinQuality', 0.005);
                    points = double(points.Location);
                    [inPoints] = inpolygon(points(:,1),points(:,2),handles.Region(i).ROIx{f},handles.Region(i).ROIy{f});
                    points = points(inPoints,:);
                    setPoints(pointTracker, points);

                end

                waitbar((f+(get(handles.frame_slider,'Max')*(i-1)))/...
                    (get(handles.frame_slider,'Max')*length(handles.Region)),h)
            end
end
        close(h)

    % update the image axes using show_image function (bottom)
    show_image(hObject,handles)

% Update handles structure
guidata(hObject, handles);


% function plot_data(hObject,handles)
% 
% axes(handles.mat_plot); cla;
% 
% plot([0 1], [0 1])
% 
% guidata(hObject,handles);
 
 
% --- Function to plot the data variables from spike mat, c3d or biodex dat
% file
% function plot_data(hObject,handles)
% 
% if isfield(handles,'D')
% 
%     % determine the variables highlighted in list
%     cur_var = get(handles.variable_list,'Value');
% 
% 
%     frame_no = round(get(handles.frame_slider,'Value'));
%     frametime = handles.Time(frame_no+handles.start_frame-1);
% 
%     axes(handles.mat_plot); cla;
% 
%     % define the colours to plot
%     %col = {'b','k','r','g','y','m','k'};
% 
%     hold on;
% 
%     % plot the data in the selected cells against their corresponding time
%     % cell
%     cellfun(@plot,{handles.D(cur_var).time},{handles.D(cur_var).value});%,col(1:length(cur_var)));
% 
%     % make plot legend - note some variable have more than one plot in them
%     % and therefore we need to take this into account
%     % only put legend on plot if the legend check box is true
%     if get(handles.legend_check,'Value')
%         cur_legend = [];
%         for i = 1:length(cur_var)
%             for j = 1:size(handles.D(cur_var(i)).value,2)
%                 if strcmp(handles.D(cur_var(i)).type,'Marker')
%                     switch j
%                         case 1
%                             M = 'X';
%                         case 2
%                             M = 'Y';
%                         case 3
%                             M = 'Z';
%                     end
%                     cur_legend = [cur_legend {[handles.D(cur_var(i)).name ' ' M]}];
%                 else cur_legend = [cur_legend {handles.D(cur_var(i)).name}];
%                 end
%                 y_range = get(handles.mat_plot,'YLim');
%                 [frametimeVar,frameVar] = min(abs(frametime - handles.D(cur_var(i)).time));
%                 text(handles.D(cur_var(i)).time(frameVar)+0.2,handles.D(cur_var(i)).value(frameVar,j),...
%                     num2str(handles.D(cur_var(i)).value(frameVar,j)),...
%                     'FontSize',10,'FontWeight','bold');
%             end
%         end
%         legend(cur_legend)
%     else legend off
%     end
% 
%     if isfield(handles,'c3d_events')
% 
%         elist = get(handles.EventList,'String');
%         enum  = get(handles.EventList,'Value');
%         cur_ev = elist{enum}; %ultrasound_tracking_v2_5
%         eventT = handles.c3d_events.(cur_ev);
% 
%         y_vals = get(handles.mat_plot,'YLim');
%         for k = 1:length(eventT)
%             eplot{k} = plot([eventT(k),eventT(k)],y_vals,':','Color',[0.4 0.4 0.4]);
%         end
% 
% 
%     end
% 
% 
%     hold off;
% 
%     guidata(hObject,handles);
% 
% end

function handles = detect_ROIs(handles,frame_no,i,j)

% if ~isnan(handles.ApoDetect_frames(1))
%     if any(handles.ApoDetect_frames == frame_no)

parms.apo.frangi = [18 20];
parms.apo.apox = round(linspace(1, handles.vidWidth, 10));

parms.apo.method = 'Hough';
parms.apo.th = 0.5;
parms.apo.show = 0;

parms.apo.deep.cut = handles.Hough.apo.deep.cut;
parms.apo.super.cut = handles.Hough.apo.super.cut;

parms.fas.filter = 0;

% filter and fit aponeurosis
[~, super_obj, deep_obj] = filter_usimage(handles.NIm,parms);
super_coef = polyfit(parms.apo.apox,super_obj, 1);
deep_coef = polyfit(parms.apo.apox,deep_obj, 1);

% optional: fix to edges
APOROI(:,1) = [1; 1; handles.vidWidth+1; handles.vidWidth+1; 1];

% make sure ROIs are on aponeuroses
APOROI([2,3],2) = polyval(deep_coef, APOROI([2,3],1));
APOROI([1,4,5],2) = polyval(super_coef, APOROI([1,4,5],1));

[handles.Region(i).APOROI{frame_no},handles.Region(i).APOROIx{frame_no}, handles.Region(i).APOROIy{frame_no}] = roipoly(handles.NIm, APOROI(:,1),APOROI(:,2));


% function to perform transformation (warp) of the points
function handles = apply_transform(handles,frame_no,prev_frame_no,i,j)

w = handles.Region(i).warp(:,:,frame_no);

%% ROI
ROIpos = transformPointsForward(w, [handles.Region(i).ROIx{prev_frame_no} handles.Region(i).ROIy{prev_frame_no}]);
[handles.Region(i).ROI{frame_no},handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no}] = roipoly(handles.NIm,ROIpos(:,1),ROIpos(:,2));

% Adjust the ROI if it goes outside the image
handles.Region(i).ROIx{frame_no}(handles.Region(i).ROIx{frame_no} > handles.vidWidth) = handles.vidWidth;
handles.Region(i).ROIy{frame_no}(handles.Region(i).ROIy{frame_no} > handles.vidHeight) = handles.vidHeight;
handles.Region(i).ROIx{frame_no}(handles.Region(i).ROIx{frame_no} < 1) = 1;
handles.Region(i).ROIy{frame_no}(handles.Region(i).ROIy{frame_no} < 1) = 1;

% store aponeurosis coefficients
APOdeepx = handles.Region(i).ROIx{frame_no}(2:3);
APOdeepy = handles.Region(i).ROIy{frame_no}(2:3);
APOsuperx = handles.Region(i).ROIx{frame_no}([1,4]);
APOsupery = handles.Region(i).ROIy{frame_no}([1,4]);

handles.super_coef{frame_no} = polyfit(APOsuperx,APOsupery, 1);
handles.deep_coef{frame_no} = polyfit(APOdeepx,APOdeepy, 1);

%% Fascicle
% loop through all fasicles defined for the region and apply
% the warp and calculate new fascicle length
if isfield(handles.Region(i).Fascicle(j),'fas_x_corr') && ~isempty(handles.Region(i).Fascicle(j).fas_x_corr{prev_frame_no})
    handles.Region(i).Fascicle(j).current_xy(1,1) = handles.Region(i).Fascicle(j).fas_x_corr{prev_frame_no}(1);
    handles.Region(i).Fascicle(j).current_xy(1,2) = handles.Region(i).Fascicle(j).fas_y_corr{prev_frame_no}(1);
    handles.Region(i).Fascicle(j).current_xy(2,1) = handles.Region(i).Fascicle(j).fas_x_corr{prev_frame_no}(2);
    handles.Region(i).Fascicle(j).current_xy(2,2) = handles.Region(i).Fascicle(j).fas_y_corr{prev_frame_no}(2);
    handles.Region(i).Fascicle(j).current_xy = transformPointsForward(w, handles.Region(i).Fascicle(j).current_xy);
    handles.Region(i).Fascicle(j).fas_x_corr{frame_no}(1) = handles.Region(i).Fascicle(j).current_xy(1,1);
    handles.Region(i).Fascicle(j).fas_y_corr{frame_no}(1) = handles.Region(i).Fascicle(j).current_xy(1,2);
    handles.Region(i).Fascicle(j).fas_x_corr{frame_no}(2) = handles.Region(i).Fascicle(j).current_xy(2,1);
    handles.Region(i).Fascicle(j).fas_y_corr{frame_no}(2) = handles.Region(i).Fascicle(j).current_xy(2,2);
else % in case there is no corrected version (normal case)
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
    
    % more book keeping
    handles.Region(i).Fascicle(j).fas_x_original{frame_no}(1) = handles.Region(i).Fascicle(j).current_xy(1,1);
    handles.Region(i).Fascicle(j).fas_y_original{frame_no}(1) = handles.Region(i).Fascicle(j).current_xy(1,2);
    handles.Region(i).Fascicle(j).fas_x_original{frame_no}(2) = handles.Region(i).Fascicle(j).current_xy(2,1);
    handles.Region(i).Fascicle(j).fas_y_original{frame_no}(2) = handles.Region(i).Fascicle(j).current_xy(2,2);
end

%% update fascicle with aponeurosis
% fit fascicle
handles.fas_coef{frame_no} = polyfit(handles.Region(i).Fascicle(j).fas_x{frame_no}, handles.Region(i).Fascicle(j).fas_y{frame_no},1);

Deep_intersect_x = round((handles.deep_coef{frame_no}(2) - handles.fas_coef{frame_no}(2)) ./ ...
           (handles.fas_coef{frame_no}(1) - handles.deep_coef{frame_no}(1)));

Super_intersect_x = round((handles.super_coef{frame_no}(2) - handles.fas_coef{frame_no}(2)) ./ ...
            (handles.fas_coef{frame_no}(1) - handles.super_coef{frame_no}(1)));

% find intersect    
% Deep_intersect_x = round(fzero(@(x) polyval(handles.deep_coef{frame_no}     -handles.fas_coef{frame_no},x),0));
% Super_intersect_x = round(fzero(@(x) polyval(handles.super_coef{frame_no}   -handles.fas_coef{frame_no},x),0));

Super_intersect_y = polyval(handles.super_coef{frame_no}, Super_intersect_x); 
Deep_intersect_y = polyval(handles.deep_coef{frame_no}, Deep_intersect_x); 

% update
handles.Region(i).Fascicle(j).fas_x{frame_no}(1) = Deep_intersect_x;
handles.Region(i).Fascicle(j).fas_y{frame_no}(1) = Deep_intersect_y;
handles.Region(i).Fascicle(j).fas_x{frame_no}(2) = Super_intersect_x;
handles.Region(i).Fascicle(j).fas_y{frame_no}(2) = Super_intersect_y;

%% Fascicle length
% calculate the length and pennation for the current frame
if isfield(handles.Region(i).Fascicle(j),'fas_x_corr') && ~isempty(handles.Region(i).Fascicle(j).fas_x_corr{frame_no})
    handles.Region(i).fas_pen_corr(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y_corr{frame_no})),abs(diff(handles.Region(i).Fascicle(j).fas_x_corr{frame_no})));
    scalar = handles.ID;%/handles.vidHeight;
    handles.Region(i).fas_length_corr(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y_corr{frame_no}).^2 + diff(handles.Region(i).Fascicle(j).fas_x_corr{frame_no}).^2);
else
    handles.Region(i).fas_pen(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y{frame_no})),abs(diff(handles.Region(i).Fascicle(j).fas_x{frame_no})));
    scalar = handles.ID;%/handles.vidHeight;
    handles.Region(i).fas_length(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y{frame_no}).^2 + diff(handles.Region(i).Fascicle(j).fas_x{frame_no}).^2);
end

% original values
handles.Region(i).fas_pen_original(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y_original{frame_no})),abs(diff(handles.Region(i).Fascicle(j).fas_x_original{frame_no})));
scalar = handles.ID;%/handles.vidHeight;
handles.Region(i).fas_length_original(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y_original{frame_no}).^2 + diff(handles.Region(i).Fascicle(j).fas_x_original{frame_no}).^2);

handles.Region(i).Fascicle(j).analysed_frames = sort([handles.Region(i).Fascicle(j).analysed_frames frame_no]);

% TimTrack function
function handles = do_TimTrack(handles, frame_no,i,j)
  % cropped image
        data = handles.ImStack(handles.crop_rect(2):(handles.crop_rect(2)+handles.crop_rect(4)-1), handles.crop_rect(1):(handles.crop_rect(1)+handles.crop_rect(3)-1), frame_no);
        
        % filter fascicle 
        % parameters
        %automatic load parms independently of the current location to
        %prevent errors
        tmp = mfilename; 
        tmp = erase(which(tmp),[tmp '.m']);
        load([tmp 'updates/parms.mat'])
        
        %cd('C:\Users\timvd\Documents\ultrasound-automated-algorithm\Parameters')
        %load('parms.mat')
        
        fas_filt = FrangiFilter2D(double(data), parms.fas.frangi);  % filtered fascicle image
        
        [n,m] = size(data);
        fas_cut = fas_filt .* handles.Region.ROI{frame_no}(1:n,1:m) .* handles.Region.fasROI;
        
        % Threshold fascicle image
        fas_thres = imbinarize(fas_cut,'adaptive','sensitivity', .1); % threshold
        
        % find TimTrack_downsample angle        
        parms.fas.Emask_radius = handles.Region.fasROI_radius;
        parms.fas.Emask_angle = handles.Region.fasROI_angle;
        parms.fas.thetares = .5;
        parms.fas.show = 0;
        parms.fas.npeaks = 10;
        
        parms.npeaks = handles.Hough.npeaks;
        parms.fas.range = handles.Hough.range;

        [Hough_angles, ws, H, G] = dohough(fas_thres,parms.fas);
        handles.Hough.Hough_angle{frame_no} = weightedMedian(Hough_angles,ws);
        handles.Hough.Hough_matrix{frame_no} = H;
        handles.Hough.Hough_gamma{frame_no} = G;
        handles.Hough.Hough_angles{frame_no} = Hough_angles;
        
        % obtain aponeurosis fits
        APOsupx = handles.Region.ROIx{frame_no}([1 4]);
        APOdeepx = handles.Region.ROIx{frame_no}(2:3);
        APOsupy = handles.Region.ROIy{frame_no}([1 4]);
        APOdeepy = handles.Region.ROIy{frame_no}(2:3);

%         super_coef = polyfit(APOsupx, APOsupy, 1);
%         deep_coef = polyfit(APOdeepx, APOdeepy, 1);
        
        % different for 1st frame and other frames
        if frame_no > 0
            [~,m,~] = size(fas_thres);
            
            
            costfunc = @(p,alpha, super_coef, deep_coef, image_width) (1-round(fzero(@(x) polyval(deep_coef -[-tand(alpha) p],x),0))).^2 + ... 
                                                          (image_width - round(fzero(@(x) polyval(super_coef-[-tand(alpha) p],x),0))).^2;

            M(1) = round(m/2);
            M(2) = mean([polyval(handles.deep_coef{frame_no}, M(1)) polyval(handles.super_coef{frame_no}, M(1))]);

            Pinit =  M(2) - M(1) * -tand(handles.Hough.Hough_angle{frame_no});
                                                      
            P = fminsearch(@(p) costfunc(p, handles.Hough.Hough_angle{frame_no}, handles.super_coef{frame_no},handles.deep_coef{frame_no},m), Pinit);

            fas_coef(1) = -tand(handles.Hough.Hough_angle{frame_no});
            fas_coef(2) = P;
           
            
        else

            % update current xy
            handles.Region(i).Fascicle(j).current_xy = [handles.Region(i).Fascicle(j).fas_x{frame_no};handles.Region(i).Fascicle(j).fas_y{frame_no}]';

            % determine midpoint on fascicle
            M = [mean([handles.Region(i).Fascicle(j).current_xy(1,1) handles.Region(i).Fascicle(j).current_xy(2,1)]), ...
                 mean([handles.Region(i).Fascicle(j).current_xy(1,2) handles.Region(i).Fascicle(j).current_xy(2,2)])];

                fas_coef(1) = -tand(handles.Hough.Hough_angle{frame_no});
                fas_coef(2) =  M(2) - M(1) * fas_coef(1);
        end
        
%         handles.fas_coef{frame_no} = fas_coef;

%         Deep_intersect_x = round(fzero(@(x) polyval(deep_coef(:)-fas_coef(:),x),0));
%         Super_intersect_x = round(fzero(@(x) polyval(super_coef(:)-fas_coef(:),x),0));

        Deep_intersect_x = round((handles.deep_coef{frame_no}(2) - fas_coef(2)) ./ ...
                   (fas_coef(1) - handles.deep_coef{frame_no}(1)));
               
        Super_intersect_x = round((handles.super_coef{frame_no}(2) - fas_coef(2)) ./ ...
                    (fas_coef(1) - handles.super_coef{frame_no}(1)));
        
        Super_intersect_y = polyval(handles.super_coef{frame_no}, Super_intersect_x); 
        Deep_intersect_y = polyval(handles.deep_coef{frame_no}, Deep_intersect_x); 
  
        
        % assign to variable        
        handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(1,1) = Deep_intersect_x;
        handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(1,2) = Deep_intersect_y;
        handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(2,1) = Super_intersect_x;
        handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(2,2) = Super_intersect_y;

        handles.Region(i).Fascicle(j).fas_x_Hough{frame_no}(1) = handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(1,1);
        handles.Region(i).Fascicle(j).fas_y_Hough{frame_no}(1) = handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(1,2);
        handles.Region(i).Fascicle(j).fas_x_Hough{frame_no}(2) = handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(2,1);
        handles.Region(i).Fascicle(j).fas_y_Hough{frame_no}(2) = handles.Region(i).Fascicle(j).fas_xy_Hough{frame_no}(2,2);
        
         if frame_no == 1
            handles.Region(i).Fascicle(j).fas_x{frame_no} = handles.Region(i).Fascicle(j).fas_x_Hough{frame_no};
            handles.Region(i).Fascicle(j).fas_y{frame_no} = handles.Region(i).Fascicle(j).fas_y_Hough{frame_no};
         end
        
        % calculate the length and pennation for the current frame
        handles.Region(i).fas_pen_Hough(frame_no,j) = atan2(abs(diff(handles.Region(i).Fascicle(j).fas_y_Hough{frame_no})),...
            abs(diff(handles.Region(i).Fascicle(j).fas_x_Hough{frame_no})));
        
        scalar = handles.ID;%/handles.vidHeight;
        
        handles.Region(i).fas_length_Hough(frame_no,j) = scalar*sqrt(diff(handles.Region(i).Fascicle(j).fas_y_Hough{frame_no}).^2 +...
            diff(handles.Region(i).Fascicle(j).fas_x_Hough{frame_no}).^2);
        
        if frame_no == 1
            handles.Region(i).fas_length_original(frame_no,j) = handles.Region(i).fas_length_Hough(frame_no,j);
%             handles.Region(i).fas_length_original(frame_no,j) = handles.Region(i).fas_length_Hough(frame_no,j);
        end

% --- Executes on button press in Auto_Detect.
function Auto_Detect_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_Detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'ImStack')

    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'));

    % show the image
    handles.IM = handles.ImStack(:,:,frame_no+handles.start_frame-1);

    ii = get(handles.no_tracked_regions,'Value');
    
    if ii == 1 || 2

        img = im2double(handles.IM);
        axes(handles.axes1); hold off

        %% Aponeurosis detection
        for jj = 1:2
            [~,apCentre] = ginputYellow(1);

            apCentre = round(apCentre/handles.vidHeight,2)*100;
            apRound(jj) = round(apCentre,-1);
        end
        
        data = img(handles.crop_rect(2):(handles.crop_rect(2)+handles.crop_rect(4)-1), handles.crop_rect(1):(handles.crop_rect(1)+handles.crop_rect(3)-1));

        % some default parameters
        range = 15;
        
        [m,n] = size(data);
        %load('parms.mat') %not necessary
        %automatic load parms independently of the current location to
        %prevent errors
        tmp = mfilename; 
        tmp = erase(which(tmp),[tmp '.m']);
        load([tmp 'updates/parms.mat'])
        
        parms.apo.apox = round(linspace(1, size(data,2), 10));
        parms.apo.super.cut = [max(apRound(1)-range, 0), apRound(1)+range] / 100;
        parms.apo.deep.cut = [apRound(2)-range, min(apRound(2)+range, 100)] / 100;
        
        % save in handles
        handles.Hough.apo.deep.cut = parms.apo.deep.cut;
        handles.Hough.apo.super.cut = parms.apo.super.cut;
        
        % filter and fit aponeurosis
        parms.fas.filter = 0;
        [~, super_obj, deep_obj] = filter_usimage(data,parms);
        
        super_coef = polyfit(parms.apo.apox,super_obj, 1);
        deep_coef = polyfit(parms.apo.apox,deep_obj, 1);
        
        handles.super_coef{frame_no} = super_coef;
        handles.deep_coef{frame_no} = deep_coef;
        
        h1 = drawline('Color','g' , 'InteractionsAllowed','all', 'Position', [1 n; polyval(super_coef, 1) polyval(super_coef, n)]');
        h2 = drawline('Color','g' , 'InteractionsAllowed','all', 'Position', [1 n; polyval(deep_coef, 1) polyval(deep_coef, n)]');
        
        apo1.line = round(h1.Position);
        apo2.line = round(h2.Position);
        
        if ii == 1
            i = 1;
            j = 1;
        elseif ii == 2
            i = 2;
            j = 1;
        end
        
        [handles.Region(i).ROI{frame_no},handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no}] = roipoly(data,...
        vertcat(apo1(1).line(1,1),apo2(1).line(1,1),apo2(1).line(2,1),apo1(1).line(2,1),apo1(1).line(1,1)),...
        vertcat(apo1(1).line(1,2),apo2(1).line(1,2),apo2(1).line(2,2),apo1(1).line(2,2),apo1(1).line(1,2)));
        
        %% fascicle ROI
        % construct ellipse region
        a = handles.Region.ROIy{frame_no}(2) - handles.Region.ROIy{frame_no}(1);
        b = handles.Region.ROIx{frame_no}(3) - handles.Region.ROIx{frame_no}(2);
        
        s  = regionprops(handles.Region.ROI{frame_no}, 'centroid', 'MinorAxisLength','MajorAxisLength');

        rhorz = n;
        rvert = polyval(deep_coef, n/2) - polyval(super_coef, n/2);
        
%         rvert = s.MinorAxisLength;
%         rhorz = s.MajorAxisLength;
        
%         [nx,ny] = size(data);
%         a = nx/2;
%         b = ny/2;
        
        % vertical radius of ellipse
        gamma = -atan2d((handles.Region.ROIy{frame_no}(3) - handles.Region.ROIy{frame_no}(2)), b);
        
        % points on ellipse relative to centroid
        th = linspace(0,2*pi);
        xc = rhorz/2 * cos(th); 
        yc = rvert/2 * sin(th); 
        
        % rotate relative points and add centroid
        xcR = xc * cosd(-gamma) - yc * sind(-gamma) + s.Centroid(1);
        ycR = xc * sind(-gamma) + yc * cosd(-gamma) + s.Centroid(2);
              
        [nx,ny] = size(data);
        [X,Y] = meshgrid(1:ny,1:nx) ;
        
        idx = inpolygon(X, Y, round(xcR), round(ycR));
         
        handles.Region.fasROI = ones(size(data));
        handles.Region.fasROI(~idx) = 0;
        
        handles.Region.fasROI_radius = [rvert rhorz];
        handles.Region.fasROI_angle = gamma;
        
        % first call to determine quadrant
        handles.Hough.range = [1 180];
        handles.Hough.npeaks = 1;
        
        handles = do_TimTrack(handles, frame_no,i,j);
        
        if handles.Hough.Hough_angle{frame_no} < 90 && handles.Hough.Hough_angle{frame_no} > 0
            handles.Hough.range = [5 80];
            handles.flip = 0;
        else
            handles.Hough.range = [100 175];
            handles.flip = 1;
        end
        
        % second call to estimate first frame
        handles.Hough.npeaks = 10;
        handles = do_TimTrack(handles, frame_no,i,j);
        
        for k = 1:2
            if frame_no == 1
                handles.Region.Fascicle.fas_x_original{frame_no}(k) =  handles.Region.Fascicle.fas_x{frame_no}(k);
                handles.Region.Fascicle.fas_y_original{frame_no}(k) =  handles.Region.Fascicle.fas_y{frame_no}(k);
            end
        end
        
        %% whole fascicle, whole roi
        if ii == 1
            i = 1;
            j = 1;
%             [handles.Region(i).ROI{frame_no},handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no}] = roipoly(handles.IM,...
%             vertcat(apo1(1).line(1,1),fasROI.position(4,1),fasROI.position(3,1),apo1(1).line(2,1),apo1(1).line(1,1)),...
%             vertcat(apo1(1).line(1,2),fasROI.position(4,2),fasROI.position(3,2),apo1(1).line(2,2),apo1(1).line(1,2)));

        elseif ii == 2
            i = 2;
            j = 1;
%             [handles.Region(i).ROI{frame_no},handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no}] = roipoly(handles.IM,...
%             vertcat(fasROI.position(1,1),apo2(1).line(1,1),apo2(1).line(2,1),fasROI.position(2,1),fasROI.position(1,1)),...
%             vertcat(fasROI.position(1,2),apo2(1).line(1,2),apo2(1).line(2,2),fasROI.position(2,2),fasROI.position(1,2)));
        end

%         [handles.Region(i).ROI{frame_no},handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no}] = roipoly(
        

%         [handles.Region(i).ROI{frame_no},handles.Region(i).ROIx{frame_no}, handles.Region(i).ROIy{frame_no}] = roipoly(handles.IM,...
%             vertcat(fasROI.position(1,1),fasROI.position(4,1),fasROI.position(3,1),fasROI.position(2,1),fasROI.position(1,1)),...
%             vertcat(fasROI.position(1,2),fasROI.position(4,2),fasROI.position(3,2),fasROI.position(2,2),fasROI.position(1,2)));
%         handles.Region(i).Fascicle(j).fas_x{frame_no}(1) = FPT(1).apo1int(1,1);
%         handles.Region(i).Fascicle(j).fas_x{frame_no}(2) = FPT(1).apo2int(1,1);
%         handles.Region(i).Fascicle(j).fas_y{frame_no}(1) = FPT(1).apo1int(1,2);
%         handles.Region(i).Fascicle(j).fas_y{frame_no}(2) = FPT(1).apo2int(1,2);

        if handles.Region(i).Fascicle(j).fas_y{frame_no}(1) < handles.Region(i).Fascicle(j).fas_y{frame_no}(2)
            handles.Region(i).Fascicle(j).fas_y{frame_no} = fliplr(handles.Region(i).Fascicle(j).fas_y{frame_no});
            handles.Region(i).Fascicle(j).fas_x{frame_no} = fliplr(handles.Region(i).Fascicle(j).fas_x{frame_no});
        end

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

    end

    clear img
    show_image(hObject,handles)

end



function TimTrack_downsample_Callback(hObject, eventdata, handles)
% hObject    handle to TimTrack_downsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimTrack_downsample as text
%        str2double(get(hObject,'String')) returns contents of TimTrack_downsample as a double



% --- Executes during object creation, after setting all properties.
function TimTrack_downsample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimTrack_downsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TTgain_Callback(hObject, eventdata, handles)
% hObject    handle to TTgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TTgain as text
%        str2double(get(hObject,'String')) returns contents of TTgain as a double


% --- Executes during object creation, after setting all properties.
function TTgain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TTgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_TimTrack.
function run_TimTrack_Callback(hObject, eventdata, handles)
% hObject    handle to run_TimTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if not defined, assume all frames
if strcmp(handles.TimTrack_downsample.String,'0') || isempty(handles.TimTrack_downsample.String)
    handles.Hough_frames = 1:handles.NumFrames;
else 
    handles.Hough_frames = 1:str2double(handles.TimTrack_downsample.String):handles.NumFrames;
end

if isfield(handles,'ImStack')
    % find current frame number from slider
    frame_no = round(get(handles.frame_slider,'Value'));

    if ~isfield(handles.Region(1).Fascicle(1),'current_xy')
        warndlg('First frame must have fascicle and ROI defined', 'Warning');
    else

        h = waitbar(0,'Please wait while processing forward...');

        for i = 1:length(handles.Region)

            % create image from movie frame
            im1 = imcrop(handles.ImStack(:,:,frame_no+handles.start_frame-1),handles.crop_rect);
            handles.CIm = im1;

            for f = handles.Hough_frames
                % Get the current image

                im2 = imcrop(handles.ImStack(:,:,handles.start_frame+f-1),handles.crop_rect);
                handles.NIm = im2;
                handles.CIm = handles.NIm;

                for j = 1:length(handles.Region(i).Fascicle)
%                     profile on
                    handles = do_TimTrack(handles, f,i,j);
%                     profile viewer
                end

                waitbar((f+(get(handles.frame_slider,'Max')*(i-1)))/...
                    (get(handles.frame_slider,'Max')*length(handles.Region)),h)
            end
        end
        close(h)
    end

    set(handles.frame_slider,'Value',frame_no);
    set(handles.frame_number,'String',num2str(frame_no));

    % update the image axes using show_image function (bottom)
    show_image(hObject,handles)
end

% Update handles structure
guidata(hObject, handles);



function Cut_off_time_Callback(hObject, eventdata, handles)
% hObject    handle to Cut_off_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cut_off_time as text
%        str2double(get(hObject,'String')) returns contents of Cut_off_time as a double


% --- Executes during object creation, after setting all properties.
function Cut_off_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cut_off_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ApoDetect_downsample_Callback(hObject, eventdata, handles)
% hObject    handle to ApoDetect_downsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ApoDetect_downsample as text
%        str2double(get(hObject,'String')) returns contents of ApoDetect_downsample as a double


% --- Executes during object creation, after setting all properties.
function ApoDetect_downsample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ApoDetect_downsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ApoGain_Callback(hObject, eventdata, handles)
% hObject    handle to ApoGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ApoGain as text
%        str2double(get(hObject,'String')) returns contents of ApoGain as a double


% --- Executes during object creation, after setting all properties.
function ApoGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ApoGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clear_correction.
function clear_correction_Callback(hObject, eventdata, handles)
% hObject    handle to clear_correction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


for j = 1:2 % superficial and deep aponeurosis
    for i = 1:handles.NumFrames
        handles.Region.Fascicle.fas_x{i}(j) = handles.Region.Fascicle.fas_x_original{i}(j);
        handles.Region.Fascicle.fas_y{i}(j) = handles.Region.Fascicle.fas_y_original{i}(j);
    end
end

for frame_no = 1:handles.NumFrames
    handles.Region(1).fas_length(frame_no,1) = handles.ID *sqrt(diff(handles.Region(1).Fascicle(1).fas_y{frame_no}).^2 + ...
                                                                diff(handles.Region(1).Fascicle(1).fas_x{frame_no}).^2);
end

show_image(hObject,handles)
        


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1

function[parms] = get_TimTrack_parms()
    % TimTrack parameters
    load('parms.mat','parms')
    % parms.cut_image = 0;

    parms.ROI = [239   936; 50   608]; % [0812]
    parms.extrapolation = 0;
    parms.fas.thetares = .5;

    parms.apo.deep.cut(1) = .35;
    parms.apo.super.cut(1) = .03;
    
    parms.show = 0;
    parms.fas.show = 0;
    parms.redo_ROI = 0;
