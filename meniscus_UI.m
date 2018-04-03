function varargout = meniscus_UI(varargin)
% MENISCUS_UI MATLAB code for meniscus_UI.fig
%      This program provides a graphical user interface to track the 
%      movement of a meniscus in a capillary.
%      
%      Contents:
%      1. Default code: For initialising the GUI
%      2. Get rotation function: Determines angle to rotate by given a line
%      3. Image tilt function: Rotates images & fills in empty spots
%      4. Low pass filter: Blurs the image to help with meniscus location
%      5. File select button: Lets user select files for analysis and runs
%         code to tilt and trim the images
%      6. Size calibration button: User selects ruler image for calibration
%      7. Rotation calibration button: Lets user specify rotation
%      8. Background removal button: Averages images and removes background
%      9. Trim button: Allows user to discard irrelevant parts of the image
%      10. Low pass button: Passes images through low pass filter
%      11. Auto detect button: Uses 'sobel type' detector to track meniscus
%      12. Manual detect button: Allows user to click along meniscus path
%      13. Plot button: Plots penetration of meniscus vs time
%      14. Radio buttons: For selecting original, bg_sub and lp images
%      15. Display buttons: Displays the selected set of images
%      16. Initial penetration button: Calculates penetration of 1st datum
%      17. Save button: Saves time, distance and parameters into txt file
%      18. Default code: Handles objects in UI
%      
%%% DEFAULT STUFF %%%
%      MENISCUS_UI, by itself, creates a new MENISCUS_UI or raises the existing
%      singleton*.
%
%      H = MENISCUS_UI returns the handle to a new MENISCUS_UI or the handle to
%      the existing singleton*.
%
%      MENISCUS_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MENISCUS_UI.M with the given input arguments.
%
%      MENISCUS_UI('Property','Value',...) creates a new MENISCUS_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before meniscus_UI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to meniscus_UI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help meniscus_UI
% Last Modified by GUIDE v2.5 23-Nov-2017 11:43:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @meniscus_UI_OpeningFcn, ...
                   'gui_OutputFcn',  @meniscus_UI_OutputFcn, ...
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

% --- Executes just before meniscus_UI is made visible.
function meniscus_UI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to meniscus_UI (see VARARGIN)

% Choose default command line output for meniscus_UI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes meniscus_UI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = meniscus_UI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%% GET ROTATION FUNCTION %%%
function [tta,crds] = getRotation(pics, hObject, eventdata, handles)
% This function determines the difference between a selected image from
% horizontal and returns tta (the angle) and crds (the coordinates of the
% points selected by the user)

image(pics{1}{1}); colormap gray; % Plot the first of the selected images

but=0; % button is a vector of integers indicating which mouse buttons
       % you pressed (1 for left, 2 for middle, 3 for right), from ginput()

ct = 1;     % Initialise click count and coordinates
crds = []; 

% This section lets the user select desired points and only stops when
% they're done
while size(crds,2)<2 && ~strcmp(crds, 'bypass') 
    % Double operators act as a short-circuit (don't test 2nd if 1st fails)
    % strcmp compares 2 strings = 1/0 if matching/not

    [a,b,but] = ginput(1); % Graphical input
    
    % Left click - Add click location to crds
    if but == 1 && size(crds,2)<2 
        crds(1,ct) = a; 
        crds(2,ct) = b; 
        ct = ct+1;
    end
    
    % Middle click - Bypass rotation step
    if but == 2 && size(crds,2)==0; 
        crds='bypass'; % bypass rotation
    
    % Right click - Undo
    elseif but == 3
        cla % Clear and reset the image
        image(pics{1}{1}); colormap gray; 
        crds = []; 
        ct = 1;
    end

    % Plot the user selected line
    hold on; 
    try
        plot(crds(1,:),crds(2,:),'-g') 
        plot(crds(1,:),crds(2,:),'og')
    catch ploterror; 
    end
    hold off;
end

% Caluclate rotation angle
pic = pics{1}{1};
s = size(pic);
[x y] = meshgrid(1:s(2),1:s(1)); % Define original coordinates
A = [crds(1,2)-crds(1,1); crds(2,2)-crds(2,1)]; % Define vector of new axis
Alen = sqrt( A(1).^2 + A(2).^2 ); 
A = A/Alen; % Make A a unit vector 
B = [0 ; -1]; % Define the old axis
tta = asin(dot(B,A)); % Define the angle between


%%% IMAGE TILT FUNCTION %%%
function [newpic] = my_imtilt(pic,tta,type)
% This program rotates an image through a given angle, and automatically 
% fills in blank points to make it the same shape.

intp = 'y'; %Pick 'y' for nearest neighbour, 
                 %'n' for none, 
                 %'a' for an average fill in blank points.

s = size(pic);
[x y] = meshgrid(1:s(2),1:s(1)); % Define original coordinates

x2 = x*cos(tta) - y*sin(tta); % Rotate xs and ys in the usual way. 
y2 = x*sin(tta) + y*cos(tta); 

x2 = round(x2); % Round so that they fall on pixels
y2 = round(y2);

% Create 's2' which defines the size of the square the rotated will fit
% onto, this will be called newpic
miny2 = min(min(y2)); 
maxy2 = max(max(y2)); 
minx2 = min(min(x2)); 
maxx2 = max(max(x2)); 
xshift = 0; 
yshift = 0; 
if miny2<0 
    yshift = -1*miny2;
end
if minx2<0 
    xshift = -1*minx2;
end
s2 = [(maxy2-miny2+2) (maxx2-minx2+2) size(pic,3)]; 

% Deals with the spaces left between pixels in the rotated image. 
% This is done by a cunning choice of background colour: 
% - Define the new background, newpic, to be zeros then gaps are all black. ('n')
% - Define the newpic to the average colour of pic then the gaps may blend in more. ('a') 
% - If we give the background a distinctive value - eg. pi - 
%   we can later fill the cells of value pi with an average value from the adjacent cells. ('y') 
if intp == 'y'
    newpic = ones(s2)*pi;
elseif intp == 'a'
    newpic = ones(s2)*mean(mean(pic));
else
    newpic = zeros(s2);
end

% This loop fills the pixels from pic into newpic according to their new x and y values. (x2 and y2) 
for m = 1:s(1) 
    for n = 1:s(2)
        newpic(y2(m,n)+yshift+1,x2(m,n)+xshift+1,:) = pic(m,n,:);
    end
end

s = size(newpic);

% This section is the 'filler-inner'. 
if intp == 'y' 
    for loop=1:size(pic,3) 
        for m = 2:s(1)-1 
            for n = 2:s(2)-1 
                if newpic(m,n,loop) == pi;
                    a1 = newpic(m+1,n,loop);
                    a2 = newpic(m-1,n,loop);
                    a3 = newpic(m,n+1,loop);
                    a4 = newpic(m,n-1,loop);
                    c1 = [newpic(m+1,n,loop),newpic(m-1,n,loop),newpic(m,n+1,loop),newpic(m,n-1,loop)]; 
                    c2 = 0; 
                    l = 0; 
                    for ct2 = 1:4 
                        if c1(ct2) ~= pi 
                            c2=c2+c1(ct2); 
                            l = l+1;
                        end
                    end
                    mc2 = c2/l; 
                    if mc2>0
                        newpic(m,n,loop) = mc2-pi;
                    end
                end
            end
        end
    end
end
newpic(newpic==pi) = 0; % Set the outsides back to 0;


%%% LOW PASS FILTER %%%
function [low_passed_pic] = lp_filtery(pic,ratio)
% This function will return the image after convolution with a box plot
% This blurs the image such that non-meniscus things should be hidden
boxKernel = ones(10,6); % Or whatever size window you want.
blurredImage = conv2(pic, boxKernel, 'same');
low_passed_pic = blurredImage;


%%% FILE SELECT BUTTON %%%
function file_selector_Callback(hObject, eventdata, handles)
% This causes the file browser to open and allows user to choose files
% All files are rotated & trimmed according to current values
clear global
% Initialise some things for later use
axes(handles.axes1);
view_time = str2double(get(handles.view_time, 'string'));
instructions = 'Loading...';
set(handles.instruction_box, 'String', instructions);

[FileName,PathName,FilterIndex] = uigetfile(...
    'C:\Users\George\OneDrive\University\Chemistry\Year 4 Chemistry\Project\Drive data\Meniscus experiments\*.tif', ...
    'MultiSelect','on');
%[FileName,PathName,FilterIndex] = uigetfile(...
%    'C:\Users\George\Desktop\*.tif', ...
%    'MultiSelect','on');

% Filename = 1xn array for the n files that you select, eg {'file1.tif'}...
% PathName = Direction to the folder you're in, eg C:\Users\George\Desktop
global file_folder
file_folder = PathName; % To include in output save file

% Create 'pics' which will be a 1xn cell array:
% - pics{i}{1} = image data for the ith image (huge array of 1 & 0)
% - pics{i}{2} = [image height, image width] (ie [no.rows, no. cols]
top_trim = str2double(get(handles.top_trim, 'string'));
bottom_trim = str2double(get(handles.bottom_trim, 'string'));
left_trim = str2double(get(handles.left_trim, 'string'));
right_trim = str2double(get(handles.right_trim, 'string'));
if iscell(FileName) == 0 % iscell asks if it's a cell array, 0 means no
    % This is for when there is one input and you have a character vector
    file_path = strcat(PathName, FileName);
    pic = imread(file_path);
    picsize = size(pic);
    pic = imcrop(pic,[left_trim, top_trim, right_trim-left_trim, bottom_trim-top_trim]);
    % Imcrop takes (pic, [xmin, ymin, width, height])
    picsize = size(pic);
    pics{1} = {pic, picsize};
else
    % For more than 1 file input
    for i =1:length(FileName)
        file_path = strcat(PathName, cell2mat(FileName(i)));
        % strcat concatatanates 2 strings and cell2mat converts to matrix
        pic = imread(file_path);
        picsize = size(pic);
        pic = imcrop(pic,[left_trim, top_trim, right_trim-left_trim, bottom_trim-top_trim]);
        picsize = size(pic);
        pics{i} = {pic, picsize};
    end
end

% Get the necessary rotation values
tta = str2double(get(handles.rot_angle, 'string')) * 2*pi / 360; 
if tta ~= 0
    for i =1:length(pics)
        pics{i}{1} = my_imtilt(pics{i}{1},tta,'pic'); % do rotations 
    end
else
    disp('Rotation function bypassed');
end

global original_pics;
original_pics = pics;

global current_pics;
clear current_pics;
current_pics = pics; % Put pics into a global variable for later use

set(handles.org_sel, 'Value', 1);
set(handles.bgsub_sel, 'Value', 0);
set(handles.lp_sel, 'Value', 0);

% Display selected images
for i = 1:length(current_pics)
    imagesc(current_pics{i}{1}); colormap gray;
    daspect([1 1 1]); %Stops the image from stretching 
    pause(view_time);
end

set(handles.instruction_box, 'String', ' ');

%%% SIZE CALIBRATION BUTTON %%%
function calibrate_size_button_Callback(hObject, eventdata, handles)
% User can select a calibration image (ie pic of ruler). They then calibrate
% image size by drawing a line over the ruler markings. 
% The number of markings are found using the number of minima in an intensity
% profile of the line drawn over the ruler.

 % User selects ruler file
[FileName,PathName,FilterIndex] = uigetfile(...
    'C:\Users\George\OneDrive\University\Chemistry\Year 4 Chemistry\Project\Drive data\Meniscus experiments\*.tif');
file_path = strcat(PathName, FileName);
pic = imread(file_path);
picsize = size(pic);
pics{1} = {pic, picsize};

% Rotate the image by currently selected amount
tta = str2double(get(handles.rot_angle, 'string')) * 2*pi / 360;
if tta == 0
    disp('Rotation function bypassed');
else
    pic = my_imtilt(pics{1}{1},tta,'pic'); % do rotation 
    imagesc(pic)
    pics{1} = {pic, picsize};
end

% User draws line
instructions = 'Draw a roughly vertical line along the markings of the ruler';
set(handles.instruction_box, 'String', instructions);
[tta, ruler_crds] = getRotation(pics);
set(handles.instruction_box, 'String', ' ');

xline = [ruler_crds(1,1),ruler_crds(1,1)]; % Forces selected line vertical
yline = ruler_crds(2,:);

% Get intensity profile along line
intensities = improfile(pic, xline, yline); 
avg_int = mean(intensities);

% Determine the number of minima in the intensity profile
no_markings = 0;
quartiles = quantile(intensities,[0.25, 0.5, 0.75]);
for i = 2:length(intensities)
    try
        prev_int = intensities(i-1);
        int = intensities(i);
        next_int = intensities(i+1);
        
        % Check if point is local min and also has comparatively low intensity
        if int <= prev_int && int <= next_int & int<quartiles(3)- 10
            % Need equality since some dips are 'flat'
            
            if no_markings == 0 % For the first min
                first_mark = i;
                i_hit = i;
                no_markings = no_markings + 1;
                
            elseif i > (i_hit + 10) % Avoids ruler markings with a 'double dip'
                no_markings = no_markings + 1;
                last_mark = i;
                i_hit = i;
            end
        end
    end
end

mm_per_px = (no_markings - 1) /(last_mark - first_mark);
%pic_width = pics{1}{2}(2) * mm_per_px
%pic_height = pics{1}{2}(1) * mm_per_px
set(handles.mm_per_px, 'String', mm_per_px);

%%% ROTATION CALIBRATION BUTTON %%%
function calibrate_rot_button_Callback(hObject, eventdata, handles)
% Allows user to select a file and then define a horizontal axis by drawing
% a line, subsequent images used will be rotated by this amount.

% User selects ruler file
[FileName,PathName,FilterIndex] = uigetfile(...
    'C:\Users\George\OneDrive\University\Chemistry\Year 4 Chemistry\Project\Drive data\Meniscus experiments\*.tif'); 

top_trim = str2double(get(handles.top_trim, 'string'));
bottom_trim = str2double(get(handles.bottom_trim, 'string'));
left_trim = str2double(get(handles.left_trim, 'string'));
right_trim = str2double(get(handles.right_trim, 'string'));
axes(handles.axes1); %set the current axes to axes1
file_path = strcat(PathName, FileName);
pic = imread(file_path);
picsize = size(pic);
pic = imcrop(pic,[left_trim, top_trim, right_trim-left_trim, bottom_trim-top_trim]);
picsize = size(pic);
pics{1} = {pic, picsize};

% Get the necessary rotation values
instructions = 'Left click 2 points to define horizontal axis, Right click to undo, Middle click to bypass.';
set(handles.instruction_box, 'String', instructions);
[tta, crds] = getRotation(pics);
set(handles.rot_angle, 'String', tta*360/(2*pi));
set(handles.instruction_box, 'String', ' ');

% Rotate the image
if tta == 0
    disp('Rotation function bypassed');
else
    pic = my_imtilt(pics{1}{1},tta,'pic'); % do rotation 
    imagesc(pic)
    pics{1} = {pic, picsize};
end


%%% BACKGROUND REMOVAL BUTTON %%%
function bg_remove_button_Callback(hObject, eventdata, handles)
% Subtract the average of all images from each image to reduce noise

% Initialise things
axes(handles.axes1); %set the current axes to axes1
view_time = str2double(get(handles.view_time, 'string'));
global original_pics;
global bgsub_pics;
bgsub_pics = original_pics;

%{
% 1ST PIC REMOVAL
% Take 1st pic away from all images (except 1st)
for i = 2:length(original_pics)
    bgsub_pics{i}{1}= double(original_pics{i}{1}) - double(original_pics{1}{1});
end
%}

%{
% CURRENT AVERAGE
% Removes 'current' average (later images not included in background)
pic_sum = double(zeros(size(original_pics{1}{1})));
for i = 1:length(original_pics)
    pic_sum = pic_sum + double(original_pics{i}{1});
    bgsub = pic_sum/i;
    bgsub_pics{i}{1} = double(original_pics{i}{1}) - bgsub;
end
%}    


% SUBTRACT PREVIOUS IMAGE
% Removes 'current' average (later images not included in background)
pic_sum = double(zeros(size(original_pics{1}{1})));
bgsub_pics{1}{1} = double(original_pics{1}{1} - original_pics{2}{1});
for i = 2:length(original_pics)
    bgsub_pics{i}{1} = double(original_pics{i}{1} - original_pics{i-1}{1});
end


%{
% FULL AVERAGE REMOVAL
% Sum all the image values & divide by no. for the average image
pic_sum = double(zeros(size(original_pics{1}{1})));
for i = 1:length(original_pics)
    pic_sum = pic_sum + double(original_pics{i}{1});
end
bgsub = pic_sum/length(original_pics);

% Take background away from all images
for i = 1:length(original_pics)
    bgsub_pics{i}{1}= double(original_pics{i}{1}) - bgsub;
end
%}

% Display the background-subtracted images
for i = 1:length(bgsub_pics)
    image(bgsub_pics{i}{1});
    daspect([1 1 1]) %Stops the image from stretching 
    pause(view_time)
end

global current_pics;
current_pics = bgsub_pics;
set(handles.org_sel, 'Value', 0);
set(handles.bgsub_sel, 'Value', 1);
set(handles.lp_sel, 'Value', 0);


%%% TRIM BUTTON %%%
function trim_button_Callback(hObject, eventdata, handles)
% Lets user discard the top & bottom of the image to reduce computation

% Let user select file and display the image
axes(handles.axes1);
[FileName,PathName,FilterIndex] = uigetfile(...
    'C:\Users\George\OneDrive\University\Chemistry\Year 4 Chemistry\Project\Drive data\Meniscus experiments\*.tif'); 
file_path = strcat(PathName, FileName);
pic = imread(file_path);
imagesc(pic); colormap gray;

% Get the user input for top, bottom, left and right boundaries
set(handles.instruction_box, 'String', 'Click an upper boundary');
[~,yu] = ginput(1);
yu = round(yu);
set(handles.top_trim, 'String', yu);

set(handles.instruction_box, 'String', 'Click a lower boundary');
[~,yl] = ginput(1);
yl = round(yl);
set(handles.bottom_trim, 'String', yl);

set(handles.instruction_box, 'String', 'Click a left boundary');
[xl,~] = ginput(1);
xl = round(xl);
set(handles.left_trim, 'String', xl);

set(handles.instruction_box, 'String', 'Click a right boundary');
[xr,~] = ginput(1);
xr = round(xr);
set(handles.right_trim, 'String', xr);

% Plot trimmed image
picsize = size(pic);
pic = imcrop(pic,[xl, yu, xr-xl, yl-yu]);
imagesc(pic); colormap gray;


%%% LOW PASS FILTER BUTTON %%%
function lp_filter_Callback(hObject, eventdata, handles)
% Calls the low pass filter function to run on current images

% Initialise things
axes(handles.axes1); 
global bgsub_pics;
global lp_pics;
lp_pics = bgsub_pics;
view_time = str2double(get(handles.view_time, 'string'));

% Define lp_pics to be low pass filtered images
for i = 1:length(bgsub_pics)
    pic = bgsub_pics{i}{1};
    lp_pics{i}{1}= lp_filtery(pic,0.01);
end

% Plot low pass filtered images
for i = 1:length(lp_pics)
    imagesc(lp_pics{i}{1}); colormap gray;
    daspect([1 1 1]); %Stops the image from stretching 
    pause(view_time);
end

global current_pics;
current_pics = lp_pics;
set(handles.org_sel, 'Value', 0);
set(handles.bgsub_sel, 'Value', 0);
set(handles.lp_sel, 'Value', 1);

%%% AUTO DETECT BUTTON %%%
function auto_detect_button_Callback(hObject, eventdata, handles)
% Function tries to detect the meniscus location automatically

% Initialise things
axes(handles.axes1); 
global current_pics;
global lp_pics;
pics = lp_pics;
view_time = str2double(get(handles.view_time, 'string'));

global men_xs_calc  % This will be the calculated positions

sobelx1 = [-1 0 1; -2 0 +2; -1 0 1];  % Horizontal sobel operator
modified = zeros(1,length(men_xs_calc)); % Will allow to remove bad points
for i = 1:length(lp_pics)
    pic = pics{i}{1};
    totals=sum(pic); % Produces row with sum of each column
    sumedge = conv2(totals,[-1 0 1],'same'); % Convolves with simplified sobel
    sumedge=abs(sumedge); % Interested in total gradient, not direction
    
    % Locate meniscus for current image
    % If loop ensures meniscus is not selected behind current location
    if i ~= 1 % Irrelevant for first image
        for x = 1:d2
            if sumedge(x) > d2 && x < d2
                sumedge(x) = 0;               
            end
        end
    end
    
    [~,d2]=max(sumedge);
    
    % Catches any values that are still somehow lower than the previous
    if i ~= 1
        if d2 < men_xs_calc(i-1)
            d2 = men_xs_calc(i-1);
            modified(i) = 1;
        end
    end
    
    % If autodetect jumps more than set # of px, then leave this point out
    if i ~= 1
        if d2 > men_xs_calc(i-1) + 13
            d2 = men_xs_calc(i-1);
            modified(i) = 1;
        end
    end
    
    if i == 1 
        if d2 > 20
            d2 = 1
        end
    end
    
    men_xs_calc(i) = d2;
    
    % Plot what has happened
    %{
    % Plot what a real sobel operator sees
    edge = conv2(pic,sobelx1,'same');
    %BW = edge(pic,'sobel','vertical'); % Plot the edges found
    hold on
    imagesc(edge); colormap gray;
    %}
    hold on
    imagesc(current_pics{i}{1}); colormap gray;
    plot(d2,20,'r*')
    hold off
    pause(view_time);
end

% Remove data points that were artificially set equal to previous value
for i = 0:length(men_xs_calc)-2
    dist = men_xs_calc(end-i);
    prev_dist = men_xs_calc(end-i-1);
    if dist == prev_dist && modified(end-i) == 1
        men_xs_calc(end-i) = NaN;
    end
end


%%% MANUAL DETECT BUTTON %%%
function manual_detect_button_Callback(hObject, eventdata, handles)
% Allows user to manually choose the meniscus location

% Initialise things
axes(handles.axes1); %set the current axes to axes1
global current_pics;
pics = current_pics;
global men_xs_manual
%clear men_xs_manual

set(handles.instruction_box, 'String', 'Click on the meniscus');

% Get user input for each image
for i = 1:length(current_pics)
    pic = pics{i}{1};
    imagesc(pic); colormap gray;
    daspect([1 1 1]) %Stops the image from stretching 
    [xpic,~] = ginput(1); 
    xpic=round(xpic);
    men_xs_manual(i)=xpic;
end
set(handles.instruction_box, 'String', 'Meniscus locations stored');

%%% PLOT BUTTON %%%
function plot_button_Callback(hObject, eventdata, handles)
% Plots penetration vs time (or t^0.5) for current data

% Initialise things
global current_pics;
pics = current_pics;
global men_xs_calc;
global men_xs_manual;

mm_per_px = str2double(get(handles.mm_per_px, 'string'));
axes(handles.axes2); %set the current axes to axes2
cla

% Calculate the timestep
no_pics = length(pics);
fps = 25;
frametime=1/fps; % get time length of each frame 
timestep=0:frametime:(frametime*(no_pics-1));

i_frames = str2double(get(handles.i_frames, 'string'));
i_time = i_frames*frametime;
timestep = timestep + i_time;

% Convert the meniscus distances into mm from px
i_pen = str2double(get(handles.i_penetration, 'string'));
dist_conv_calc = men_xs_calc * mm_per_px + i_pen;
dist_conv_manual = men_xs_manual * mm_per_px + i_pen;

% Plot desired graph
hold on
sqrt_check = get(handles.sqrt_checkbox, 'Value');
if sqrt_check == 1  % t^0.5 selected
    try
        plot(sqrt(timestep),dist_conv_calc, 'b--o')
    end
    
    try
        plot(sqrt(timestep),dist_conv_manual, 'r--o')
    end
    
    xlabel('t^{1/2} / s^{1/2}') % x-axis label
else                % t selected
    try
        plot(timestep,dist_conv_calc, 'b--o')
    end
    
    try
        plot(timestep,dist_conv_manual, 'r--o')
    end
    
    xlabel('Time / s') % x-axis label
end
%xlim([timestep(1), timestep(end)]) % Define x axis limits
hold off
ylabel('Penetration / mm') % y-axis label

% Add legend if both auto and manual are being used
if ~isempty(dist_conv_calc) && ~isempty(dist_conv_manual)
    legend('Auto-detection','Manual-detection','Location','east')
end

% Put data into global variables for saving
global calc_dist
calc_dist = dist_conv_calc - i_pen;
global man_dist
man_dist = dist_conv_manual - i_pen;
global times
times = timestep - i_time;


%%% ORIGINAL IMAGE RADIO BUTTON %%%.
function org_sel_Callback(hObject, eventdata, handles)
% Changes selection to original images
global current_pics;
global original_pics;
current_pics = original_pics;
set(handles.org_sel, 'Value', 1);
set(handles.bgsub_sel, 'Value', 0);
set(handles.lp_sel, 'Value', 0);


%%% BACKGROUND SUBTACTED RADIO BUTTON %%%
function bgsub_sel_Callback(hObject, eventdata, handles)
% Changes selection to background-subtracted images
global current_pics;
global bgsub_pics;
current_pics = bgsub_pics;
set(handles.org_sel, 'Value', 0);
set(handles.bgsub_sel, 'Value', 1);
set(handles.lp_sel, 'Value', 0);


%%% LOW-PASS RADIO BUTTON %%%
function lp_sel_Callback(hObject, eventdata, handles)
% Changes selection to low-pass filtered images
global current_pics;
global lp_pics;
current_pics = lp_pics;
set(handles.org_sel, 'Value', 0);
set(handles.bgsub_sel, 'Value', 0);
set(handles.lp_sel, 'Value', 1);


%%% DISPLAY BUTTON %%%
function disp_button_Callback(hObject, eventdata, handles)
% Plays through the currently selected images

axes(handles.axes1);
gca;
global current_pics;
view_time = str2double(get(handles.view_time, 'string'));

for i = 1:length(current_pics)
    imagesc(current_pics{i}{1}); colormap gray;
    daspect([1 1 1]); %Stops the image from stretching 
    pause(view_time);
end

%%% INITIAL PENETRATION BUTTON %%%
function i_penetration_button_Callback(hObject, eventdata, handles)
% Calculates the initial penetration from when data collecting begins,
% using a user drawn line

% Open image , trim and rotate it
[FileName,PathName,FilterIndex] = uigetfile(...
    'R:\Lab_Data\George Lewis\Meniscus experiments\*.tif'); 

% Trim
top_trim = str2double(get(handles.top_trim, 'string'));
bottom_trim = str2double(get(handles.bottom_trim, 'string'));
axes(handles.axes1); %set the current axes to axes1
file_path = strcat(PathName, FileName);
pic = imread(file_path);
picsize = size(pic);
pic = imcrop(pic,[0, top_trim, picsize(2), bottom_trim-top_trim]);
picsize = size(pic);
pics{1} = {pic, picsize};

% Rotate
tta = str2double(get(handles.rot_angle, 'string')) * 2*pi / 360; 
if tta ~= 0
    for i =1:length(pics)
        pics{i}{1} = my_imtilt(pics{i}{1},tta,'pic'); % do rotations 
    end
else
    disp('Rotation function bypassed');
end

% Display trimmed and rotated image
imagesc(pics{1}{1}); colormap gray;
daspect([1 1 1]); %Stops the image from stretching 

% Get user to draw line
instructions = 'Draw line from left capillary end to first data point';
set(handles.instruction_box, 'String', instructions);


but=0; % button is a vector of integers indicating which mouse buttons
       % you pressed (1 for left, 2 for middle, 3 for right), from ginput()

ct = 1;     % Initialise click count and coordinates
crds = []; 

% This section lets the user select desired points and only stops when
% they're done
while size(crds,2)<2 && ~strcmp(crds, 'bypass') 
    [a,b,but] = ginput(1); % Graphical input
    
    % Left click - Add click location to crds
    if but == 1 && size(crds,2)<2 
        crds(1,ct) = a; 
        crds(2,ct) = b; 
        ct = ct+1;
    end
    
    % Right click - Undo
    if but == 3
        cla % Clear and reset the image
        image(pics{1}{1}); colormap gray; 
        crds = []; 
        ct = 1;
    end

    % Plot the user selected line
    hold on; 
    try
        plot(crds(1,:),crds(2,:),'-g') 
        plot(crds(1,:),crds(2,:),'og')
    catch ploterror; 
    end
    hold off;
end

% Calculate initial penetration in mm
dist = crds(1,2) - crds(1,1)
mm_per_px = str2double(get(handles.mm_per_px, 'string'));
dist_mm = dist*mm_per_px
set(handles.i_penetration, 'String', dist_mm);

set(handles.instruction_box, 'String', ' ');


%%% SAVE BUTTON %%%
function save_button_Callback(hObject, eventdata, handles)
% Saves the data distance & time data in user specified location as a
% comma delimited .txt file with parameters in header 
% Note that original times and distancecs are used, without accounting for
% the initial penetration.

% Initiate many things
set(handles.instruction_box, 'String', 'Saving...');
global calc_dist;
global man_dist;
global times;
mm_per_px = str2double(get(handles.mm_per_px, 'string'));
i_frames = str2double(get(handles.i_frames, 'string'));
i_pen = str2double(get(handles.i_penetration, 'string'));
top_trim = str2double(get(handles.top_trim, 'string'));
bottom_trim = str2double(get(handles.bottom_trim, 'string'));
left_trim = str2double(get(handles.left_trim, 'string'));
right_trim = str2double(get(handles.right_trim, 'string'));
tta = str2double(get(handles.rot_angle, 'string'));
parameters = [mm_per_px, i_frames, i_pen, top_trim, bottom_trim, left_trim, right_trim, tta];
parameter_strings = {'mm_per_px', 'i_frames', 'i_pen', 'top_trim', 'bottom_trim', 'left_trim', 'right_trim', 'tta'};

% Ask user for save destination
[FileName,PathName] = uiputfile(...
    'C:\Users\George\OneDrive\University\Chemistry\Year 4 Chemistry\Project\Code\Data generated\*.txt');
full_path = fullfile(PathName,FileName);

% Populate the file
fid = fopen(full_path, 'w');  % Open a file for writing into
fprintf(fid, ['# Meniscus tracking data' 13 10]); % [13 10] is a line break
fprintf(fid, ['# Data from directory ']);
global file_folder
fprintf(fid, ['%s ' 13 10], file_folder);
% Write in parameters and their values
fprintf(fid, '# ');
for i = 1:length(parameters) % For the titles
    if i == length(parameters)
        fprintf(fid, ['%s ' 13 10], parameter_strings{i});
    else
        fprintf(fid, '%s, ', parameter_strings{i});
    end
end
fprintf(fid, '# ');
for i = 1:length(parameters) % For the values
    if i == length(parameters)
        fprintf(fid, ['%.3f' 13 10], parameters(i));
    else
        fprintf(fid, '%.3f, ', parameters(i));
    end
end

% Fill in any empty data with zeros
if isempty(calc_dist)
    calc_dist = zeros(length(times));
end
if isempty(man_dist)
    man_dist = zeros(length(times));
end

% Write in the times and distances
fprintf(fid, ['# times, calc_dist, man_dist' 13 10]);
for i = 1:length(times)
    fprintf(fid, ['%f, %f, %f' 13 10], [times(i), calc_dist(i), man_dist(i)]);
end

fclose(fid);  % Close the file
set(handles.instruction_box, 'String', ' ');


%%% DEFAULT CODE FOR GUI OBJECTS %%%
% --- Executes during object creation, after setting all properties.
function bottom_trim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bottom_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function view_time_Callback(hObject, eventdata, handles)
% hObject    handle to view_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of view_time as text
%        str2double(get(hObject,'String')) returns contents of view_time as a double

% --- Executes during object creation, after setting all properties.

function view_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mm_per_px_Callback(hObject, eventdata, handles)
% hObject    handle to mm_per_px (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mm_per_px as text
%        str2double(get(hObject,'String')) returns contents of mm_per_px as a double

% --- Executes during object creation, after setting all properties.

function mm_per_px_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mm_per_px (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rot_angle_Callback(hObject, eventdata, handles)
% hObject    handle to rot_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rot_angle as text
%        str2double(get(hObject,'String')) returns contents of rot_angle as a double

% --- Executes during object creation, after setting all properties.
function rot_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rot_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function top_trim_Callback(hObject, eventdata, handles)
% hObject    handle to top_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of top_trim as text
%        str2double(get(hObject,'String')) returns contents of top_trim as a double

% --- Executes during object creation, after setting all properties.
function top_trim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to top_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bottom_trim_Callback(hObject, eventdata, handles)
% hObject    handle to bottom_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bottom_trim as text
%        str2double(get(hObject,'String')) returns contents of bottom_trim as a double


% --- Executes on button press in sqrt_checkbox.
function sqrt_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to sqrt_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sqrt_checkbox


function left_trim_Callback(hObject, eventdata, handles)
% hObject    handle to left_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of left_trim as text
%        str2double(get(hObject,'String')) returns contents of left_trim as a double


% --- Executes during object creation, after setting all properties.
function left_trim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to left_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function right_trim_Callback(hObject, eventdata, handles)
% hObject    handle to right_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of right_trim as text
%        str2double(get(hObject,'String')) returns contents of right_trim as a double


% --- Executes during object creation, after setting all properties.
function right_trim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to right_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function i_penetration_Callback(hObject, eventdata, handles)
% hObject    handle to i_penetration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_penetration as text
%        str2double(get(hObject,'String')) returns contents of i_penetration as a double


% --- Executes during object creation, after setting all properties.
function i_penetration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_penetration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function i_frames_Callback(hObject, eventdata, handles)
% hObject    handle to i_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_frames as text
%        str2double(get(hObject,'String')) returns contents of i_frames as a double


% --- Executes during object creation, after setting all properties.
function i_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end