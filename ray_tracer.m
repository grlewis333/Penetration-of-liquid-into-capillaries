function varargout = ray_tracer(varargin)
% RAY_TRACER MATLAB code for ray_tracer.fig
%      This code is for a ray tracing GUI geared towards modelling light
%      that passes through capillaries. It is capable of modelling point 
%      sources of light at any given point in the capillary for any given
%      refractive indices. It can also analyse the detection of these rays.
%      
%      Contents:
%      1. Ray generator: creates point source with angular spread of rays
%      2. Plot button: propagates rays through capillary
%      3. Analyse button: analyses frequency of detection in capillary
%
%      Note that analysis of rays reaching the detector which are limited
%      by numerical aperture is not yet fully incorporated into the GUI and
%      the program is currently setup to analyse a specific instance.
%
%      Default help code:
%      RAY_TRACER, by itself, creates a new RAY_TRACER or raises the existing
%      singleton*.
%
%      H = RAY_TRACER returns the handle to a new ray_tracer or the handle to
%      the existing singleton*.
%
%      RAY_TRACER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAY_TRACER.M with the given input arguments.
%
%      RAY_TRACER('Property','Value',...) creates a new RAY_TRACER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ray_tracer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ray_tracer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ray_tracer

% Last Modified by GUIDE v2.5 03-Apr-2018 19:15:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ray_tracer_OpeningFcn, ...
                   'gui_OutputFcn',  @ray_tracer_OutputFcn, ...
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


% --- Executes just before ray_tracer is made visible.
function ray_tracer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ray_tracer (see VARARGIN)
%WRITE STUFF AFTER THIS

% Choose default command line output for ray_tracer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ray_tracer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ray_tracer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%% Ray generator %%%
function rays = ray_generator(xii, yii, no_rays)
% This function will create a user defined number of rays with evenly 
% spaced angles, coming from a point source
% Generates rays in form {{x1,y1,theta1},{x2,y2,theta2}...}
rays = {};
angular_spread = pi/2;
angles = linspace(-angular_spread/2, angular_spread/2, no_rays);
% Note the pol2cart function used later defines theta counter c/w from +x

if no_rays == 1
    for i = 1:length(angles) % loops through to populate rays
        rays{1,i} = [xii,yii,0.0];
    end
else
    for i = 1:length(angles) % loops through to populate rays
        item = angles(i);
        rays{1,i} = [xii,yii,item];
    end
end
    

%%% Plot button %%%
function plot_button_Callback(hObject, eventdata, handles)
% Plots rays propagating out from a point source within the capillary

% Clear anything already on the graph
axes(handles.axes1); 
cla reset; 

% Reset error bar
error_message = ' ';
set(handles.error_bar, 'String', error_message);

hold on % To ensure all lines are plotted on same axes

% Access user defined inputs from GUI 
xii = str2double(get(handles.xi, 'string')); % Starting point
yii = str2double(get(handles.yi, 'string'));

no_rays = str2double(get(handles.no_rays, 'string')); % Number of rays

int_d = str2double(get(handles.int_d, 'string')); % Capillary diameters
ext_d = str2double(get(handles.ext_d, 'string'));

square_switch = get(handles.sq_cap_check, 'Value') % Is it a square cap?

n1 = str2double(get(handles.n_int, 'string')); % Refractive indices
n2 = str2double(get(handles.n_cap, 'string'));
n3 = str2double(get(handles.n_ext, 'string'));

x4 = str2double(get(handles.x_det, 'string')); % Detector specification
h_det = str2double(get(handles.detector_height, 'string'));
fov = str2double(get(handles.fov, 'string'));

% Call function to generate rays
rays = ray_generator(xii, yii, no_rays);

% Propagate the rays to the inner capillary wall
ray_innerWall = {};

% Calculate intersection of ray with inner capillary wall
%{
% Generate array with positions at inner wall {[x1, y1], [x2, y2]...}
% Note this way disperses rays over an angular spread *always* defined at 
% the origin, the alternative way below generates the same angular spread 
% but always from the same given point.

for i = 1:length(rays) 
    ray = rays{1,i};
    theta = ray(1,3);
    [x,y] = pol2cart(theta,int_d/2);
    ray_innerWall{1,i} = [x,y,theta];   
end 
%}

% Better way of generating positions at inner wall
if square_switch == 1; % For square capillaries
    for i = 1:length(rays) 
    ray = rays{1,i};
    theta = ray(1,3);
    m = tan(theta); % Caluclate straight line parameters of ray
    c = yii-m*xii;
    x = int_d/2 % x position is always the same for square cap
    y = m*x+c % y position is simply calculated from intersection of ray
    ray_innerWall{1,i} = [x,y,theta];
    end
    
else
for i = 1:length(rays) % For circular capillaries
    ray = rays{1,i}; 
    theta = ray(1,3);
    m = tan(theta); % Caluclate straight line parameters of ray
    c = yii-m*xii;
    % A, B and C define coefficients for eqn describing intersection of 
    % line (ray) with circle (capillary)
    A = m^2+1;  
    B = 2*m*c;
    C = -(int_d/2)^2+c^2;
    x = (-B + sqrt((B^2-4*A*C)/2*A))/(2*A); % x and y from intersection
    y = m*x + c;
    ray_innerWall{1,i} = [x,y,theta];
end
end

% Plot rays from origin to inner capillary
for i = 1:length(rays)  
    iray = rays{1,i};
    xi = iray(1,1);
    yi = iray(1,2);
    fray = ray_innerWall{1,i};
    xf = fray(1,1);
    yf = fray(1,2);
    
    % RGB can be used to define a color tuple (55, 55, 55 is grey)
    % Must be normalised to 1 for Matlab to like it
    % More colors at http://www.rapidtables.com/web/color/RGB_Color.htm
    R = 55;
    G = 55;
    B = 55;
    
    plot([xi,xf],[yi,yf], 'color', 1/(R+G+B)*[R,G,B])   
end

% Plot the capillary
if square_switch == 1 % If square
    for i = 1:2
        if i == 1 % Plot inner wall
            x = int_d/2;
            y1 = -int_d/2;
            y2 = int_d/2;
        else      % Plot outer wall
            x = ext_d/2;
            y1 = -ext_d/2;
            y2 = ext_d/2;
        end
        plot([x,x],[y1,y2],'b-'), axis equal
        plot([0,x],[y1,y1],'b-')
        plot([0,x],[y2,y2],'b-')
    end
    
else   % Plot circular capillaries
    for i = 1:2
        if i == 1 % Plot inner wall
            r = int_d/2;
        else      % Plot outer wall
            r = ext_d/2;
        end

        % Point 1 and 2 are the bottom and top points of the capillary
        x1 = 0;
        x2 = 0;
        y1 = -r;
        y2 = r;

        d = sqrt((x2-x1)^2+(y2-y1)^2); % Distance between points
        a = 0; % Perpendicular bisector angle (atan(-(x2-x1),(y2-y1)) if tilt)
        b = asin(d/2/r); % Half arc angle
        c = linspace(a-b,a+b); % Arc angle range
        e = sqrt(r^2-d^2/4); % Distance, center to midpoint
        x = (x1+x2)/2-e*cos(a)+r*cos(c); % Cartesian coords. of arc
        y = (y1+y2)/2-e*sin(a)+r*sin(c);
        plot(x,y,'b-'), axis equal
    end
end

% Propagate from inner to outer edge
if square_switch == 1 % For square capillaries
    for i = 1:length(rays)
    % Ray1 is the initial point
    ray1 = rays{1,i};
    x1 = ray1(1,1);
    y1 = ray1(1,2);
    
    % Ray2 is the point at the inner wall
    ray2 = ray_innerWall{1,i};
    x2 = ray2(1,1);
    y2 = ray2(1,2);
    
    
    
    
    % Theta1 and 2 are the initial and 1st modified angles
    theta1 = atan((y2-y1)/(x2-x1));
    theta2 = asin(n1*sin(theta1)/n2); % Snell's law
       
    % Gradient and intercept of ray leaving the inner edge (towards outer)
    m = tan(theta2);
    c = y2 - m*x2;
    
    % Point 3 is where the rays meet the outer capillary edge
    x3 = ext_d/2;
    y3 = m*x3+c;
    rays3{1,i} = [x3,y3,theta2];
    
    
    % RGB can be used to define a color tuple (55, 55, 55 is grey)
    % Must be normalised to 1 for Matlab to like it
    % More colors at http://www.rapidtables.com/web/color/RGB_Color.htm
    R = 55;
    G = 55;
    B = 55;
    
    % If internal reflection occurs, theta3 will have an imaginary
    % component
    theta3 = asin(n2*sin(theta2)/n3);
    if imag(theta3) ~= 0
         R = 1;
         G = 0;
         B = 0;
         error_message = 'Total internal reflection is occuring for red lines';
         set(handles.error_bar, 'String', error_message);
    end
    
    NA_switch = 0
    if get(handles.detector_analysis_switch,'Value') ~= 0
        if abs(y2) > fov/2/1000 % Is the ray limited by NA from reaching the 
                                % detector?
             R = 1;
             G = 0;
             B = 0;
             error_message = 'NA';
             set(handles.error_bar, 'String', error_message);
             theta3=1i
            NA_switch = 1 % This will prevent the ray being 'counted'
        end
    end
    
    plot([x2,x3],[y2,y3], 'color', 1/(R+G+B)*[R,G,B]);
    end
    
else   % Plot inner to outer for circular capillaries
for i = 1:length(rays)
    NA_switch = 0
    % Ray1 is the initial point
    ray1 = rays{1,i};
    x1 = ray1(1,1);
    y1 = ray1(1,2);
    
    % Ray2 is the point at the inner wall
    ray2 = ray_innerWall{1,i};
    x2 = ray2(1,1);
    y2 = ray2(1,2);
    
    % Theta1 and 2 are the initial and 1st modified angles
    theta1 = atan((y2-y1)/(x2-x1));
    theta2 = asin(n1*sin(theta1)/n2); % Snell's law
       
    % Gradient and intercept of ray leaving the inner edge (towards outer)
    m = tan(theta2);
    c = y2 - m*x2;
    
    % Coefficients of quadratic for line/circle intersection, defined at 
    % https://math.stackexchange.com/questions/228841/how-do-i-calculate-the-intersections-of-a-straight-line-and-a-circle
    A = m^2 + 1;
    B = 2*(m*c);
    C = (c^2 - (ext_d/2)^2);
    
    % Point 3 is where the rays meet the outer capillary edge
    x3 = (-B + (B^2 - 4*A*C)^0.5)/(2*A);
    y3 = m*(-B + (B^2 - 4*A*C)^0.5)/(2*A) + c;
    rays3{1,i} = [x3,y3,theta2];
    
    
    % RGB can be used to define a color tuple (55, 55, 55 is grey)
    % Must be normalised to 1 for Matlab to like it
    % More colors at http://www.rapidtables.com/web/color/RGB_Color.htm
    R = 55;
    G = 55;
    B = 55;
    
    % If internal reflection occurs, theta3 will have an imaginary
    % component
    theta3 = asin(n2*sin(theta2)/n3);
    if imag(theta3) ~= 0
         R = 1; % Make the line red
         G = 0;
         B = 0;
         error_message = 'Total internal reflection is occuring for red lines';
         set(handles.error_bar, 'String', error_message);
    end
    
    plot([x2,x3],[y2,y3], 'color', 1/(R+G+B)*[R,G,B]);
end
end

% Propagate to the detector
for i = 1:length(rays)
    % Point 3 is where the ray meets the outer edge
    ray3 = rays3{1,i};
    x3 = ray3(1,1);
    y3 = ray3(1,2);
    
    % As before, theta2 is the first modified angle, 3 is the second
    theta2 = ray3(1,3);
    theta3 = asin(n2*sin(theta2)/n3);
    
    % Gradient and intercept of ray leaving outer edge
    m = tan(theta3);
    c = y3 - m*x3;
    
    % y intersection of ray with the detector
    y4 = m*x4 + c;
    if imag(y4) == 0
        y4s{i} = y4;
    end
    
    % Stop rays limited by NA/internal reflection
    NA_switch = 0;
    ray2 = ray_innerWall{1,i};
    y2 = ray2(1,2);
    if get(handles.detector_analysis_switch,'Value') ~= 0
        if abs(y2) > fov/2/1000
            NA_switch = 1;
            y4 = h_det;
            y4s{i} = y4;
        end
    end
    if imag(theta3) ~= 0
        y4 = h_det;
    end
    
    
    % RGB can be used to define a color tuple (55, 55, 55 is grey)
    % Must be normalised to 1 for Matlab to like it
    % More colors at http://www.rapidtables.com/web/color/RGB_Color.htm
    R = 55;
    G = 55;
    B = 55;
    
    % Plot so long as internal reflection hasn't occurred 
    if imag(theta3) == 0
        if NA_switch ~= 1 % Stop in case limited by NA
        plot([x3,x4],[y3,y4], 'color', 1/(R+G+B)*[R,G,B]);
        end
    end
end

% Set axis limits and labels
%axis([0 x4 -h_det/2 h_det/2])%(-ext_d/2 - ext_d/20) (ext_d/2 + ext_d/20)])
axis([0 x4 (-ext_d/2 - ext_d/20) (ext_d/2 + ext_d/20)])
xlabel('x / mm')
ylabel('y / mm')

% Determine number of rays hitting detector
if get(handles.detector_analysis_switch,'Value') ~= 0 % If analysis 'on'
    plot([x4,x4],[-h_det/2,h_det/2], 'color', 'red')
    hits = 0;
    if length(y4s) == 0
        hits = 0;
    else
        for i = 1:length(y4s)
            y4 = y4s{i};

            if y4 < h_det/2 & y4 > -h_det/2 % Check if ray is hitting detector
                hits = hits + 1; % Score as a 'hit' if it is
                y4;
            end
        end
    end
    global hit_count % Make global for use in analysis code
    hit_count = hits
end
hold off

%%% Analyse button %%%
function analyse_button_Callback(hObject, eventdata, handles)
% This function creates a heatmap showing the relative detection of rays
% initiating from a series of evenly spaced points within the capillary

% Obtain inputs from GUI
square_switch = get(handles.sq_cap_check, 'Value'); % see if square or circ
set(handles.detector_analysis_switch, 'Value', 1); % Switch analysis on
int_d = str2double(get(handles.int_d, 'string')); % Check diameter

% Initialise some other variables
division = 0.01; % Value determines spacing between the analysed points
xs = -int_d/2+division:division:(int_d/2-division); % Generate x,y for analysis points
ys = -int_d/2+division:division:int_d/2-division;

all_hits = []; 
j = 0;
results={};
xvalues = [];
yvalues = [];
cdata = [];

% Cycle through points and determine hit count for each
for i = 1:length(xs)
    x = xs(i);
    for i2 = 1:length(ys)
        y = ys(i2);
        if square_switch == 1 % If square
            j = j+1;
            set(handles.xi, 'String', x); % Change ray start in GUI
            set(handles.yi, 'String', y);
            plot_button_Callback(hObject, eventdata, handles); % Call func.
            global hit_count;
            hits = hit_count;
            results{1,j} = [x,y,hits]; % Put hit-count into results array
        
        else   % If circular
            if sqrt(x^2+y^2)<int_d/2 % Discard points outside inner cap.
                j = j+1;
                set(handles.xi, 'String', x); % Change ray start in GUI
                set(handles.yi, 'String', y);
                plot_button_Callback(hObject, eventdata, handles); % Call
                global hit_count; 
                hits = hit_count;
                results{1,j} = [x,y,hits]; % Hit-count into results array
            end
        end
    end
end

% Extract data from results array into simple lists
for i = 1:length(results)
    res = results{1,i};
    x = res(1,1);
    y = res(1,2);
    hits = res(1,3);
    xvalues = [xvalues, x];
    yvalues = [yvalues, y];
    cdata = [cdata, hits];
end
cdata = cdata*100/max(cdata)


% Plot results
colormap hot % choose a colormap of your liking
pointsize = 25;
scatter(xvalues*1000, yvalues*1000, pointsize, cdata, 'filled')
daspect([1 1 1]);

% This will display the make the figure also display in a separate window
% so that it can be saved
figure, scatter(xvalues*1000, yvalues*1000, pointsize, cdata, 'filled')%,'MarkerEdgeColor',[0.1 0.11 0.12])
xlabel('Optical axis / µm')
ylabel('Vertical axis / µm')
xlim([-100,100])
ylim([-100,100])
%colormap(flipud(jet))
colormap jet
pointsize = 25;
colorbar()


%%% DEFAULT CODE FOR GUI OBJECTS %%%
function int_d_Callback(hObject, eventdata, handles)
% hObject    handle to int_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_d as text
%        str2double(get(hObject,'String')) returns contents of int_d as a double

% --- Executes during object creation, after setting all properties.
function int_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ext_d_Callback(hObject, eventdata, handles)
% hObject    handle to ext_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ext_d as text
%        str2double(get(hObject,'String')) returns contents of ext_d as a double


% --- Executes during object creation, after setting all properties.
function ext_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ext_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xi_Callback(hObject, eventdata, handles)
% hObject    handle to xi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xi as text
%        str2double(get(hObject,'String')) returns contents of xi as a double


% --- Executes during object creation, after setting all properties.
function xi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yi_Callback(hObject, eventdata, handles)
% hObject    handle to yi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yi as text
%        str2double(get(hObject,'String')) returns contents of yi as a double

% --- Executes during object creation, after setting all properties.
function yi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x_det_Callback(hObject, eventdata, handles)
% hObject    handle to x_det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_det as text
%        str2double(get(hObject,'String')) returns contents of x_det as a double

% --- Executes during object creation, after setting all properties.
function x_det_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n_int_Callback(hObject, eventdata, handles)
% hObject    handle to n_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_int as text
%        str2double(get(hObject,'String')) returns contents of n_int as a double

% --- Executes during object creation, after setting all properties.
function n_int_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n_cap_Callback(hObject, eventdata, handles)
% hObject    handle to n_cap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_cap as text
%        str2double(get(hObject,'String')) returns contents of n_cap as a double

% --- Executes during object creation, after setting all properties.
function n_cap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_cap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n_ext_Callback(hObject, eventdata, handles)
% hObject    handle to n_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_ext as text
%        str2double(get(hObject,'String')) returns contents of n_ext as a double

% --- Executes during object creation, after setting all properties.
function n_ext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function no_rays_Callback(hObject, eventdata, handles)
% hObject    handle to no_rays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of no_rays as text
%        str2double(get(hObject,'String')) returns contents of no_rays as a double

% --- Executes during object creation, after setting all properties.
function no_rays_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_rays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in detector_analysis_switch.
function detector_analysis_switch_Callback(hObject, eventdata, handles)
% hObject    handle to detector_analysis_switch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detector_analysis_switch

function detector_height_Callback(hObject, eventdata, handles)
% hObject    handle to detector_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of detector_height as text
%        str2double(get(hObject,'String')) returns contents of detector_height as a double

% --- Executes during object creation, after setting all properties.
function detector_height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to detector_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in sq_cap_check.
function sq_cap_check_Callback(hObject, eventdata, handles)
% hObject    handle to sq_cap_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sq_cap_check



function NA_Callback(hObject, eventdata, handles)
% hObject    handle to NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NA as text
%        str2double(get(hObject,'String')) returns contents of NA as a double


% --- Executes during object creation, after setting all properties.
function NA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fov_Callback(hObject, eventdata, handles)
% hObject    handle to fov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fov as text
%        str2double(get(hObject,'String')) returns contents of fov as a double


% --- Executes during object creation, after setting all properties.
function fov_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
