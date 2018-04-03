%% Find and read file - plot in Fourier space
clear
    % Allow user to specify file
    [FileName,PathName,FilterIndex] = uigetfile(...
        'C:\Users\George\OneDrive\University\Chemistry\Year 4 Chemistry\Project\Drive data\HPIV\Capillary test\Last week\T8\holo vid 2\*.tif');
    file_path = strcat(PathName, FileName);

%% Plot in Fourier space
    % Read image
    im = imread(file_path);
    imagesc(im); colormap gray
    
    fim = fftshift(fft2(im)); % fftshift makes zero-frequency component
                              % be in the centre. fft2 = 2D FT
                              % When masked, this will be R*O
    
    image(abs(fim)/1e2)
    axis equal tight;
%% Determine location phase information and mask everything else
% Input centre of the off-axis phase information   
Cx = 268
Cy = 114
R=100*0.8; % R is radius of the phase information region
% according to the theory, the higher frequency the better 
% resolution, but we found with a ratio slightly less than 1
% (normally use 0.8),we got better signal-to-noise ratio. 
    
s = size(im);
[mx my] = meshgrid(1:s(2),1:s(1)); % Create arrays the size of the image
% (Meshgrid is important for interpolation which reduces aliasing?)

% Create 'mask' (open circle around phase info)
mask = (mx-Cx).^2 + ((s(2)/s(1))*(my-Cy)).^2 <R.^2;
% Eqn of circle translated to phase info for values within radius

imagesc(abs(fim).*mask); axis equal; axis tight % Multply Fourier image & mask
% A.*B is element by element multiplication of A and B
% A.^B raises each element of A to the power B

%% Extraction of object phase term (masking + demodulation)
fim2 = fim.*mask;
im2 = ifft2(fftshift(fim2)); % This is extracted and recentred cross term
                             % i.e. R* O

% Create plane wave
% s is simply (560,760) and Cx - 760/2, Cy - 560/2 is just a conversion
% from 'pixel coordinate' to spatial frequency (i.e. wrt the centre of the
% Fourier image).
% The factors of 1/560 and 1/760 are there since a Discrete Fourier
% transform with 560 x 760 samples has been performed (search DFT for info)
fx_c = (Cx - s(2)/2 - 1)/s(2);
fy_c = (Cy - s(1)/2 - 1)/s(1);
ref = exp(1i*(2*pi*(fx_c)*mx + 2*pi*(fy_c)*my));

demod = im2.*conj(ref); %|R|^2 O
imagesc(abs(demod))

fdemod = fftshift(fft2(demod)); 
imagesc(abs(fdemod));
%% Propagate through surface
glass_capillary = 0 % 'switch' for this section
if glass_capillary == 1
    effpixel = 8.3e-6/50;
    thick = 10e-6; % thickness in meter
    n_cs = 1.5; % refractive index
    [A,B] = refocuslite_companion(fdemod, 0.633e-6/n_cs, effpixel);
    demod = refocuslite(fdemod,thick, A, B);
    fdemod = fftshift(fft2(demod));
    imagesc(abs(fdemod));
end

%% Propagate and display images over selected planes
effpixel = 8.3e-6/50.8 % effective pixel is actual pixel size in camera
                       % (from operation manual), divided by magnification
                       % (calculated from image of calibration dots)
                       
% Calculate some of the values needed in eqn's of propagation
through_water = 1;
if through_water == 1;
    [A B] = refocuslite_companion(fdemod, 0.633e-6/1.33, effpixel);
else
    [A B] = refocuslite_companion(fdemod, 0.633e-6, effpixel);
end

% Display reconstructed image focused over range of z planes
ii = 1;
for z=-100:2:100% Display planes focussed in this z range - z in um
    % Propagate to new plane
    zplane = refocuslite(fdemod,z*1e-6, A, B);
    intensity=(abs(zplane)).^2; % x, y, unit in pixels
    fig = figure(1)
    %imagesc(angle(zplane)); axis equal tight; colormap gray; colorbar
    imagesc(abs(zplane)); axis equal tight; colormap gray; colorbar
    text(400,400,['\color{yellow}',num2str(z),' \mum']);
    
    % Saves images in current directory (for later creating video)
    %{
    name = char(strcat('im',string(ii),'.tif'));
    saveas(fig,name);
    ii = ii + 1;
    %}
end
