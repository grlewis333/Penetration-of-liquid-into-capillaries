%%% 3D holography particle identification code
clear all;
goodholos = [0]%,1,2,3,4,5,6,7,8];

xyzs = [];
level_flag = 1; % Set = 1 for smoothing of images

% Load in the 'ideal particle' file
load('C:\Users\George\OneDrive\University\Chemistry\Year 4 Chemistry\Project\Copy of Drive\From Lisong\Holography Experiments\2015-04-10-drop-08-50um nozzle-3um PS-100H-A-2-60deg-good one\matlab\ideal_particle_2um.mat','P')

% Choose file/folder containing files
[FileName,PathName,FilterIndex] = uigetfile(...
    'R:\Lab_Data\George Lewis\HPIV\*.tif'); 
file_path = strcat(PathName, FileName);

%%
global all_xyzs
all_xyzs = {};
for N = 1:length(goodholos)
   % Reading in
    disp(['processing hologram No.',num2str(goodholos(N)),'...'])  

    % Cycle through holograms if looking at more than 1
    if length(goodholos) > 1
        % Specify file naming format here
        FileName = strcat('holo', string(goodholos(N)), '.tif');
        file_path = char(strcat(PathName, FileName));
    end
    
    im1 = imread(file_path);
    fim1 = fftshift(fft2(im1)); % Fourier transform and centre
    
    % Smooth Background - to level image intensity
    if level_flag == 1
        % smoothing with low pass
        Nlp = 5;    % Half width of low pass filter
        siz = size(im1);
        [lpx, lpy] = meshgrid(-siz(2)/2:siz(2)/2-1,-siz(1)/2:siz(1)/2-1);
        lpmask = lpx.^2 + lpy.^2 < Nlp^2;
        sm = ifft2(fftshift(fftshift(fft2(im1)).*lpmask));
        sm = sm.^2;
        sm = sm-min(sm(:));
        sm = 1+sm./max(sm(:));
        clear lpx lpy lpmask
        imagesc(sm)
    else
         sm = ones(size(im1));
    end
%% Input centre of the off-axis phase information 
% This should come from a reference hologram with no object in place
Cx= 240;
Cy=132;
R=100*0.8;

%% Mask everything but cross-phase term
s = size(im1);
[mx, my] = meshgrid(1:s(2),1:s(1));

mask = (mx-Cx).^2 + ((s(2)/s(1))*(my-Cy)).^2 <R.^2;
imagesc(abs(fim1).*mask); axis equal; axis tight; colormap gray;

second_mask = 1;
if second_mask
    mask2 = (mx-114).^2 + ((s(2)/s(1))*(my-47)).^2 >30.^2;
    mask = mask & mask2;
    imagesc(abs(fim1).*mask)
end

%% Extract field of object
fim2 = fim1.*mask;
im2 = ifft2(fftshift(fim2));

%  create plane wave method-1
% fref = zeros(s);
% fref(f_cent(2),f_cent(1)) = 1;
% fref(f_cent(2),f_cent(1)) = fim(f_cent(2),f_cent(1));
% ref = ifft2(fftshift(fref)); 

%  create plane wave method-2
fx_c = (Cx - s(2)/2 - 1)/s(2);
fy_c = (Cy - s(1)/2 - 1)/s(1);
ref = exp(1i*(2*pi*(fx_c)*mx + 2*pi*(fy_c)*my));
%imagesc(real(ref))
demod = im2.*conj(ref);
imagesc(abs(demod))

fdemod = fftshift(fft2(demod));
imagesc(abs(fdemod));
%     cla;    
%     imagesc(abs(fim1).*mask); axis equal; 
 %% Prepare values for use later
effpixel = 8.3e-6/50.8; % Camera pixel dimension 8.3 microns, 
                        % image magnification 48.
[A, B] = refocuslite_companion(fdemod, 0.633e-6/1.33, effpixel);  

zplane = refocuslite(fdemod,0e-6, A, B);% z varied in m
intensity=(abs(zplane)).^2; % x, y, unit in pixels
cla;
imagesc(abs(zplane)./sm); axis equal tight; colormap gray; colorbar

% Get images to display later for checking
zplane = refocuslite(fdemod,0e-6, A, B);
org_im = abs(zplane)./sm;

zplane = refocuslite(fdemod,-112e-6, A, B);
org_im2 = abs(zplane)./sm;
imagesc(org_im2)

%% Refocus object field for each z to generate 3D matrix block1=560x760 px
% then detect particles for each zs to generate 2D matrix as 
% particles = [x y z] = particle number x 3

z_step=0.7; % Specify distance betweeen z planes in um
zs = (-135:z_step:-90);
Nz = length(zs);
% Note, if you get a 'Subscripted assignment dimension mismatch' error
% then try a smaller zs range

% Just quickly make sure that Nz is a good number.
while mod(Nz,4)
    Nz = Nz+1;
end
zs = linspace(min(zs),max(zs),Nz);
z_step = zs(2)-zs(1);

% block1 is an array holding the image of each refocused z plane
block1 = zeros(size(fdemod,1),size(fdemod,2),Nz,'single');
[A,B] = refocuslite_companion(fdemod, 0.633e-6/1.33, effpixel);
block1(1) = 1+1i;
particles=[];

% Create reconstructed images for each z
for ct = 1:Nz
    block1(:,:,ct) = imcomplement(refocuslite(fdemod,zs(ct)*1e-6, A, B)); 
end

Np = size(P); % P is defined in the initially imported 'ideal particle'
sb = size(block1);
Pbig = zeros(sb,'single');
if Np(3) > sb(3)
    Pbig((sb(1)/2+1-Np(1)/2:sb(1)/2+Np(1)/2),(sb(2)/2+1-Np(2)/2:sb(2)/2+Np(2)/2),:) ...
        = P(:,:,(Np(3)/2+1-sb(3)/2:Np(3)/2+sb(3)/2));
else
    Pbig((sb(1)/2+1-Np/2:sb(1)/2+Np/2),(sb(2)/2+1-Np/2:sb(2)/2+Np/2),(1:sb(3))) = P;
end

Pbig = fftn(Pbig);
fblock = fftn(single(bsxfun(@times,block1,max(sm(:))./sm)));
Pbig = conj(Pbig).*(fblock);
clear fblock
Pbig = ifftn(Pbig);
Pbig = fftshift(Pbig);

I = sum(abs(Pbig),3);
I = bpass(I,1,11);  %filter, second argument: lengthscale of noise in pixels 
                    %third argument: integer length in pixels larger than a typical object
I = I./max(I(:));

pks_2D=pkfnd(I,0.5,7);  % second argument: intensity level-threhold, the third: the size of average blobs

hold on
plot(pks_2D(:,1),pks_2D(:,2),'oy')
hold off

particles=particles;
cnt = pks_2D;
pks_z = zeros(length(pks_2D),1);
%% Ensure program doesn't crash
% If <1 particle detected, next section will fail without this
 try
     for ct = 1:length(pks_2D)
         pks_2D(ct,2)
     end
 catch
     continue
 end
 %%
        
for ct = 1:length(pks_2D)
    i1 = pks_2D(ct,2);
    i2 = pks_2D(ct,1);
    st1 = i1-7;
    st2 = i2-7;
    en1 = i1+7;
    en2 = i2+7;
    if st1 < 1; st1 = 1; end
    if st2 < 1; st2 = 1; end
    if en1 > sb(1); en1 = sb(1); end
    if en2 > sb(2); en2 = sb(2); end
    L1 = st1:en1;
    L2 = st2:en2;
    rod = block1(L1,L2,:);
%     [iall] = find(real(rod) == min(real(rod(:))));
    [iall] = find(abs(rod) == max(abs(rod(:))));
    [su1,su2,su3] = ind2sub(size(rod),iall);
    pks_z(ct) = min(su3);
end

% Plot particles with image of desired focal plane
imagesc(sum(abs(block1),3))
imagesc(org_im2); colormap gray;
hold on
plot3(pks_2D(:,1),pks_2D(:,2),pks_z,'ob')
hold off

particles=[particles; [round(cnt(:,1)) round(cnt(:,2)) pks_z]]; %note x,y are not integer

%% Make a 3D plot of all the detected particles
     % xyzs=[x y z holo_number]
     xyzs = []
     xyzs = cat(1,xyzs,[particles(:,1),particles(:,2),particles(:,3),ones(size(particles(:,1)))*goodholos(N)]);
     cla;
     % Plot in um
     %plot3(xyzs(:,1)*effpixel*1e6,xyzs(:,2)*effpixel*1e6,zs(xyzs(:,3)),'b*','markersize',3);
     % Plot in px
     plot3(xyzs(:,1),xyzs(:,2),zs(xyzs(:,3)),'ro','markersize',3);

     set(gca,'YDir','reverse');
                xlabel('x /\mum');
                ylabel('y /\mum');
                zlabel('z /\mum');
                axis normal;
     disp(['  ',num2str(round(200*N/length(goodholos))/2),' % completed.']);
     hold on
     imagesc(org_im2); colormap gray;
     hold off
end

