%%%% 3D holographic PTV code Sept 2014

clear all

level_flag = 1;
auto_centre_flag = 0;

folder = 'R:\Lab_Data\George Lewis\From Lisong\Holography Experiments\2015-04-10-drop-08-50um nozzle-3um PS-100H-A-2-60deg-good one\';
filename = ['water-',num2str(30),'.tif'];
im1 = imread([folder,filename]);

fim1 = fftshift(fft2(im1));

% Smooth Background - to level image intensity
if level_flag
    % smoothing with low pass - best and easiest way I have found
    Nlp = 5;    % Half width of low pass filter
    siz = size(im1);
    [lpx lpy] = meshgrid(-siz(2)/2:siz(2)/2-1,-siz(1)/2:siz(1)/2-1);
    lpmask = lpx.^2 + lpy.^2 < Nlp^2;
    sm = ifft2(fftshift(fftshift(fft2(im1)).*lpmask));
    sm = sm.^2;
    sm = sm-min(sm(:));
    sm = 1+sm./max(sm(:));
    clear lpx lpy lpmask
else
    sm = ones(size(im1));
end
%%
if auto_centre_flag
    DC_w=200;
    LP_n=5;
    thresh_m=1.5;
    flag_plot=1;
    
    [Cx,Cy,R]=find_synthcentre(fim1,DC_w,LP_n,thresh_m,flag_plot)
    f_cent=round([Cx Cy]); % centre of spectrum of interference term
    R=R*0.6;
else
    Cx= 176;
    Cy=96;
    R=100*0.8;
end
%%
s = size(im1);
[mx my] = meshgrid(1:s(2),1:s(1));

mask = (mx-Cx).^2 + ((s(2)/s(1))*(my-Cy)).^2 <R.^2;
imagesc(abs(fim1).*mask); axis equal; axis tight; colormap gray;

second_mask = 1;
if second_mask
    mask2 = (mx-114).^2 + ((s(2)/s(1))*(my-47)).^2 >30.^2;
    mask = mask & mask2;
    imagesc(abs(fim1).*mask)
end

%%
fim2 = fim1.*mask;
im2 = ifft2(fftshift(fim2));

%  create plane wave method-2
fx_c = (Cx - s(2)/2 - 1)/s(2);
fy_c = (Cy - s(1)/2 - 1)/s(1);
ref = exp(1i*(2*pi*(fx_c)*mx + 2*pi*(fy_c)*my));
demod = im2.*conj(ref);
imagesc(abs(demod))

fdemod = fftshift(fft2(demod));
imagesc(abs(fdemod));

%%
effpixel = 8.3e-6/50.8; % Camera pixel dimension 8.3microns, image magnification 48.
[A B] = refocuslite_companion(fdemod, 0.633e-6, effpixel); % wavelength 0.633um

zplane = refocuslite(fdemod,0e-6, A, B);% z varied in m
intensity=(abs(zplane)).^2; % x, y, unit in pixels
cla;
imagesc(abs(zplane)./sm); axis equal tight; colormap gray; colorbar
%disp('pausing...press any key to continue'); pause;

%% refocussing object field to each z section to generate 3D matrix block1=560x760xzs

z_step=0.1969;  % <-set this to exactly equal to the one in the tracking code
zs = (-30:z_step:30); % um ; aware a=55*1e-6 and b=a*1e6!=55
Nz = length(zs);

if mod(Nz,2)
    Nz = Nz-1;
    zs = zs(1:end-1);
end


block1 = zeros(size(fdemod,1),size(fdemod,2),Nz,'single');
[A,B] = refocuslite_companion(fdemod, 0.633e-6, effpixel);
block1(1) = 1+1i;
particles=[];

for ct = 1:Nz
    block1(:,:,ct) = refocuslite(fdemod,zs(ct)*1e-6, A, B); % reconstructed image
end
%%
Pc = [266,277];
Nrod = 8;
rod = block1(Pc(2)-Nrod/2:Pc(2)+Nrod/2-1,Pc(1)-Nrod/2:Pc(1)+Nrod/2-1,:);
plot(abs(squeeze(sum(sum((rod),1),2))))

pkind = find(rod == max(rod(:)),1,'first');
[i1 i2 i3] = ind2sub(size(rod),pkind);
p1 = polyfit(i1-1:i1+1,squeeze(abs(rod(i1-1:i1+1,i2,i3))).',2);
i1s = -p1(2)/(2*p1(1)); 
% pk1 = polyval(p1,i1s);
p2 = polyfit(i2-1:i2+1,squeeze(abs(rod(i1,i2-1:i2+1,i3))),2);
i2s = -p2(2)/(2*p2(1));
% pk2 = polyval(p2,i2s);
p3 = polyfit(i3-1:i3+1,squeeze(abs(rod(i1,i2,i3-1:i3+1))).',2);
i3s = -p3(2)/(2*p3(1));
% pk3 = polyval(p3,i3s);

frod = fftshift(fftn(rod));
lr = size(rod,3);
[fx,fy,fz] = meshgrid(-Nrod/2:Nrod/2-1,-Nrod/2:Nrod/2-1,-lr/2:lr/2-1);
fx = fx/Nrod;
fy = fy/Nrod;
fz = fz/lr;
k1 = 2*pi*(i1s-(Nrod/2+1));
k2 = 2*pi*(i2s-(Nrod/2+1));
k3 = 2*pi*(i3s-(lr/2+1));
fringes = exp(1i*k1*fy + 1i*k2*fx + 1i*k3*fz);
frod = frod.*fringes;
rod2 = ifftn(fftshift(frod));

P = rod2;
pfolder = 'R:\Lisong Yang\AndyFunctions\2015-04\27\';
pfilename = 'ideal_particle_2um.mat';
save([pfolder,pfilename],'P');
disp(['Model particle saved as: ',pfolder,pfilename])