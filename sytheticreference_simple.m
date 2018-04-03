clear
    im = imread('3um-698-holo.tif'); 
    imagesc(im); colormap gray
    im=im(:,:,1);

    fim = fftshift(fft2(im));
    image(abs(fim)/1e2)
    axis equal tight;
%%
auto_centre_flag=0;
if auto_centre_flag

    DC_w=200;
    LP_n=5;
    thresh_m=1.3;
    flag_plot=1;

    [Cx,Cy,R]=find_synthcentre(fim,DC_w,LP_n,thresh_m,flag_plot)
    f_cent=round([Cx Cy]); % centre of spectrum of interference term
else 
%     Cx = 145;
%     Cy = 395;
%     R=100*0.8; % note: if you use left-lower corner, varied the construct
%                % z- sign

    Cx = 623;
    Cy = 173;%104;
    R=100*0.8; % according to the theory, the higher frequency the better 
              % resolution, but we found with a ratio slightly less than 1
              % (normally use 0.8),we got better signal-to-noise ratio. 
end

s = size(im);
[mx my] = meshgrid(1:s(2),1:s(1));

mask = (mx-Cx).^2 + ((s(2)/s(1))*(my-Cy)).^2 <R.^2;
imagesc(abs(fim).*mask); axis equal; axis tight

%%

fim2 = fim.*mask;
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
%%
effpixel = 8.3e-6/50.8;
[A B] = refocuslite_companion(fdemod, 0.633e-6, effpixel);
for z=50:1:80
zplane = refocuslite(fdemod,z*1e-6, A, B);
intensity=(abs(zplane)).^2; % x, y, unit in pixels
figure(1)
%imagesc(angle(zplane)); axis equal tight; colormap gray; colorbar
imagesc(abs(zplane)); axis equal tight; colormap gray; colorbar
text(100,100,['\color{yellow}',num2str(z),' \mum']);
%pause(0.1)
end
