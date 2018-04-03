
list = [];
a=[];
for N = 16:120;
    
      folder = 'Y:\Lab_Data\Lisong Yang\AndyFunctions\2015-04\10\drop-08-50um nozzle-3um PS-100H-A-2-60deg-good one\';
      filename = ['water-',num2str(N),'.tif'];
      im1 = imread([folder,filename]);
    %im1 = imread(['Y:\Lab_Data\Lisong Yang\AndyFunctions\2014-09-22\drop-good04-30um nz\3um-',num2str(N),'.tif']);
%     im1= imread(['Y:\Lab_Data\Lisong Yang\AndyFunctions\2015-04\10\drop-08-50um nozzle-3um PS-100H-A-2-60deg-good one\water-',num2str(N),'.tif']);
   
%     imagesc(im1); colormap gray
    
    % im = double(im)-imtot;
    fim1 = fftshift(fft2(im1));
     fim1(:,1*end/4:end) = 0;
     fim1(1*end/4:end,:) = 0;
 
    %lis image(abs(fim1)/11e1); pause
    list(end+1) = sum(abs(fim1(:)));
    if list(end) > 4.0e7
        a=[a N];
        disp(['Good (',num2str(N),').']);
         image(abs(fftshift(fft2(im1)))/11e1); 
         colormap(gray);
    end
     % pause;%(0.01)
end

% Ns = 1:length(list);
% Ns(list > 2.4e7)