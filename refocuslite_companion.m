function [sqrtR valid] = refocuslite_companion(fplane, lambda, px)

% Anderson McKeague. 2010
% Two dimensional refocussing.
%   refocuslite_companion(fplane, lambda, px)
%   fplane the same size as the hologram 
%   returns the sqrtR and valid for use in refocuslite.
%
% Run this function first to precalculate sqrtR and valid
% fplane - is the fftshift(fft2(*whatever-you-want-to-refocus*))
% lambda - is the wavelength of the light
% px     - is the effective pixel size in *whatever-you-want-to-refocus*
%
% This is propagation via a convolution with a phase map. Find details in
% Goodman's Fourier Optics.

sn = size(fplane);
% Make sure fplane and sn are compatible (ie both single or both double)
if  isa(fplane,'single')
    sn = single(sn); % Conversion to single rather than double precision
                     % takes up less memory space since it is accurate to
                     % ~10E-7 instead of ~10E-16 or so.
end

[fx fy] = meshgrid(-sn(2)/2:sn(2)/2-1,-sn(1)/2:sn(1)/2-1);
fx = fx./(sn(2)*px); fy = fy./(sn(1)*px);

sqrtR = ((1/lambda).^2)-(fx.^2)-(fy.^2);
valid = sqrtR<0; % Either gives 1 or 0 if +ve or -ve
sqrtR = sqrt(sqrtR); % Parameter in propagation eqn (p119 A. McKeague thesis)












