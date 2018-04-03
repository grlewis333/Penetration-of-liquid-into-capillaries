function [newplane] = refocuslite(fplane, distance, sqrtR, valid)

% Anderson McKeague. 2010
% refocuslite faster two dimensional refocussing.
%   refocuslite(fplane, distance, sqrtR, valid)
%   fplane must be fftshift(fft2(*whatever-you-want-to-refocus*))
%   sqrtR and valid come from refocuslite_companion
%   Returns a refocussed plane. No fancy pixel resizing.
%
% Run refocuslite_companion first to precalculate sqrtR and valid
% fplane   - is the fftshift(fft2(*whatever-you-want-to-refocus*))
% distance - is the distance the plane is to be propagated
% sqrtR    - comes from refocuslite_companion
% valid    - comes from refocuslite_companion
%
% This is propagation via a convolution with a phase map. Find details in
% Goodman's Fourier Optics.

phasemap = exp((1j).*2*pi.*distance.*sqrtR);
phasemap(valid)=0;
fnewplane = fplane.*phasemap; % Multiply in fourier space = convolution
fnewplane = ifftshift(fnewplane); % Undoes shift to centre
newplane = ifft2(fnewplane); 











