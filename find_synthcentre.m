function [Cx,Cy,R] = find_synthcentre(fim,DC_w,LP_n,thresh_m,flag_plot)

%%%%%%%%%%%%%%%%%%
% This function finds the centre of a cross term in the Fourier transform
% of a hologram.
% fim       - the Fourier transform of a hologram
% DC_w      - the radius of the DC terms - better to overestimate this
% LP_n      - the width of the low pass filter
% thresh_m  - the threshold is the mean*thresh_m
% flag_plot - 1 for a final plot, 0 for no plot

% This is fairly crude and I've only tested it for one hologram but good
% parameters seem to be:
% DC_w = 200;
% LP_n = 5;
% thresh_m = 1;

% The hologram I used to develop it was:
% 'R:\Lisong Yang\AndyFunctions\DF-ref-holo-24-03-2015-01.tif'

% The function works by first removing the really bright stuff - the DC
% terms and the horizontal and vertical bright stripes.
% The data is then low-pass filtered
% Then a threshold is calculated and applied.
% The points from the top of the image are retained - a single cross term
% A convex hull around these points is calculated.
% Only points from above the equator of the convex hull are retained - the
% DC masking takes a bite out of the cross term.
% A circle fitting routine (that I am fairly proud of) fits a circle to
% these points and gives the coordinates of the centre.
%%%%%%%%%%%%%%%%%%

if flag_plot
    fim_save = fim;
end

s = size(fim);
[mx,my] = meshgrid(1:s(2),1:s(1));
mx_1 = mx - s(2)/2 - 1;
my_1 = my - s(1)/2 - 1;

% Remove very bright bits - DC and the stripes at fx = 0 and fy = 0
DC_mask = mx_1.^2 + my_1.^2 < DC_w.^2;
DC_mask = (DC_mask | abs(mx_1) < 5) | abs(my_1) < 5;
thresh = mean(abs(fim(:))).*thresh_m;
fim(DC_mask) = 0;
%image(abs(fim)/1e2);axis equal

%%
% Low pass filter
kern = pascal(LP_n,1);
kern = abs(kern(end,:));
kern = kern.'*kern;
kern = kern./sum(kern(:)); % Roughly gaussian kernal that integrates to 1

LP = conv2(abs(fim),kern,'same');   % Convolve

%%
% Outlier removal - this isn't really necessary
TOP = my_1 < 0;
BTM = my_1 > 0;
valid_TOP = (LP>thresh) & TOP;
cent_TOP = mean([mx(valid_TOP),my(valid_TOP)],1);
outliermask_TOP = sqrt( (mx - cent_TOP(1)).^2 + (my - cent_TOP(2)).^2 ) < DC_w*2/3;

valid_BTM = (LP>thresh) & BTM;
cent_BTM = mean([mx(valid_BTM),my(valid_BTM)],1);
outliermask_BTM = sqrt( (mx - cent_BTM(1)).^2 + (my - cent_BTM(2)).^2 ) < DC_w*2/3;

outliermask = outliermask_TOP | outliermask_BTM;


%%
filtered = LP>thresh & outliermask;

TOP = my_1 < 0;

XT = mx(filtered & TOP);
YT = my(filtered & TOP);
KT = convhull(XT,YT);

midy = mean(YT);

XT = XT(KT);
YT = YT(KT);

valid = YT < midy;
XT = XT(valid);
YT = YT(valid);

% Fit a circle to the top of the hull
x = reshape(XT,1,length(XT));
y = reshape(YT*s(2)/s(1),1,length(YT));

F1 = 2*x;
F2 = 2*y;
F3 = ones(size(y));

F = [F1',F2',F3'];
a = F\(x.^2+y.^2)';

Cx = a(1);
Cy = a(2);
R = sqrt(a(3) + (Cx.^2 + Cy.^2));

tta = linspace(0,2*pi,100);
rho = R*ones(size(tta));
[x,y] = pol2cart(tta,rho);
x = x+Cx;
y = y+Cy;

Cy = Cy*s(1)/s(2);

if flag_plot
image(mx(1,:),my(:,1),abs(fim_save)/1e2)
hold on
plot(x,y*s(1)/s(2),'-g',Cx,Cy,'og'); axis equal tight
hold off
end
