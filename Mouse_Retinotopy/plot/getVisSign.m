function [VSgnMap, VS_Thr] = getVisSign(map_hor, map_vert, reszfact)
% INPUTS:
%     map_horz - Map of horizontal retinotopy
%     map_vert - Map of vertical retinotopy
%     pixpermm = 1mm / pixSize(in mm) of the retinotopy images
% OUTPUTS:
%     VisSgnMap is a the visual field sign map
%     VS_Thr is threshold of VisSgnMap 

%%
pcopxsz  = 0.0065; % pixel size in PCO edge
% if nargin <4
%     pixpermm = 1/pcopxsz; %The size of each pixel is 0.0063 in PCO images
% end
if nargin <3
    reszfact = 1; % resize factor used for retinotopy maps
end
xdom = (0:size(map_hor,2)-1)*pcopxsz/reszfact;
ydom = (0:size(map_hor,1)-1)*pcopxsz/reszfact;

% mmperpix = 1/pixpermm;
% xdom = (0:size(map_horz,2)-1)*mmperpix;
% ydom = (0:size(map_horz,1)-1)*mmperpix;

%% create gaussian filter for smoothing
GF = fspecial('gaussian',size(map_hor), 10); % Gaussian Sigma = 10
GF = GF/sum(GF(:));

%% Plot absolute horz and vert retinotopic maps
% map_horz = ifft2( fft2(map_horz).*abs(fft2(GF)) );
% map_vert = ifft2( fft2(map_vert).*abs(fft2(GF)) );
figure, clf
subplot(2,2,1)
imagesc(xdom, ydom, map_hor); % [min(map_horz(:)) max(map_horz(:))]
axis image; colorbar; colormap jet;
title('Azm (deg)')

subplot(2,2,2)
imagesc(xdom, ydom, map_vert);% [min(map_vert(:)) max(map_vert(:))]
axis image; colorbar; colormap jet;
title('Elv (deg)')

%% Calculate visual field sign map VSgnMap
[hdx, hdy] = gradient(map_hor);
[vdx, vdy] = gradient(map_vert);

graddir_horz = atan2(hdy, hdx);
graddir_vert = atan2(vdy, vdx);

%[FX, FY] = gradient (F)
%   The 1st output FX is always the gradient along the 2nd dimension of F, going across columns. 
%   The 2nd output FY is always the gradient along the 1st dimension of F, going across rows. 
% vdiff = graddir_vert-graddir_horz;
vdiff   = exp(-1i*graddir_horz) .* exp(1i*graddir_vert); % Should be vert-horz, but see the comment above.
VSgnMap = sin(angle(vdiff)); % Visual field sign map
id      = isnan(VSgnMap);
id      = find(isnan(VSgnMap));
VSgnMap(id) = 0;

%% Spatial filtering (smoothing) the VSgnMap
% VSgnMap = imfilter( VSgnMap , GF, 'conv');
VSgnMap = ifft2( fft2(VSgnMap).*abs(fft2(GF)) );

%% Plotting smoothed VSgnMap
subplot(2,2,3); 
imagesc(xdom, ydom, VSgnMap); axis image%, [-1 1]
colorbar; colormap jet;
title('Visual Field Sign')

%% Apply threshold to the smoothed VSgnMap to create discrete patches
AbsVFS  = abs(VSgnMap);
Thresh  = 2 * std(VSgnMap(:));
VS_Thr  = (sign(AbsVFS - Thresh/2)+1)/2;  %threshold VSgnMap at (+-) 1.5 std

id      = find(VS_Thr);
imdum   = VS_Thr.*VSgnMap; 
imdum(id)= imdum(id)+1.1;

%% Plotting thresholded patches of VSgnMap
%% Plotting smoothed VSgnMap
subplot(2,2,4); 
% imagesc(xdom, ydom, imdum);
ploteccmapMA(imdum, [.1 2.1], xdom, ydom);
colorbar; 
axis image
StrThr = num2str(Thresh);
title(['Threshold: ' StrThr])

%% plot vfs threshold on the phase maps
subplot(221)
hold on
contour(xdom,ydom,VS_Thr,1,'k','linewidth',2)

subplot(222)
hold on
contour(xdom,ydom,VS_Thr,1,'k','linewidth',2)




%% Log
%{
- Mohammad 2015; wrote the file
%}



