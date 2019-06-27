function map = hv(m)
%HV    "Color opponent" color map
%   HV(M) returns an M-by-3 matrix containing an RGBY colormap.
%   HV, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   An HV colormap has approximate constant "luminance" and varies from
%   red to yellow to green to light blue, and finally back to red. The map
%   is particularly useful for displaying periodic functions.  
%
% 2005-12 Matteo Carandini

if nargin < 1
    m = size(get(gcf,'colormap'),1);
end

h = (0:m-1)'/max(m,1);
map = squeeze(hv2rgb( h, ones(m,1)));

% 
% figure; 
% plot(map(:,1),'r'); hold on
% plot(map(:,2),'g'); hold on
% plot(map(:,3),'b'); hold on
