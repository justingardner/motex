function [jointmap, amplim, legendmap] = bpMakeJointMap(ampmap, angmap, amplim, mycolormap, rectLegend)
% bpMakeJointMap combines an intensity map and an angle map
% 
% jointmap = bpMakeJointMap(ampmap,angmap)
% takes as input a matrix of amplitude values ampmap (must be >=0), and a
% matrix of phases angmap (must be in [-pi,pi]) and returns a jointmap,
% which is a 3D map, one for each gun (R, G, B). You can show it by using
% imshow.
%
% [jointmap, amplim] = bpMakeJointMap(ampmap,angmap)
% returns also amplim, the clipping value used in the ampmap. 
% 
% [jointmap, amplim, legendmap] = bpMakeJointMap(ampmap,angmap)
% returns also legendmap, a sort of pinwheel legend for the joint map.
% Again, use imshow to see it. 
% 
% bpMakeJointMap(ampmap, angmap, amplim) allows you to specify amplim
% (DEFAULT: the maximum value of ampmap).
% 
% bpMakeJointMap(ampmap, angmap, [], mycolmap) allows you to specify the
% desired colormap function, e.g. 'jet(128)' (DEFAULT: hsv(128)).
%
% bpMakeJointMap(ampmap, angmap, [], [], rectLegend) outputs a rectangular
% instead of a circular legend, if rectLegend is set to 1.
% Useful when using a non-circular colormap.
%
% EXAMPLE:
%
% foo = 2*(rand(100,200)-0.5) + i*2*(rand(100,200)-0.5);
% [jointmap, amplim, legendmap] = bpMakeJointMap( abs(foo) , angle(foo) );
% figure; 
% subplot(2,1,1); imshow(jointmap);
% subplot(2,1,2); imshow(legendmap); title(sprintf('rad = %2.2f',amplim));
 
% Matteo Carandini
% 2004-10 MC added amplim
% 2005-12 MC removed luminance confound (replaced colormap hack with Matlab function hsv2rgb)
% 2005-12 MC improved by calling the new function hv2rgb rather than hsv2rgb
% 2006-03 MC added as output legendmap.
% 2006-09 MC inverted the whole thing!! Now absence of activity is white.
% 2007-12 LB added as input the colormap
% 2008-01 LB added the option for a rectangular legendmap

if any(size(ampmap)~=size(angmap))
    error('Must have same sizes');
end

if nargin < 5
    rectLegend = 0;
end

if nargin < 4
    mycolormap = [];
end

if nargin < 3
    amplim = max(max(ampmap));
end

if any( angmap(:)<-pi | angmap(:)>pi )
     error('The angle map must be in radians between -pi and pi');
end

%-------------------------------------------

[nx, ny] = size(ampmap);

nanmap = isnan(ampmap);

ampmap(nanmap) = 0;

ampmap = ampmap/amplim;
ampmap = min(1,ampmap); % ensures that we don't scale by more than 1
ampmap = max(0,ampmap);

%new: rotate by 180 deg
if rectLegend == 0
    angmap = angle( exp(i*angmap)*exp(i*pi) );
end


angmap = (angmap+pi)/(2*pi);  % now it goes from 0 to 1

jointmap = hv2rgb( angmap, ampmap, mycolormap );

% to fix a strange problem that should not occur
jointmap((jointmap>1)) = 1;

% new: invert
if rectLegend == 0
    jointmap = 1 - jointmap;
end

if nargout == 3
    if ~rectLegend % circular legend
        n = 50;
        [xx, yy] = meshgrid(-n:n,-n:n);
        AA = sqrt(xx.^2+yy.^2)/n; AA(AA>1)=0;
        PP = angle(xx+i*yy);
    else
        AA = repmat((1:-0.02:0)', [1 51]);
        PP = repmat(-pi:pi/25:pi, [51, 1]); 
    end
    legendmap = bpMakeJointMap( AA*amplim, PP, amplim, mycolormap, rectLegend);
end


