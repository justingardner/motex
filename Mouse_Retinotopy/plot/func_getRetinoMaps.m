function [map_hor, map_vert, VFS, VFS_thr]=func_getRetinoMaps(CmplxMaps,visual_field,resFac)

%% Direct and Reverse maps for altitude and azimuth
trn = 1; % trial number
ang0 = angle(CmplxMaps{trn,1}{1});
ang2 = angle(CmplxMaps{trn,1}{2});
ang1 = angle(CmplxMaps{trn,2}{1});
ang3 = angle(CmplxMaps{trn,2}{2});

% phase elev/azim/sign maps (orig author: Ian Nauhaus)
% find delay as the angle between the vectors
delay_hor  = angle( exp(1i*ang0)+exp(1i*ang2) );
delay_vert = angle( exp(1i*ang1)+exp(1i*ang3) );

% make delay go from 0 to pi and 0 to pi, instead of 0 to pi and 0 to -pi;
% delay can't be negative. if the delay vector is in the bottom two 
% quadrants, it is assumed that it started at -180. the delay always pushes
% the vectors counter clockwise. (This is simply mod(val, pi), MA20150727)
delay_hor  = delay_hor  + (pi/2)*(1-sign(delay_hor));
delay_vert = delay_vert + (pi/2)*(1-sign(delay_vert));

% angle-correct the delay maps according to parameter delaythr
delaythr = 0;
angind = delay_hor < delaythr;
delay_hor(angind) = abs(delay_hor(angind) - pi);
angind = delay_vert < delaythr;
delay_vert(angind) = abs(delay_vert(angind) - pi);

% use delay vector to calculate retinotopy.
map_hor   = 0.5*(angle(exp(1i*(ang0-delay_hor)))  - angle(exp(1i*(ang2-delay_hor))));
map_vert  = 0.5*(angle(exp(1i*(ang1-delay_vert))) - angle(exp(1i*(ang3-delay_vert))));

% radians to degrees
map_hor   = (map_hor  /(2*pi)) * visual_field(1);
map_vert  = (map_vert /(2*pi)) * visual_field(2);

% Smooth the maps if needed and then get the visual field sign map
% GF = fspecial('gaussian',size(map_hor), 3);
% GF = GF/sum(GF(:));
% map_hor = ifft2( fft2(map_hor).*abs(fft2(GF)) );
% map_vert = ifft2( fft2(map_vert).*abs(fft2(GF)) );

[VFS, VFS_thr] = getVisSign(map_hor, map_vert, resFac);


%% Log
%{
- Mohammad 2015; wrote the file
%}