function h = sidetitle(ax, strtitle, side)
% SIDETITLE places a title on the right side of an axes object
% 
% h = sidetitle( s ) places the title s on the right side of the current axes
% object. returns handle h to text object.
%
% sidetitle( ax, s ) lets you specify the axes object
%
% sidetitle( ax, s, side) prints the title on the specified side (default:
% right). Valid options are 'right' (default), 'left', 'leftOutside'
%
% part of matteobox
%
% 2003-11 Matteo Carandini
% 2004-11 VB returns handle.
% 2008-06 LB adds the option to display title on the left side

if nargin < 3
    side = 'right';
end
if nargin < 2, strtitle = ax; ax = gca; end

oldax = gca;

axes(ax);

switch side
    case 'right'
        xpos = 1;
    case 'left'
        xpos = 0;
    otherwise
        error('<sidetitle> Unknown option %s', side);
end

h = text(xpos,0.5,strtitle,'units','norm','hori','left','vert','middle');

axes(oldax);