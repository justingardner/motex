function rgb = hv2rgb(h,v, cmap)
%HV2RGB Convert hue-value pairs to red-green-blue.
%
%   RGB = HV2RGB(H,V) converts the image H,V to the 
%   equivalent RGB image stored in the 3-D array (RGB).
%
%   As the hue H varies from 0 to 1, the resulting color varies from red,
%   through yellow, green, light blue, back to red. The "luminance" is kept
%   approximately constant. As the value V varies from 0 to 1, the
%   luminance increases.
%
%   RGB = HV2RGB(H,V, cmap) uses the colormap cmap for the conversion to
%   the 3-D RGB array. cmap contains a string with the desired colormap function, e.g. 'jet(128)'.
%
% 2005-12 Matteo Carandini
% 2006-09 Matteo Carandini simplified code - colors now consistent with HSV
% 2007-09 LB added argument for passing a string containing the appropriate
%           colormap

if nargin < 3
    cmap = [];
end

if ~all(size(h)==size(v))
    error('H and V must have the same dimensions');
end

if any( h(:)<0 | h(:)>1 | v(:)<0 | v(:)>1 )
    error('Values in H and V must be between 0 and 1');
end

v(isnan(v)) = 0;

if isempty(cmap)
    rgb = ind2rgb( gray2ind( mat2gray( h, [0 1] ), 128 ), hsv(128) );
else
    rgb = ind2rgb( gray2ind( mat2gray( h, [0 1] ), 128 ), eval(cmap));
end

rgb = rgb .*repmat( v, [ 1 1 3 ] );
    
%     
% [nx, ny] = size(h);
% 
% chan1 =  cos( h* 2*pi ); % the red-green channel
% chan2 =  sin( h* 2*pi ); % the blue-yellow channel
% 
% % This would be more similar to what we do with HSV:
% % chan1 =  cos( h* 2*pi  ); % the red-green channel
% % chan2 =  sin( h* 2*pi + pi ); % the blue-yellow channel
% 
% rgby = zeros(4,nx,ny);
% rgby(1,:,:) = max(0, chan1).^2;
% rgby(2,:,:) = max(0,-chan1).^2;
% rgby(3,:,:) = max(0, chan2).^2;
% rgby(4,:,:) = max(0,-chan2).^2;
% 
% rgb = zeros(nx,ny,3);
% 
% % rgb(:,:,1) = sum( repmat([ 1 0 0.4 1 ]',[1 nx ny]) .* rgby, 1 );
% % rgb(:,:,2) = sum( repmat([ 0 1 0.4 1 ]',[1 nx ny]) .* rgby, 1 );
% % rgb(:,:,3) = sum( repmat([ 0 0  1  0 ]',[1 nx ny]) .* rgby, 1 );
% 
% rgb(:,:,1) = sum( repmat([ 1 0  0  1 ]',[1 nx ny]) .* rgby, 1 );
% rgb(:,:,2) = sum( repmat([ 0 1 0.5 1 ]',[1 nx ny]) .* rgby, 1 );
% rgb(:,:,3) = sum( repmat([ 0 0  1  0 ]',[1 nx ny]) .* rgby, 1 );
% 
% rgb = rgb .* repmat(v,[1 1 3]);

return

%---------------------------------------------------
%       code that went into designing the function
%---------------------------------------------------

n = 64;

thetas = linspace(0,2*pi,n+1);
thetas = thetas(1:end-1);

chan1 =  cos( thetas );
chan2 =  sin( thetas );

figure; 
plot(thetas, chan1, 'k-' ); hold on
plot(thetas, chan2, 'k--' );

rgby(1,:) = max(0, chan1).^2;
rgby(2,:) = max(0,-chan1).^2;
rgby(3,:) = max(0, chan2).^2;
rgby(4,:) = max(0,-chan2).^2;

figure; 
plot( thetas, rgby(1,:), 'r' ); hold on
plot( thetas, rgby(2,:), 'g' );
plot( thetas, rgby(3,:), 'c' );
plot( thetas, rgby(4,:), 'y' );

rgb = [ 1 0 0; 0 1 0; 0.4 0.4 1; 1 1 0 ]' * rgby;

lum_contribs = [ 1 1 0.3 ];

lum = lum_contribs * rgb

figure; 
plot( thetas, rgb(1,:), 'r' ); hold on
plot( thetas, rgb(2,:), 'g' );
plot( thetas, rgb(3,:), 'b' );
plot( thetas, lum,      'k' );



