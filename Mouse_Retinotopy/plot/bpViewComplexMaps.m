function [figs, jointImg, max_amp] = bpViewComplexMaps( ComplexMaps, RowTagNumbers, ColTagNumbers, roi, yscale, plotflag )
% bpViewComplexMaps plots amplitudes and phases 
%
% [figs, jointImg, max_amp] = bpViewComplexMaps( ComplexMaps )
%
% bpViewComplexMaps( ComplexMaps, RowTagNumbers) lets you specify an array
% of tags, one for each row. They must be numbers. (DEFAULT: [], which means
% [1:nrows]).
%
% bpViewComplexMaps( ComplexMaps, RowTag Numbers, ColTagNumbers) lets you specify an array
% of tags, one for each column. They must be numbers. (DEFAULT: [], which means
% [1:ncols]).
%
% bpViewComplexMaps( ComplexMaps, [], [], roi) lets you specify the region
% of interest ROI (DEFAULT: [])
%
% bpViewComplexMaps( ComplexMaps, [], [], [], yscale ) lets you specify the maximum amplitude.
% Can also be set to 'matchy' to scale them all to the 99th percentile, or
% to 'byrow' (or []) to do so row by row (which is the default).
%
% bpViewComplexMaps( ComplexMaps, [], [], [], [], 'complex') shows only the
% joint map (DEFAULT: 'all')
% bpViewComplexMaps( ComplexMaps, [], [], [], [], 'amp') shows only the
% amplitude map (DEFAULT: 'all')
% bpViewComplexMaps( ComplexMaps, [], [], [], [], 'phase') shows only the
% phase map (DEFAULT: 'all')
% 
% 
% 2005-09 MC derived from PlotAmpPhaseMaps
% 2005-12 MC substituted HSV colormap with HV colormap, fixed a bug with
% minima in amplitude plots.
% 2006-03 AB added jointImg output
% 2006-09 MC from bpViewAbsPhase
% 2007-11 LB added option to specify which map to look at

%% Parse the inputs

if nargin < 6
    plotflag = 'all';
end

if nargin < 5
    yscale = 'byrow';
end

if ~iscell(ComplexMaps) && ndims(ComplexMaps)==3
    [ncols, nrows,ns] = size(ComplexMaps);
    foo = cell(ns,1);
    for is = 1:ns
        foo{is} = ComplexMaps(:,:,is);
    end
    ComplexMaps = foo;
end

if ~iscell(ComplexMaps) && ndims(ComplexMaps)==2
    ComplexMaps = {ComplexMaps};
end

if nargin < 4 || isempty(roi)
    roi = ones(size(ComplexMaps{1}));
end

if nargin < 3
    ColTagNumbers = [];
end

if nargin < 2
    RowTagNumbers = [];
end

%% Basic definitions
    
[ncols, nrows] = size(ComplexMaps);
nImg = ncols * nrows;

if isempty(RowTagNumbers)
    RowTagNumbers = 1:nrows;
end
if isempty(ColTagNumbers)
    ColTagNumbers = 1:ncols;
end

figs = [];

%% Crop according to ROI 

rows = any(roi');
cols = any(roi);

for iImg = 1:nImg
    ComplexMaps{iImg} = ComplexMaps{iImg}(rows,cols);
end

newroi = roi(rows,cols);

%% Rotate if appropriate

% [nx, ny] = size(ComplexMaps{1});
% 
% RotateAngle = 0;
% 
% if (nrows>ncols)&&(nx>ny), RotateAngle = -90; end
% if (nrows<ncols)&&(nx<ny), RotateAngle = -90; end
% 
% if RotateAngle~=0
%     for iImg = 1:nImg
%         ComplexMaps{iImg} = imrotate( ComplexMaps{iImg}, RotateAngle );
%     end
%     newroi = imrotate(newroi, RotateAngle);
% end

%% Extract amplitude and phase

sigAmp   = {};
sigPhase = {};
for irow = 1:nrows
    for icol = 1:ncols
        sigAmp  {icol,irow} =   abs(ComplexMaps{icol,irow});
        sigPhase{icol,irow} = angle(ComplexMaps{icol,irow});
    end
end

%% 

max_amp = zeros(1,nrows);
if isempty(yscale) || ischar(yscale)
    for irow = 1:nrows
        allAmps = [sigAmp{:,irow}];
        max_amp(irow) = prctile(allAmps(:),99);
    end
    if ischar(yscale) && strcmp(yscale,'matchy')
        max_amp(:) = max(max_amp);
    end
else
    max_amp(:) = yscale;
end

% this is for the amplitude plot
min_amp = 0; % don't go above 0
for irow = 1:nrows
    allAmps = [sigAmp{:,irow}];
    min_amp(irow) = min(allAmps(:));
end
if ischar(yscale) && strcmp(yscale,'matchy')
    min_amp(:) = min(min_amp);
end
  
jointlegendmap = {};
for irow = 1:nrows
    for icol = 1:ncols
        jointmap{icol,irow} = bpMakeJointMap(sigAmp{icol,irow},sigPhase{icol,irow},max_amp(irow));

        n = 50;
        [xx, yy] = meshgrid(-n:n,-n:n);
        AA = sqrt(xx.^2+yy.^2)/n; AA(AA>1)=0;
        PP = angle(xx+i*yy);
        jointlegendmap{irow} = bpMakeJointMap( AA*max_amp(irow), PP, max_amp(irow));

    end
end

%% The joint map 
if strcmp(plotflag, 'complex') || strcmp(plotflag, 'all')
    JointFig = figure;
    figs = [figs JointFig];

    jointImg = {};
    roi3d = repmat( newroi, [1, 1, 3] );

    OutofROIGrayLevel = 1;

    ax_joint = zeros(ncols,nrows);
    for icol = 1:ncols
        for irow = 1:nrows
            ax_joint(icol,irow) = gridplot(nrows, ncols+1, irow, icol+1);
            myimg = jointmap{icol,irow} .* roi3d + OutofROIGrayLevel*(1-roi3d);
            im = imshow( myimg );
            set(im,'buttondownfcn','bpPlaceCross(''w+'')');
            jointImg{irow,icol} = myimg;
        end
    end

    if strcmp(yscale,'matchy')
        rowlist = 1;
    else
        rowlist = 1:nrows;
    end

    for irow = rowlist
        gridplot(nrows, ncols+1, irow, 1);
        imgobj = imshow( jointlegendmap{irow} );
        axis xy
        xlabel('-pi/2');
        ylabel('pi');
        title(sprintf('%2.3g',max_amp(irow)));
        set(imgobj,'userdata','noscalebar'); % this way there will not be a scale bar on it
    end

    % scalebar(cmap_ax, 'Phase','ylim',[-pi pi],'ytick',pi*[-1:1/2:1],'yticklabel',{'-pi' '-pi/2' '0' 'pi/2' 'pi'});

    for irow = 1:nrows
        axes( ax_joint(end,irow) );
        sidetitle(num2str(RowTagNumbers(irow)));
    end
    for icol = 1:ncols
        axes( ax_joint(icol,1) );
        title(num2str(ColTagNumbers(icol)));
    end

    set( ax_joint,'dataaspectratio',[1 1 1]);
    set( ax_joint(2:end),'xticklabel',[],'yticklabel',[]);

    set(gcf,'pointer','crosshair');
end
%% The amplitude map

if strcmp(plotflag, 'amp') || strcmp(plotflag, 'all')
    AmpFig = figure;
    figs = [figs AmpFig];

    ax_amp = zeros(ncols,nrows);
    for icol = 1:ncols
        for irow = 1:nrows
            ax_amp(icol,irow) = gridplot(nrows, ncols+1, irow, icol+1);
            myimg = sigAmp{icol,irow};
            myimg(~newroi(:)) = 0;
            im = imagesc( myimg );
            set(im,'buttondownfcn','bpPlaceCross(''r+'')');
            caxis([ min_amp(irow) max_amp(irow)]);
        end
    end
    colormap gray
    % colormap(1-gray(128))
    for irow = 1:nrows
        cmap_ax = gridplot(nrows, ncols+1, irow, 1);
        scalebar(cmap_ax, 'Amplitude', 'ylim', [ min_amp(irow) max_amp(irow)]);
    end
    for irow = 1:nrows
        axes( ax_amp(end,irow) );
        sidetitle(num2str(RowTagNumbers(irow)));
    end
    for icol = 1:ncols
        axes( ax_amp(icol,1) );
        title(num2str(ColTagNumbers(icol)));
    end

    set(ax_amp,'dataaspectratio',[1 1 1],'xtick',[],'ytick',[],'box','off','xcolor', 'w', 'ycolor', 'w');
end
 
%% The phase map

if strcmp(plotflag, 'phase') || strcmp(plotflag, 'all')
    PhaseFig = figure;
    figs = [figs PhaseFig];

    ax_phase = zeros(ncols,nrows);
    for icol = 1:ncols
        for irow = 1:nrows
            ax_phase(icol,irow) = gridplot(nrows, ncols+1, irow, icol+1);
            my_img = sigPhase{icol,irow};
            my_img( ~newroi(:) ) = NaN;
            im = imagesc(my_img);
            set(im,'buttondownfcn','bpPlaceCross(''w+'')');
            caxis([-pi pi]);
        end
    end
    colormap(hv);
    cmap_ax = gridplot(nrows, ncols+1, 1, 1);
    scalebar(cmap_ax, 'Phase','ylim',[-pi pi],'ytick',pi*(-1:1/2:1),'yticklabel',{'-pi' '-pi/2' '0' 'pi/2' 'pi'});
    for irow = 1:nrows
        axes( ax_phase(end,irow) );
        sidetitle(num2str(RowTagNumbers(irow)));
    end
    for icol = 1:ncols
        axes( ax_phase(icol,1) );
        title(num2str(ColTagNumbers(icol)));
    end

    set(ax_phase,'dataaspectratio',[1 1 1]);
    set(ax_phase(2:end),'xticklabel',[],'yticklabel',[]);
end