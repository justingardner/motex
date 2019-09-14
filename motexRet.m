% motexRet.m
%
%      usage: motexRet(v)
%         by: justin gardner
%       date: 09/10/19
%    purpose: To be run on mrTools mouse retinotopy data sets that have been run through motex2mrtools
%
%             cd ~/data/M190724_s5r1_2;
%             v = newView;
%             motexRet(v);
%
function retval = motexRet(v,varargin)

% check arguments
if nargin < 1
  help motexRet
  return
end

% parse arguments
%getArgs(varargin,{'preProcessedDirs',{'/Volumes/tex/IMAGING/WIDEFIELD/'},'resFac=1'});
getArgs(varargin,{'preProcessedDirs',{'~/data/motex/raw/retinotopy'},'resFac=1','corAnalName=corAnal','saveCorAnalOverlays=0','savePreProcessedOverlays=0','makeROIs=0'});

% check that the view sturcture has a retinotopy scan
[tf v] = checkView(v);
if ~tf,return,end

% check retinotopy scans in averages
[v retinotopyInfo] = checkRetinotopyScans(v,corAnalName);
if ~any(retinotopyInfo.isret),return,end

% now try to load pre-processed retinotopy information
retinotopyInfo = loadPreProcessedRetinotopy(retinotopyInfo,preProcessedDirs);

% compute the maps from the pre-processed info
if savePreProcessedOverlays || makeROIs
  retinotopyInfo = makePreProcessedMaps(retinotopyInfo,resFac);
end

% now run retinotopy through the automatic visual field sign code
if saveCorAnalOverlays
  if ~isempty(retinotopyInfo.corAnal) && ~isempty(retinotopyInfo.preProcessed)
    % get visual_field which I do not yet know how to compute from pre_processed file
    visual_field = retinotopyInfo.preProcessed{1}.visual_field;
    % compute visual field sign
    [v retinotopyInfo.corAnal] = computeVisualFieldSign(v,retinotopyInfo.corAnal,visual_field,resFac);
    % load map overlays into corAnal analysis
    if ~isempty(retinotopyInfo.corAnal)
      % add map overlay computed for h/v field map from correlation analysis
      addMapOverlays(v,retinotopyInfo.corAnal);
    end
  end
end

% make pre-processed maps and save into MLR
if savePreProcessedOverlays && ~isempty(retinotopyInfo.preProcessed)
  % save the first pre-processed maps to MLR
  addMapOverlays(v,retinotopyInfo.preProcessed{1});
end

if makeROIs
  % convert into rois for MLR
  rois = makeMLRRois(v,retinotopyInfo);
  % and save the rois
  saveROI(v,rois);
end

% now figure out what the stimImages should be for each
% corAnal scan so that we can run pRF analysis
addStimImages(v,retinotopyInfo);

keyboard
%%%%%%%%%%%%%%%%%%%%%%%
%    addStimImages    %
%%%%%%%%%%%%%%%%%%%%%%%
function addStimImages(v,retinotopyInfo)

% get the stim images
stimImages = motexGetRetinotopyStimImages;
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%
%    addMapOverlays    %
%%%%%%%%%%%%%%%%%%%%%%%%
function v = addMapOverlays(v,maps)

% make hMap overlay
o.name = 'hMap';
o.groupName = viewGet(v,'GroupName');
o.params = viewGet(v,'overlayParams');
o.range = [min(maps.hMap(:)) max(maps.hMap(:))];
for iScan = 1:viewGet(v,'nScans')
  o.data{iScan} = maps.hMap;
end
[tf o] = isoverlay(o);
% and add to analysis
v = viewSet(v,'newOverlay',o);

% make vMap overlay
o = [];
o.name = 'vMap';
o.groupName = viewGet(v,'GroupName');
o.params = viewGet(v,'overlayParams');
o.range = [min(maps.vMap(:)) max(maps.vMap(:))];
for iScan = 1:viewGet(v,'nScans')
  o.data{iScan} = maps.vMap;
end
[tf o] = isoverlay(o);
v = viewSet(v,'newOverlay',o);

% make visualFieldSign overlay
o = [];
o.name = 'visualFieldSign';
o.groupName = viewGet(v,'GroupName');
o.params = viewGet(v,'overlayParams');
o.range = [min(maps.visualFieldSign(:)) max(maps.visualFieldSign(:))];
for iScan = 1:viewGet(v,'nScans')
  o.data{iScan} = maps.visualFieldSign;
end
[tf o] = isoverlay(o);
v = viewSet(v,'newOverlay',o);

% make thresholdMap
o = [];
o.name = 'thresholdMap';
o.groupName = viewGet(v,'GroupName');
o.params = viewGet(v,'overlayParams');
o.range = [min(maps.thresholdMap(:)) max(maps.thresholdMap(:))];
for iScan = 1:viewGet(v,'nScans')
  o.data{iScan} = maps.thresholdMap;
end
[tf o] = isoverlay(o);
v = viewSet(v,'newOverlay',o);

% save the analysis
saveAnalysis(v,'corAnal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeVIsualFieldSign    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v corAnal] = computeVisualFieldSign(v,corAnal,visual_field,resFac)

% make sure we have all the scans we need
missingDirs = setdiff([0 90 180 270],corAnal.dirs);
if ~isempty(missingDirs)
  disp(sprintf('(motexRet:computeVisualFieldSign) Need to have corAnal run on directions: %s',num2str(missingDirs)));
  return
end

% convert overlays into complex numbers
corAnal.complexMap{1}{1} = makeComplexMap(first(find(corAnal.dirs==0)),corAnal);
corAnal.complexMap{1}{2} = makeComplexMap(first(find(corAnal.dirs==180)),corAnal);
corAnal.complexMap{2}{1} = makeComplexMap(first(find(corAnal.dirs==90)),corAnal);
corAnal.complexMap{2}{2} = makeComplexMap(first(find(corAnal.dirs==270)),corAnal);

% filter the maps
if 0
  corAnal.mapFilterSigma = 10;
  for i = 1:2
    for j = 1:2
      corAnal.complexMap{i}{j} = spatialGaussianFilter(corAnal.complexMap{i}{j},corAnal.mapFilterSigma);
    end
  end
end

% check for functions form Moha
if isempty(which('func_getRetinoMaps'))
  disp(sprintf('(motexRet:computeVisualFieldSign) Colud not find function func_getRetinoMaps'));
  return
end

% run the retino maps function
[corAnal.hMap, corAnal.vMap, corAnal.visualFieldSign, corAnal.thresholdMap] = func_getRetinoMaps(corAnal.complexMap,visual_field,resFac);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    spatialGaussianFilter  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = spatialGaussianFilter(map,filterSigma)

% make x and y corrdinates in pixel dimensions with
% 0,0 in the center
[ylen xlen] = size(map);
filterX = (0:(xlen-1)) - (xlen-1)/2;
filterY = (0:(ylen-1)) - (ylen-1)/2;
[filterX filterY] = meshgrid(filterX,filterY);

% compute gaussian of appropriate filterSigma (in pixels)
% on this grid to make gaussian filter
filter = exp(-(filterX.^2+filterY.^2)./filterSigma^2);

% now compute convolution in the frequency domain
% FIX, FIX, FIX, this does not seem to work as expected - why not?
%map = ifft2(fft2(map).*fft2(filter));
% this works, but for some reason blocks off a corner of the image and is slow
map = conv2(map,filter,'same');

%%%%%%%%%%%%%%%%%%%%%%%%
%    makeComplexMap    %
%%%%%%%%%%%%%%%%%%%%%%%%
function complexMap = makeComplexMap(scanNum,corAnal)

% convert amplitude and phase back to a complex number
[a b] = pol2cart(squeeze(corAnal.ph(scanNum,:,:)),squeeze(corAnal.amp(scanNum,:,:)));
complexMap = complex(a,b);

%%%%%%%%%%%%%%%%%%%%%
%    makeMLRRois    %
%%%%%%%%%%%%%%%%%%%%%
function rois = makeMLRRois(v,retinotopyInfo)

rois = {};
% cycle through all pre-processed files
for iPreProcessed = 1:length(retinotopyInfo.preProcessed)
  % shortcut to structure
  p = retinotopyInfo.preProcessed{iPreProcessed};
  % set color table for rois
  colors = hsv(p.regions.n);
  for iROI = 1:p.regions.n
    % set name
    if length(retinotopyInfo.preProcessed) > 1
      rois{end+1}.name = sprintf('%s_%i_%i Area %i',p.filename,p.sessioNum,p.runNum,iROI);
    else
      rois{end+1}.name = sprintf('Area %i',iROI);
    end
    % set other info
    rois{end}.voxelSize = viewGet(v,'scanVoxelSize');
    rois{end}.xform = viewGet(v,'scansform');
    % set coords
    rois{end}.coords = p.regions.coords{iROI};
    rois{end}.coords(3,:) = 1;
    rois{end}.color = colors(iROI,:);
    % make into an roi
    [tf rois{end}] = isroi(rois{end});
  end
end

%%%%%%%%%%%%%%%%%
%    saveROI    %
%%%%%%%%%%%%%%%%%
function saveROI(v,roi)

% get roi dir
roiDIR = viewGet(v,'roiDir');
if ~isdir(roiDIR),mkdir(roiDIR);end

% make into cell array
roi = cellArray(roi);

% for each roi, go and save
for iROI = 1:length(roi)
  % get the save name
  savename = fixBadChars(roi{iROI}.name,[],{'.','_'});
  eval(sprintf('%s = roi{iROI};',savename));
  save(fullfile(roiDIR,savename),savename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makePreProcessedMaps    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retinotopyInfo = makePreProcessedMaps(retinotopyInfo,resFac)

% nothing to do if nothing found
if isempty(retinotopyInfo.preProcessed)
  disp(sprintf('(motexRet:makePrePRocessedMaps) No pre-processed retinotopy information to run'));
  return
end

% check for functions form Moha
if isempty(which('func_getRetinoMaps'))
  disp(sprintf('(motexRet:makePrePRocessedMaps) Colud not find function func_getRetinoMaps'));
  return
end

% now run it to get maps
for iPreProcessed = 1:length(retinotopyInfo.preProcessed)
  % shortcut to structure
  p = retinotopyInfo.preProcessed{iPreProcessed};
  % run the retino maps function
  [p.hMap, p.vMap, p.visualFieldSign, p.thresholdMap] = func_getRetinoMaps(p.CmplxMaps,p.visual_field,resFac);
  % find connected regions
  connComp = bwconncomp(p.thresholdMap);
  p.regions.n = connComp.NumObjects;
  % now convert these into x,y positions
  for iConnComp = 1:connComp.NumObjects
    [x y] = ind2sub(connComp.ImageSize,connComp.PixelIdxList{iConnComp});
    p.regions.nCoords(iConnComp) = length(x);
    % note that x and y are flipped for MLR
    p.regions.coords{iConnComp} = [y(:)';x(:)'];
  end
  % reorder by size (since V1 is probably the largest connected region
  [p.regions.nCoords sortorder] = sort(p.regions.nCoords,2,'descend');
  p.regions.coords = {p.regions.coords{sortorder}};
  % pack structure back
  retinotopyInfo.preProcessed{iPreProcessed} = p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadPreProcessedRetinotopy    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retinotopyInfo = loadPreProcessedRetinotopy(retinotopyInfo,preProcessedDirs)

% make sure preProcessedDirs is a cell array
preProcessedDirs = cellArray(preProcessedDirs);

% check through preProcessedDirs for a pre-processed-retinotopy
disppercent(-inf,sprintf('(motex:loadPreProcessedRetinotopy) Loading pre-processed retinotopy'));
retinotopyInfo.preProcessed = {};
for iScan = 1:retinotopyInfo.n
  % default to not found
  retinotopyInfo.preProcessedDir{iScan} = '';
  retinotopyInfo.overlayFigname{iScan} = '';
  retinotopyInfo.retinotopyFigname{iScan} = '';
  retinotopyInfo.preProcessedNum(iScan) = nan;
  retinotopyInfo.preProcessedFilename{iScan} = '';
  % no go search for it through all the possible paths
  for iDir = 1:length(preProcessedDirs)
    % look for pre-processed dir
    preProcessedFilename = fullfile(preProcessedDirs{iDir},retinotopyInfo.filename{iScan},'Analyzed',num2str(retinotopyInfo.sessionNum(iScan)),'AnalyzedRet.mat');

    % keep it if it exists
    if isfile(preProcessedFilename)
      retinotopyInfo.preProcessedFilename{iScan} = preProcessedFilename;
    end
    % now look for overlay figures
    overlayFigname = fullfile(preProcessedDirs{iDir},retinotopyInfo.filename{iScan},num2str(retinotopyInfo.sessionNum(iScan)),'overlay.fig');
    if isfile(overlayFigname)
      retinotopyInfo.overlayFigname{iScan} = overlayFigname;
    end
    % now look for retinotopy figures
    retinotopyFigname = fullfile(preProcessedDirs{iDir},retinotopyInfo.filename{iScan},num2str(retinotopyInfo.sessionNum(iScan)),'retinotopy.fig');
    if isfile(retinotopyFigname)
      retinotopyInfo.retinotopyFigname{iScan} = retinotopyFigname;
    end
  end
  % tell user if not found
  if isempty(retinotopyInfo.preProcessedFilename{iScan})
    disp(sprintf('(motexRet:loadPrePRocessedRetinotopy) No pre-processed retinotopy found for %s sessionNum: %i runNum: %i',retinotopyInfo.filename{iScan},retinotopyInfo.sessionNum(iScan),retinotopyInfo.runNum(iScan)));
  else
    % load it up if it is a new one
    [isLoaded loadNum] = ismember(retinotopyInfo.preProcessedFilename{iScan},{retinotopyInfo.preProcessedFilename{1:(iScan-1)}});
    if ~isLoaded
      % load
      retinotopyInfo.preProcessed{end+1} = load(retinotopyInfo.preProcessedFilename{iScan});
      % set info about this retinotopy
      retinotopyInfo.preProcessed{end}.filename = retinotopyInfo.filename{iScan};
      retinotopyInfo.preProcessed{end}.sessionNum = retinotopyInfo.sessionNum(iScan);
      retinotopyInfo.preProcessed{end}.runNum = retinotopyInfo.runNum(iScan);
      % keep a poitner to which one is for this scan
      retinotopyInfo.preProcessedNum(iScan) = length(retinotopyInfo.preProcessed);
    else
      % make the pointer point to the appropritate location
      retinotopyInfo.preProcessedNum(iScan) = retinotopyInfo.preProcessedNum(loadNum);
    end
  end
  disppercent(iScan/retinotopyInfo.n);
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkRetinotopyScans    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v retinotopyInfo] = checkRetinotopyScans(v,corAnalName)

% set groupName
retinotopyInfo.groupName = viewGet(v,'groupName');

% check that we have retinotoyp scans in the current 
% group (which was set by checkVIew)
disppercent(-inf,sprintf('(motexRet:checkRetinotopyScans) Checking for retinotopy stimfiles'));
nScans = viewGet(v,'nScans');
for iScan = 1:nScans
  % first, assume that it is a retinotopy scan
  retinotopyInfo.isret(iScan) = true;
  retinotopyInfo.dir(iScan) = nan;
  retinotopyInfo.contrast(iScan) = nan;
  retinotopyInfo.runInfo{iScan} = [];
  retinotopyInfo.filename{iScan} = '';
  % get the stimfile
  stimfile = viewGet(v,'stimfile',iScan);
  for iStimfile = 1:length(stimfile)
    % check for runInfo
    if ~isfield(stimfile{iStimfile},'runInfo')
      disp(sprintf('(motexRet:checkRetinotopyScans) Missing runInfo field for scan %i',iScan));
      retinotopyInfo.isret(iScan) = false;
      break;
    end
    % check for retinotopy filed
    if ~isfield(stimfile{iStimfile}.runInfo,'retinotopy')
      disp(sprintf('(motexRet:checkRetinotopyScans) Missing retinotopy field in runInfo for scan %i',iScan));
      retinotopyInfo.isret(iScan) = false;
      break;
    end
    % check that they all have the same direction of movement
    % by looking first to see which run this started out as
    filenameSplit = strsplit(stripext(getLastDir(stimfile{iStimfile}.filename)),'_');
    fileNum = str2num(filenameSplit{end});
    % get session and run while we are at it
    retinotopyInfo.filename{iScan} = getLastDir(fileparts(fileparts(stimfile{iStimfile}.runInfo.dataPath)));
    retinotopyInfo.sessionNum(iScan) = str2num(filenameSplit{end-2});
    retinotopyInfo.runNum(iScan) = str2num(filenameSplit{end-1});
    % now get the direction of motion
    thisDir = stimfile{iStimfile}.runInfo.retinotopy.dir(fileNum);    
    thisContrast = stimfile{iStimfile}.runInfo.retinotopy.contrast(fileNum);    
    % now either keep it or check it against previous value
    if isnan(retinotopyInfo.dir(iScan))
      retinotopyInfo.dir(iScan) = thisDir;
    elseif retinotopyInfo.dir(iScan) ~= thisDir
      disp(sprintf('(motexRet:checkRetinotopyScans) Direction of motion for scan %i (%i) does not match previous scans %i',iScan,thisDir,retinotopyInfo.dir(iScan)));
      retinotopyInfo.isret(iScan) = false;
      break;
    end      
    % same for contrast value
    if isnan(retinotopyInfo.contrast(iScan))
      retinotopyInfo.contrast(iScan) = thisContrast;
    elseif retinotopyInfo.contrast(iScan) ~= thisContrast
      disp(sprintf('(motexRet:checkRetinotopyScans) Contrast for scan %i (%i) does not match previous scans %i',iScan,thisContrast,retinotopyInfo.contrast(iScan)));
      retinotopyInfo.isret(iScan) = false;
      break;
    end      
  end
  % grab the runInfo - should be the same for all stimfies so
  % just choose the last one
  retinotopyInfo.runInfo{iScan} = stimfile{end}.runInfo;
  % update disppercent
  disppercent(iScan/nScans);
end
disppercent(inf);

% number of scans
retinotopyInfo.n = length(retinotopyInfo.isret);

% now check to see if there is a pre-computed retinotopy analysis
retinotopyInfo.corAnal = [];
corAnalFullName = setext(fullfile(viewGet(v,'datadir'),'corAnal',corAnalName),'mat');
if ~isfile(corAnalFullName)
  disp(sprintf('(motexRet:checkRetinotopyScans) No %s in %s',corAnalName,viewGet(v,'GroupName')));
  return
end

% ok, load the analysis
v = loadAnalysis(v,corAnalName);
retinotopyInfo.corAnal.name = corAnalName;

% now figure out which overlays are the important ones
retinotopyInfo.corAnal.scans = find((retinotopyInfo.isret) & (retinotopyInfo.contrast ~= 0));
retinotopyInfo.corAnal.dirs = retinotopyInfo.dir(retinotopyInfo.corAnal.scans);

% and grab the overlays for those
retinotopyInfo.corAnal.co = [];
retinotopyInfo.corAnal.ph = [];
retinotopyInfo.corAnal.amp = [];
for iScan = retinotopyInfo.corAnal.scans
  retinotopyInfo.corAnal.co(end+1,:,:) = viewGet(v,'overlayData',iScan,viewGet(v,'overlayNum','co'));
  retinotopyInfo.corAnal.ph(end+1,:,:) = viewGet(v,'overlayData',iScan,viewGet(v,'overlayNum','ph'));
  retinotopyInfo.corAnal.amp(end+1,:,:) = viewGet(v,'overlayData',iScan,viewGet(v,'overlayNum','amp'));
end

%%%%%%%%%%%%%%%%%%%
%    checkView    %
%%%%%%%%%%%%%%%%%%%
function [tf v] = checkView(v)

tf = false;

% check that the passed in structure is a view
if ~isview(v)
  disp(sprintf('(motexRet:checkView) Passed in structure is not a view'));
  return
end

% check for averages group
averagesGroupNum = viewGet(v,'groupNum','Averages');
if isempty(averagesGroupNum),return,end
v = viewSet(v,'curGroup',averagesGroupNum);

% so far so good
tf = 1;

