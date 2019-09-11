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
getArgs(varargin,{'preProcessedDirs',{'~/data/motex/raw/retinotopy'},'resFac=1'});

% check that the view sturcture has a retinotopy scan
[tf v] = checkView(v);
if ~tf,return,end

% check retinotopy scans in averages
[v retinotopyInfo] = checkRetinotopyScans(v);
if ~any(retinotopyInfo.isret),return,end

% now try to load pre-processed retinotopy information
retinotopyInfo = loadPreProcessedRetinotopy(retinotopyInfo,preProcessedDirs);

% make pre-processed maps
retinotopyInfo = makePreProcessedMaps(retinotopyInfo,resFac);

% convert into rois for MLR
rois = makeMLRRois(v,retinotopyInfo);

% and save the rois
saveROI(v,rois);
keyboard

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
%    p.regions.coords{iConnComp} = [x(:)';y(:)'];
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
function [v retinotopyInfo] = checkRetinotopyScans(v)

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

