% motex2mrtools.m
%
%      usage: motex2mrtools(dataDir)
%         by: justin gardner
%       date: 07/19/19
%    purpose: Convert data formats from Benucci lab to mrTools
%       e.g.: motex2mrtools('M190621_MA')
%
%             To extract a particulare session / run
%             motex2mrtools('M190718_RN','sessionNum=1','runNum=1');
%
%             This function expects data to have been downloaded form the servers
%             using getMotexData to the directory ~/data/motex/raw and will
%             make the mrtools directory into ~/data/motex. This can be overwritten
%             by doing
%
%             motex2mrtools('M190718_RN','dataDir=from/path','toPath=to/path');
%
%             The program also needs to know about what stimulus was given which is
%             set by the stimulusType variable (which defaults to miniblock). You
%             can also do a very basic import in which you specify the times of
%             different events. stimTimes is cell array of scalars (or arrays) that
%             specify the time in seconds after the start of the trial when 
%             each one of the stimNames stimuli occur. SO, for example, if you had
%             the stimulus: 'texture' at 2.5s and 'phase-scramble' at 7.5s you would do
%             
%             motex2mrtools('M190621_MA','sessionNum=1','runNum=1','stimulusType=manual','stimTimes',{2.5 7.5},'stimNames',{'texture' 'phase-scramble'});
%
function d = motex2mrtools(dataDir,varargin)

d = [];

% check arguments
if nargin < 1
  help motex2mrtools
  return
end

% process other arguments
getArgs(varargin,{'dataPath=~/data/motex/raw'},'suppressUnknownArgMessage',true);

% get info about raw data
d = getMotexRawInfo(dataDir,'dataPath',dataPath);
if isempty(d), return, end

% get the log information
[d tf] = getMotexLogs(d,varargin{:});
if ~tf, return, end

% get the stimulusInfo
[d tf] = getMotexStimulusInfo(d,varargin{:});
if ~tf, return, end

% create the stimvol matrix
[d tf] = getMotexStimvol(d,varargin{:});
if ~tf, return, end

% write out to a new mrTools session
[d tf] = motexMakeSession(d,varargin{:});
if ~tf, return, end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexMakeSession    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d tf] = motexMakeSession(d,varargin)

% default to failure
tf = false;

getArgs(varargin,{'toPath=~/data/motex'},'suppressUnknownArgMessage',true);

% sessionName
nameSplit = strsplit(d.dataDir,'_');
sessionName = nameSplit{1};
% try to find operator initials
operator = 'unknown';
if (length(nameSplit) > 1) 
  initialLoc = regexp(nameSplit{2},'[a-zA-Z]+');
  if ~isempty(initialLoc)
    operator = nameSplit{2}(initialLoc(1):end);
  end
end

% display to user what we found
disp(sprintf('(motex2mrtools:motexMakeSession) Session name: %s Operator: %s',sessionName,operator));

% setup sesion name
sessionPath = fullfile(toPath,sessionName);
if isdir(sessionPath)
  % new name for session
  newSessionPath = fullfile(toPath,sprintf('delete_me_%s_%s_%s',sessionName,datestr(now,'yyyymmdd'),datestr(now,'hhmmss')));
  if askuser(sprintf('(motex2mrtools:motexMakeSession) Session already exists, move session to %s and start over',getLastDir(newSessionPath)));
    % move
    movefile(sessionPath,newSessionPath);
  end
end

if ~askuser(sprintf('(motex2mrtools:motexMakeSession) Create session directory: %s',sessionPath)),return,end
mkdir(sessionPath);

% make an anatomy folder
anatomyPath = fullfile(sessionPath,'Anatomy');
if ~isdir(anatomyPath)
  mkdir(anatomyPath);
end

% make an Etc folder
etcPath = fullfile(sessionPath,'Etc');
if ~isdir(etcPath)
  mkdir(etcPath);
end

% make directory for Raw
rawPath = fullfile(sessionPath,'Raw');
if ~isdir(rawPath)
  mkdir(rawPath);
end

% make tSeries path
tSeriesPath = fullfile(rawPath,'TSeries');
if ~isdir(tSeriesPath)
  mkdir(tSeriesPath);
end

% descriptions
description = {};
stimfileName = {};
sessionDescription = '';
runConcatAndEventRelated = {};
startScanNum = 1;endScanNum = 0;
% cycle through each session and run we are doing and save raw data
for iSession = d.sessionNum
  for iRun = d.runNum{iSession}
    % get runInfo
    runInfo = d.runInfo{iSession}{iRun};
    % add anatomy image  if there is not one already in there
    if ~isfile(fullfile(anatomyPath,'anatomy.nii'));
      mlrImageSave(fullfile(anatomyPath,'anatomy.nii'),runInfo.image,runInfo.hdr);
    end

    if ~isempty(sessionDescription)
      sessionDescription = sprintf('%s %s',sessionDescription,runInfo.stimulusInfo.stimulusType);
    else
      sessionDescription = runInfo.stimulusInfo.stimulusType;
    end
    
    % tell user what we are doing
    disppercent(-inf,sprintf('(motex2mrtools:motexMakeSession) Making session with %i raw camera files',runInfo.nFiles));

    % cycle through each image
    for iCamera = 1:runInfo.nFiles
      
      % keep track of scan number
      endScanNum = endScanNum+1;

      % filename
      filename = fullfile(runInfo.dataPath,sprintf('%s_%i.mat',d.dataDir,iCamera));
      % check for file
      if ~isfile(filename)
	disp(sprintf('(motex2mrtools) Could not find file %s',filename));
	return
      end
      
      % load the file
      [data hdr] = motexCamera2nifti(filename,'spoof3d');
      
      % make output filename
      outputFilename = fullfile(tSeriesPath,sprintf('%s_%i_%i_%04i.nii',d.dataDir,sessionNum,runNum,iCamera));
      
      % write file
      mlrImageSave(outputFilename,data,hdr);
      
      % make a description
      description{end+1} = runInfo.stimulusInfo.stimulusType;
      
      % add trial-wise description if there is one
      if isfield(runInfo,'description')
	description{end} = sprintf('%s%s',description{end},runInfo.description{iCamera});
      end     

      % if we have stimvols
      if isfield(runInfo,'stimvols')
	% make the stimFilename
	stimfileName{end+1} = fullfile(etcPath,sprintf('stimfile_%s_%i_%i_%04i.mat',d.dataDir,sessionNum,runNum,iCamera));
	
	% grab the stimvol for that run
	stimvol = runInfo.stimvols.stimvol{iCamera};

	% if there are labels then
	if isfield(runInfo.stimvols,'labels')
	  % save with labels
	  stimNames = runInfo.stimvols.labels;
	  save(stimfileName{end},'stimvol','stimNames','runInfo');
	else
	  % otherwise just save the stimvol
	  save(stimfileName{end},'stimvol','runInfo');
	end
      end

      % save the runInfo in Etc
      save(fullfile(etcPath,sprintf('runInfo_%i_%i.mat',iSession,iRun)),'runInfo');
      
      % update disppercent
      disppercent(iCamera/runInfo.nFiles);
    end
    disppercent(inf);
    
    % keep track of whether to run concat and event-related and which scans to do it over
    if isfield(runInfo,'stimvols')
      runConcatAndEventRelated{end+1} = [startScanNum endScanNum];
    end
    startScanNum = endScanNum+1;

  end
end

try
  % switch to directory
  curpwd = pwd;
  cd(sessionPath);
  % make sure we are quit out of mrTools
  mrQuit(0);
  % get default params
  [sessionParams groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1');
  % set description
  groupParams.description = description;
  % add session fields
  sessionParams.description = sessionDescription;
  sessionParams.operator = operator;
  sessionParams.magnet = '';
  sessionParams.coil = '';
  sessionParams.pulseSequence = '';
  % now init 
  mrInit(sessionParams,groupParams,'makeReadme=0');
  % set all the stimfiles
  v = newView;
  if length(stimfileName) > 0
    for iStimfile = 1:length(stimfileName)
      viewSet(v,'stimfileName',getLastDir(stimfileName{iStimfile}),iStimfile,'Raw');
    end
  end
  % save the full d structure
  save(fullfile(etcPath,'rawInfo.mat'),'d');
  % run concatenation and event-related as needed
  for iConcatAndEventRelated = 1:length(runConcatAndEventRelated)
    % do concatenation
    [v params] = concatTSeries(v,[],'justGetParams=1','defaultParams=1');
    params.filterType = 'Detrend only';
    params.description = 'Concatenation of [x...x]';
    params.scanList = runConcatAndEventRelated{iConcatAndEventRelated}(1):runConcatAndEventRelated{iConcatAndEventRelated}(2);
    params = mlrFixDescriptionInParams(params);
    v = concatTSeries(v,params);
  end
  % run event-related on the concatenated scans above
  if ~isempty(runConcatAndEventRelated)
    % do event-related
    v = viewSet(v,'curGroup','Concatenation');
    v = viewSet(v,'curScan',viewGet(v,'nScans'));
    [v params] = eventRelated(v,[],'justGetParams=1','defaultParams=1');
    params.scanParams{1}.hdrlen = 2;
    v = eventRelated(v,params);
  end
  deleteView(v);
  mrQuit(0);
catch
  % some error, switch back to original path
  cd(curpwd);
  return
end

% switch back to original path
cd(curpwd);

% success
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%
%    getMotexStimvol    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [d tf] = getMotexStimvol(d,varargin)

% default to failure
tf = false;

% parse arguments
getArgs(varargin,{'stimulusType=miniblock'},'suppressUnknownArgMessage',true);

% cycle through each session and run we are doing
for iSession = d.sessionNum
  for iRun = d.runNum{iSession}
    % shortcut to runInfo
    runInfo = d.runInfo{iSession}{iRun};

    % set information about what type of stimulus is here
    runInfo.stimulusInfo.stimulusType = stimulusType;
    
    switch (stimulusType)
      case {'manual'}
        [tf runInfo] = getMotexStimvolManual(runInfo,varargin{:}); 
	if ~tf, return, end
      case {'miniblock'}
        [tf runInfo] = getMotexStimvolMiniblock(runInfo); 
	if ~tf, return, end
      otherwise 
        disp(sprintf('(getMotexStimvol) Unknown stimulus type: %s',stimulusType));
	return
    end
    
    % set runInfo in structure
    d.runInfo{iSession}{iRun} = runInfo;
  end
end

% success
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getMotexStimvolManual    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf runInfo] = getMotexStimvolManual(runInfo,varargin)

% default failure
tf = false;

getArgs(varargin,{'stimTimes=[]','stimNames={}'},'suppressUnknownArgMessage',true);

% make sure that the camera on and photoDiode trial numbers match
if runInfo.cameraOnN ~= runInfo.photoDiodeTimes.nTrials
  disp(sprintf('(motex2mrtools:getMotexStimvolManual) Number of camera acq and photo didode start times do not match (%i ~= %i)',runInfo.cameraOnN,runInfo.photoDiodeTimes.nTrials));
  return
end

% figure out what camera frame these stimTimes corrspond to.
for iTrial = 1:runInfo.photoDiodeTimes.nTrials
  % find when the stimulus started this trial
  trialStartTime = runInfo.photoDiodeTimes.trialStartTime(iTrial);
  % find out the times of the camera acquisition
  cameraTimes = runInfo.acqMeanTime(runInfo.whichCameraOn == iTrial);
  for iStimType = 1:length(stimTimes)
    for iStim = 1:length(stimTimes{iStimType})
      % the time of this stimulus
      stimTime = trialStartTime+stimTimes{iStimType}(iStim);
      % find the closest match for this event
      [timeDiff cameraFrame] = min(abs(stimTime-cameraTimes));
      % make sure we are not more than half a camera frame away
      if timeDiff > (runInfo.cameraFrameLen/2)
	disp(sprintf('(motex2mrtools:getMotexStimvolManual) Stimulus time: %f on trial %i is %f from a camera acquistion which is more than half the frame time',stimTimes{iStimType}(iStim),iTrial,timeDiff));
	keyboard
      end
      % record this in stimvol
      stimvols.stimvol{iTrial}{iStimType}(iStim) = cameraFrame;
    end
  end
end

% the labels
if ~isempty(stimNames)
  stimvols.labels = stimNames;
end

% pack into structure to return
runInfo.stimvols = stimvols;

% success
tf = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getMotexStimvolMiniblock    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf runInfo] = getMotexStimvolMiniblock(runInfo)

% default failure
tf = false;

% first make sure that the stimulusInfo file has miniBlock in it
if ~isfield(runInfo.stimulusInfo,'miniBlockNumber')
  disp(sprintf('(motex2mrtools:getMotexStimvolMiniblock) Missing miniBlockNumber field in stimulusInfo'));
  return
end

% compute number of miniblocks
runInfo.nMiniblocks = length(unique(runInfo.stimulusInfo.miniBlockNumber));
% Check to see if the number of trials match how many mini-blocks there were
if runInfo.pFile.dims(1) ~= runInfo.nMiniblocks
  disp(sprintf('(motex2mrtools:getMotexStimvolMiniblock) Number of miniblocks %i does not match number of trials in pFile %i',runInfo.pFile.dims(1),runInfo.nMiniblocks));
  return
end

% check that the frame numbers match between photoDiode (which should have a trigger for each image)
% and the pFile which specifies the sequence order
if ~isequal(runInfo.photoDiodeTimes.trialStartIndex,runInfo.pFile.trials(:,2)')
  disp(sprintf('(motex2mrtools:getMotexStimvolMiniblock) Photodiode triggers do not match image number in pFile'));
  if ~askuser('(motex2mrtools:getMotexStimvolMiniblock) Do you want to continue')
    return
  end
end

% check match between readme and phtoDiodes
nPhotoDiode = length(runInfo.photoDiodeTimes.allPhotoDiodeStartTime);
nReadmeStimuli = length(runInfo.readme.frameNum);
if nPhotoDiode ~= nReadmeStimuli
  disp(sprintf('(motex2mrtools:getMotexStimvolMiniblock) Photodiode triggers do not match readme: %i ~= %i',nPhotoDiode,nReadmeStimuli));
  if ~askuser('(motex2mrtools:getMotexStimvolMiniblock) Do you want to continue')
    return
  end
end
% set these numbers into the strucutre
runInfo.nImages = nReadmeStimuli;
runInfo.nImagesPerMiniblock = runInfo.nImages/runInfo.nMiniblocks;

% convert the filenames in the stimulusInfo to ones without the fullfile
for iFilename = 1:length(runInfo.stimulusInfo.texFilename)
  [runInfo.stimulusInfo.texFilePath{iFilename} runInfo.stimulusInfo.texFilename{iFilename} ext] = fileparts(runInfo.stimulusInfo.texFilename{iFilename});
  runInfo.stimulusInfo.texFilename{iFilename} = setext(runInfo.stimulusInfo.texFilename{iFilename},ext);
end

% ok turn readme into stimulus numbers for convenience
stimNumsRaw = nan(1,nReadmeStimuli);
for iReadmeStimuli = 1:nReadmeStimuli
  if strcmp(runInfo.readme.type{iReadmeStimuli},'gray')
    stimNumsRaw(iReadmeStimuli) = 0;
  else
    stimNumsRaw(iReadmeStimuli) = find(strcmp(runInfo.readme.filename{iReadmeStimuli},runInfo.stimulusInfo.texFilename));
  end
end
% reshape into trials
stimNumsRaw = reshape(stimNumsRaw,runInfo.nImagesPerMiniblock,runInfo.nMiniblocks)';

% now reorder according to how the trials were run
stimNumsInRunOrder = stimNumsRaw(runInfo.protocol.seqnums,:);

% now create some simplified versions of this matrix, first
% one that only has the onsets of the stimuli and one that
% has the envelope in time (over the flickering of the stimulus)
stimNumsInRunOrderOnsetOnly = nan(size(stimNumsInRunOrder));
stimNumsInRunOrderEnvelope = nan(size(stimNumsInRunOrder));
curStimulus = [];curStimulusFirstFrame = nan;curStimulusLastFrame = nan;
for iFrame = 1:runInfo.nImagesPerMiniblock
  % set the frames to 0
  stimNumsInRunOrderOnsetOnly(:,iFrame) = 0;
  stimNumsInRunOrderEnvelope(:,iFrame) = 0;

  % march through frame by frame and see when the stimulus changes
  if any(stimNumsInRunOrder(:,iFrame)~=0)
    % see if we are at a new presentation of a stimulus
    if ~isequal(stimNumsInRunOrder(:,iFrame),curStimulus)
      % if this is not the first stimulus in the sequence then set the last stimulus
      if ~isempty(curStimulus)
	stimNumsInRunOrderOnsetOnly(:,curStimulusFirstFrame) = curStimulus;
	stimNumsInRunOrderEnvelope(:,curStimulusFirstFrame:curStimulusLastFrame) = repmat(curStimulus,1,curStimulusLastFrame-curStimulusFirstFrame+1);
      end
      % new current stimulus
      curStimulus = stimNumsInRunOrder(:,iFrame);
      curStimulusFirstFrame = iFrame;
      curStimulusLastFrame = iFrame;
    else
      % same stimulus, then update lastframe
      curStimulusLastFrame = iFrame;
    end
  end
  % set the last stimulus
  if ~isempty(curStimulus)
    stimNumsInRunOrderOnsetOnly(:,curStimulusFirstFrame) = curStimulus;
	stimNumsInRunOrderEnvelope(:,curStimulusFirstFrame:curStimulusLastFrame) = repmat(curStimulus,1,curStimulusLastFrame-curStimulusFirstFrame+1);
  end
end

% make trial-by-trial stimvols
for iTrial = 1:runInfo.nMiniblocks
  % get stim nums for this trial
  stimNumsThisTrial = stimNumsInRunOrderOnsetOnly(iTrial,:);
  % get time points when the camera was on this trial
  cameraTimes = runInfo.acqMeanTime(find(runInfo.whichCameraOn==iTrial));
  % get time points when each stimulus came on
  imageStartIndex = runInfo.photoDiodeTimes.trialStartIndex(iTrial);
  imageEndIndex = imageStartIndex+runInfo.nImagesPerMiniblock-1;
  % if the trial ended early then we need to fix things
  % so check what the start of the next trial is
  if iTrial < runInfo.photoDiodeTimes.nTrials
    nextImageStartIndex = runInfo.photoDiodeTimes.trialStartIndex(iTrial+1);
  else
    nextImageStartIndex = length(runInfo.photoDiodeTimes.allPhotoDiodeStartTime)+1;
  end
  % see if we have gone off into the next trial
  if imageEndIndex >= nextImageStartIndex
    % set how many photo diode triggers we have
    nPhotoDiodeThisTrial = nextImageStartIndex-imageStartIndex;
    disp(sprintf('(motex2mrtools) Trial %i has %i photoDiodeTimes when it should have %i. Assuming that the trial ended early',iTrial,nPhotoDiodeThisTrial,imageEndIndex-imageStartIndex+1));
    imageEndIndex = nextImageStartIndex-1;
    % remove and stimNums that happen after the end
    stimIndexes = find(stimNumsThisTrial);
    missingStimIndexes = find(stimIndexes>nPhotoDiodeThisTrial);
    if ~isempty(missingStimIndexes)
      disp(sprintf('(motex2mrtools) Stimulus that happened at frame %s needs to be dropped',mlrnum2str(stimIndexes(missingStimIndexes),'sigfigs=0')));
      % set these to zero
      stimNumsThisTrial(stimIndexes(missingStimIndexes)) = 0;
    end
  end
  % now get all the image times we can
  imageTimes = runInfo.photoDiodeTimes.allPhotoDiodeStartTime(imageStartIndex:imageEndIndex);
  % now get times of each video frame where an image was shown in trial
  imageTimes = imageTimes(find(stimNumsThisTrial~=0));
  % convert these imageTimes into cameraFrames
  for iImage = 1:length(imageTimes)
    [timediff,cameraFrame(iImage)] = min(abs(imageTimes(iImage)-cameraTimes));
    % the difference in time between stimulus and camera should be less
    % than half a frame
    if timediff>(runInfo.cameraFrameLen/2)
      disp(sprintf('(motex2mrtools) Image %i in trial %i happened %fs from a camera frame which is more than half a camera frame length',iImage,iTrial,timediff));
    end
  end
  % get just the imageNums for the non gray frames
  imageNums = stimNumsThisTrial(stimNumsThisTrial~=0);
  
  % now we are ready to make stimvols. We make them sorted for
  % the following variables
  makeStimvolsFor = {'texGenType','texFolderName','texFamily'};
  for iSortType = 1:length(makeStimvolsFor)
    thisStimvols = {};
    % make the labels for what each type is
    stimvols.(makeStimvolsFor{iSortType}).labels = unique(runInfo.stimulusInfo.(makeStimvolsFor{iSortType}));
    % find out label of each image in this sequence
    imageLabels = {runInfo.stimulusInfo.(makeStimvolsFor{iSortType}){imageNums}};
    % now go through each one of the labels
    for iLabel = 1:length(stimvols.(makeStimvolsFor{iSortType}).labels);
      % and find which one of the images match
      imageMatch = find(strcmp(imageLabels,stimvols.(makeStimvolsFor{iSortType}).labels{iLabel}));
      % and put the camera frames in here
      stimvols.(makeStimvolsFor{iSortType}).stimvol{iTrial}{iLabel} = cameraFrame(imageMatch);
    end
  end
end

% make a description
description = {};
for iTrial = 1:runInfo.nMiniblocks
  % find the names of the stimuli shown on this miniblock
  [~,~,stimNums] = find(stimNumsInRunOrderOnsetOnly(iTrial,:));
  stimNames = {runInfo.stimulusInfo.texFilename{stimNums}};
  % and create a string with all of the names
  description{end+1} = '';
  for iStim = 1:length(stimNames)
    description{end} = sprintf('%s:::%s',description{end},stripext(stimNames{iStim}));
  end
end

% set the stimvols that will be default used
stimvols.stimvol = stimvols.texFolderName.stimvol;
stimvols.labels = stimvols.texFolderName.labels;

% pack into structure to return
runInfo.stimvols = stimvols;
runInfo.stimNums.raw = stimNumsRaw;
runInfo.stimNums.inOrder = stimNumsInRunOrder;
runInfo.stimNums.inOrderOnsetOnly = stimNumsInRunOrderOnsetOnly;
runInfo.stimNums.inOrderEnvelope = stimNumsInRunOrderEnvelope;
runInfo.description = description;

% success
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getMotexStimulusInfo    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d tf] = getMotexStimulusInfo(d,varargin)

% arguments
getArgs(varargin,{'stimulusType=miniblock','stimDir=/Volumes/GoogleDrive/My Drive/docs/2019/motex/Expt_stimuli'},'suppressUnknownArgMessage',true);

% default failure
tf = false;

% if this is a manual stimulus type then don't do anything and return
if ~isempty(findstr('manual',stimulusType))
  disp(sprintf('(motex2mrtools:getMotexStimulusInfo) Manual stimulus type, not trying to load any stimulus info'));
  tf = true;
  return
end

% get stimulusInfo for each run we are doing
for iSession = d.sessionNum
  for iRun = d.runNum{iSession}
    % load the stimulusInfo 
    stimfile = fullfile(stimDir,stimulusType,'stimulusInfo.mat');
    if ~isfile(stimfile)
      disp(sprintf('(motex2mrtools:getMotexStimuli) Could not find file: %s',stimfile));
      return
    end

    % load stimfile
    s = load(stimfile);
    if ~isfield(s,'e')
      disp(sprintf('(motex2mrtools:getMotexStimuli) Missing expermient variable e in stimfile: %s',stimfile));
      return
    end

    % save the stimulusInfo
    d.runInfo{iSession}{iRun}.stimulusInfo = s.e;

    % check for readme
    readmeFile = fullfile(stimDir,stimulusType,'readme.txt');
    if ~isfile(readmeFile)
      disp(sprintf('(motex2mrtools:getMotexStimuli) Could not find readme file: %s',readmeFile));
      return
    end
      
    % load the readme
    fReadme = fopen(readmeFile);
    if fReadme == 0
      disp(sprintf('(motex2mrtools:getMotexStimuli) Could not open readme file: %s',readmeFile));
      return
    end
    
    % read the file
    s = textscan(fReadme,'%d %s %s');
    fclose(fReadme);
    
    % save the readme
    d.runInfo{iSession}{iRun}.readme.frameNum = s{1};
    d.runInfo{iSession}{iRun}.readme.type = s{2};
    d.runInfo{iSession}{iRun}.readme.filename = s{3};
    
  end
end

% success
tf = true;
