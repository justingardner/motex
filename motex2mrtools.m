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
function d = motex2mrtools(dataDir,varargin)

d = [];

% check arguments
if nargin < 1
  help motex2mrtools
  return
end

% process other arguments
getArgs(varargin,{'dataPath=~/data/motex'},'suppressUnknownArgMessage',true);

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

%%%%%%%%%%%%%%%%%%%%%%%%%
%    getMotexStimvol    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [d tf] = getMotexStimvol(d,varargin)

% default to failure
tf = false;

% parse arguments
getArgs(varargin,{'stimuli=stimuli-miniblock'},'suppressUnknownArgMessage',true);

% cycle through each session and run we are doing
for iSession = d.sessionNum
  for iRun = d.runNum{iSession}
    % shortcut to runInfo
    runInfo = d.runInfo{iSession}{iRun};
    
    switch (stimuli)
      case {'stimuli-miniblock'}
        [tf runInfo] = getMotexStimvolMiniblock(runInfo); 
	if ~tf, return, end
      otherwise 
        disp(sprintf('(getMotexStimvol) Unknown stimulus type: %s',stimuli));
	return
    end
    
    % set runInfo in structure
    d.runInfo{iSession}{iRun} = runInfo;
  end
end

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
  keyboard
  return
end

% check match between readme and phtoDiodes
nPhotoDiode = length(runInfo.photoDiodeTimes.allPhotoDiodeStartTime);
nReadmeStimuli = length(runInfo.readme.frameNum);
if nPhotoDiode ~= nReadmeStimuli
  disp(sprintf('(motex2mrtools:getMotexStimvolMiniblock) Photodiode triggers do not match readme: %i ~= %i',nPhotoDiode,nReadmeStimuli));
  return
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

% in stimulusInfo, make numeric version of texFamily, texGenType
texFamilyNames = unique(runInfo.stimulusInfo.texFamily);
texGenTypeNames = unique(runInfo.stimulusInfo.texGenType);
texFolderNames = unique(runInfo.stimulusInfo.texFolderName);
% now get numeric value of each stimulus
for iStimulus = 1:length(runInfo.stimulusInfo.texFamily)
  texFamilyNum(iStimulus) = find(strcmp(runInfo.stimulusInfo.texFamily{iStimulus},texFamilyNames));
  texGenTypeNum(iStimulus) = find(strcmp(runInfo.stimulusInfo.texGenType{iStimulus},texGenTypeNames));
  texFolderNum(iStimulus) = find(strcmp(runInfo.stimulusInfo.texFolderName{iStimulus},texFolderNames));
end

% now link each stimulus presentation with the corresponding camera frame
for iFrame = 1:nPhotoDiode
  [timeDiff stimToCameraFrame(iFrame)] = min(abs(runInfo.photoDiodeTimes.allPhotoDiodeStartTime(iFrame)-runInfo.acqMeanTime));
  if timeDiff > runInfo.cameraFrameLen
    disp(sprintf('(motex2mrtools:getMotexStimulusInfo) Frame %i (%i of miniblock) is  %0.5fs away from a camera acquisition period which is longer than one camera frame',iFrame,rem(iFrame,runInfo.nImagesPerMiniblock),timeDiff));
  end
end

% now we can compute various stimvols
% Fix, Fix, Fix, start here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getMotexStimulusInfo    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d tf] = getMotexStimulusInfo(d,varargin)

% arguments
getArgs(varargin,{'stimuli=stimuli-miniblock','stimDir=/Volumes/GoogleDrive/My Drive/docs/2019/motex/Expt_stimuli'},'suppressUnknownArgMessage',true);

% default failure
tf = false;

% get stimulusInfo for each run we are doing
for iSession = d.sessionNum
  for iRun = d.runNum{iSession}
    % load the stimulusInfo 
    stimfile = fullfile(stimDir,stimuli,'stimulusInfo.mat');
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
    readmeFile = fullfile(stimDir,stimuli,'readme.txt');
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
