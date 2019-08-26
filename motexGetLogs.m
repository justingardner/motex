% motexGetLogs
%
%      usage: [d tf] = motexGetLogs(d)
%         by: justin gardner
%       date: 07/19/19
%    purpose: gets log info for a session, called after using motexGetRawInfo
%             tf returns whether successful
% 
%       e.g.: d = motexGetRawInfo('M190718_RN')
%             d = motexGetLogs(d);
%
%             To only get logs for a specific session / run
%             d = motexGetLogs(d,'sessionNum=1','runNum=1');
function [d tf] = motexGetLogs(d,varargin)

% default failure
tf = false;

% check arguments
if nargin < 1
  help motexGetLogs
  d = [];
  return
end

% process other arguments
[argNames argVals args] = getArgs(varargin,{'sessionNum=inf','runNum=inf'});

% check session and run num existense
[tf d] = motexSessionRunCheck(d,sessionNum,runNum);
if ~tf,return,end

% for each session
for iSession = d.sessionNum
  % for each run
  for iRun = d.runNum{iSession}
    % sync to AI 
    [tf d] = motexSyncToAI(d,iSession,iRun,args{:});
    if ~tf,return,end

    % read stimulus information
    [tf d] = motexGetStimulusInfo(d,iSession,iRun,args{:});
    if ~tf,return,end
  end
end

% success
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexSessionRunCheck    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf d] = motexSessionRunCheck(d,sessionNum,runNum)

% default to failed test
tf = false;

% check if asked for session exists
if ~isinf(sessionNum)
  if any(sessionNum < 1) | any(sessionNum > d.nSessions)
    dispHeader(sprintf('(motexGetLogs) Session %s does not exist',mlrnum2str(sessionNum,'sigfigs=0')));
    return
  end
else
  sessionNum = 1:d.nSessions;
end
% add sessions
d.sessionNum = sessionNum;

% check to see if asked for run exists
if ~isinf(runNum)
  for iSession = sessionNum
    if ~all(ismember(runNum,1:d.nRuns(iSession)))
      dispHeader(sprintf('(motexGetLogs) Run %s does not exist in session %s',mlrnum2str(runNum,'sigfigs=0'),mlrnum2str(sessionNum,'sigfigs=0')));
      return
    end
  end
  % keep run nums
  d.runNum{1:length(d.sessionNum)} = runNum;
else
  % figure out runs (i.e. get all runs for each session
  for iSession = sessionNum
    d.runNum{iSession} = 1:d.nRuns(sessionNum(iSession));
  end
end

% passed test
tf = true;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexGetStimulusInfo %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf d] = motexGetStimulusInfo(d,iSession,iRun,varargin)

% default to fail
tf = false;

% function to look into logs to get stimulus info

% check arguments
getArgs(varargin,{'sessionNum',1:d.nSessions,'runNum=inf'},'suppressUnknownArgMessage',true);

% shortcut to runInfo
runInfo = d.runInfo{iSession}{iRun};

% load the protocol
p = load(runInfo.protocolPath);
if ~isfield(p,'Protocol')
  disp(sprintf('(motexGetLogs:motexGetStimulusInfo) No Protocol variable in %s',runInfo.protocolPath));
  return
else
  runInfo.protocol = p.Protocol;
end

% load the p file
fPFile = fopen(runInfo.pFilePath,'r');
if fPFile == 0
  disp(sprintf('(motexGetLogs:motexGetStimulusInfo) Could not open %s',runInfo.pFilePath));
else
  % read the header
  header = textscan(fPFile,'%s',2);
  runInfo.pFile.filename = header{1}{1};
  runInfo.pFile.type = header{1}{2};
  dims = textscan(fPFile,'%d %d %d',2);
  runInfo.pFile.dims = [dims{1}(1) dims{2}(1)];
  trials = textscan(fPFile,repmat('%d',1,runInfo.pFile.dims(2)),runInfo.pFile.dims(1));
  runInfo.pFile.trials = cell2mat(trials);
end
fclose(fPFile);

% update runInfo in structure
d.runInfo{iSession}{iRun} = runInfo;

% success
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%
%    motexSyncToAI    %
%%%%%%%%%%%%%%%%%%%%%%%
function [tf d] = motexSyncToAI(d,iSession,iRun,varargin)

% default to fail
tf = false;

% function will find sync times in analog data files for
% each one of the images

% check arguments
getArgs(varargin,{'aiTraceTime=1','aiTracePhotoDiode=2','aiTraceCameraOn=3','aiTraceImageAcq=4','ttlCutoff=2.5'},'suppressUnknownArgMessage',true);

% shortcut to runInfo
runInfo = d.runInfo{iSession}{iRun};
% check for log file
if ~isfile(runInfo.logPath)
  disp(sprintf('(motexGetLogs:motexSyncToAI) Could not find log %s',runInfo.logPath));
  return
end

% load the log
disppercent(-inf,sprintf('(motexGetLogs:motexSyncToAI) Loading log %s',runInfo.logPath));
log = load(runInfo.logPath);
disppercent(inf);
if ~isfield(log,'AIdata')
  disp(sprintf('(motexGetLogs:motexSyncToAI) Log file %s does not contain AIdata',runInfo.logPath));
  return
end

% now figure out every time that the photoDiode trace went high
runInfo.photoDiodeTimes = motexGetPhotoDiodeTimes(log.AIdata(:,aiTraceTime),log.AIdata(:,aiTracePhotoDiode),log.AIdata(:,aiTraceImageAcq));

% figure out times that the camera was on
cameraOnEdges = motexGetEdges(log.AIdata(:,aiTraceCameraOn),ttlCutoff);
runInfo.cameraOnTimes = log.AIdata(cameraOnEdges.rising,aiTraceTime);
runInfo.cameraOffTimes = log.AIdata(cameraOnEdges.falling,aiTraceTime);
runInfo.cameraOnN = length(runInfo.cameraOnTimes);

% figure out each image acquisition time
acqEdges = motexGetEdges(log.AIdata(:,aiTraceImageAcq),ttlCutoff);
runInfo.acqOnTimes = log.AIdata(acqEdges.rising,aiTraceTime);
runInfo.acqOffTimes = log.AIdata(acqEdges.falling,aiTraceTime);

% compute mean time of camera acquisition
runInfo.acqMeanTime = runInfo.acqOnTimes + (runInfo.acqOffTimes-runInfo.acqOnTimes)/2;

% figure out how many images per each cameraOnTime
runInfo.whichCameraOn = zeros(1,length(runInfo.acqOnTimes));
for iCameraOn = 1:length(runInfo.cameraOnTimes)
  % find out all the acqOnTimes that start after cameraOn and before cameraOff
  runInfo.whichCameraOn(find((runInfo.acqOnTimes >= runInfo.cameraOnTimes(iCameraOn)) & (runInfo.acqOnTimes <= runInfo.cameraOffTimes(iCameraOn)))) = iCameraOn;
  runInfo.acqPerCameraOn(iCameraOn) = sum(runInfo.whichCameraOn==iCameraOn);
end
    
% compute how long each frame is
runInfo.cameraFrameLen = median(runInfo.acqOffTimes - runInfo.acqOnTimes);
runInfo.cameraFreq = 1/runInfo.cameraFrameLen;
    
% display some information
disp(sprintf('(motexGetLogs:syncToAI) %s Session %i Run %i: %i Camera trials (%s images per trial) %0.2f Hz',d.dataDir,iSession,iRun,runInfo.cameraOnN,mlrnum2str(unique(runInfo.acqPerCameraOn),'sigfigs=0'),runInfo.cameraFreq));
  
% check to make sure photodiode trials match camera triggers
if runInfo.cameraOnN ~= runInfo.photoDiodeTimes.nTrials
  disp(sprintf('(motexGetLogs:syncToAI) Camera on triggers do not match photoDiode trials: %i ~= %i',runInfo.cameraOnN,runInfo.photoDiodeTimes.nTrials));
  if ~askuser('Is this Ok'),keyboard,end
end
    
% update runInfo in structure
d.runInfo{iSession}{iRun} = runInfo;

% success
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexGetPhotoDiodeTimes    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function photoDiodeTimes = motexGetPhotoDiodeTimes(t,photoDiode,cameraTrace);

% default return value
photoDiodeTimes = [];

% make into arrays
t = t(:);photoDiode = photoDiode(:);

% get min and max
minPhotoDiode = min(photoDiode(:));
maxPhotoDiode = max(photoDiode(:));

% check that there were some reasonable events
if maxPhotoDiode < 1
  disp(sprintf('(motexGetLogs:motexGetPhotoDiodeTimes) The photoDiode trace does not have any values above 1'));
  return
end

% calculate a cutoff to find each video frame
cutoff = minPhotoDiode + (maxPhotoDiode-minPhotoDiode)*0.2;
edges = motexGetEdges(photoDiode,cutoff);

% get length of edges
len = median(edges.len(2:end-1));

% remove any last edges for which there is no ending event
edgesNotFullyRecorded = find((edges.rising+len) > length(photoDiode));
if ~isempty(edgesNotFullyRecorded)
  edges.n = min(edgesNotFullyRecorded)-1;
  edges.rising = edges.rising(1:edges.n);
  edges.falling = edges.falling(1:edges.n);
  edges.completeN = edges.n;
  edges.len = edges.len(1:edges.n);
end

% get the height of each edge by putting each edge
% into a matrix and taking median over the edge
for iLen = 1:len
  edgeHeight(iLen,:) = photoDiode(edges.rising+iLen-1);
end
edgeHeight = median(edgeHeight)';

% sometimes there are a few stray max edgeHeights, so to
% be robust to these we chose the max and min as the 95% and 5%
% highest / lowest edges
edgeHeightSorted = sort(edgeHeight);
maxEdgeHeight = edgeHeightSorted(round(0.95*length(edgeHeight)));
minEdgeHeight = edgeHeightSorted(round(0.05*length(edgeHeight)));

% get a cutoff
cutoff = minEdgeHeight + (maxEdgeHeight-minEdgeHeight)*0.6;

% now find edges of these - this will be when the
% video trigger was on
videoTrigger = motexGetEdges(edgeHeight,cutoff);

% get mean and standard deviation of the no trigger
% use this as a stricter cutoff 
nullVideoTriggerHeight = edgeHeight(1:videoTrigger.rising(1)-15);
nullMean = mean(nullVideoTriggerHeight);
nullStd = std(nullVideoTriggerHeight);
strictCutoff = nullMean + 4*nullStd;

% now push the start time back to the very beginning of the event
% by moving forward until we are a few standard deviations above the mean
% of the no events
for iVideoTrigger = 1:videoTrigger.n
  % get the rising edge
  rising = videoTrigger.rising(iVideoTrigger);
  % go backwards while we are still above criteria
  stepsBack = 0;
  maxStepsBack = 10;  
  while (edgeHeight(rising) > strictCutoff) & (stepsBack<maxStepsBack)
    rising = rising -1;
    stepsBack = stepsBack + 1;
  end
  videoTrigger.expandedRising(iVideoTrigger) = rising;
end

% now convert these times back to the original time series
rising = edges.rising(videoTrigger.expandedRising);
falling = edges.falling(videoTrigger.falling);

% make a new trace for visualization
trigTrace = zeros(1,length(photoDiode));
for iTrig = 1:length(rising)
  trigTrace(rising(iTrig):falling(iTrig)) = 1;
end

% now analyze to see how many trials and frames
% first figure out the biggest jump between trigs
samplesBetweenTriggers = rising(2:end)-falling(1:end-1);
% find a big jump between triggers
jumpSizes = unique(samplesBetweenTriggers);
[~,jumpLoc] = max(diff(jumpSizes));
trialJumpSize = jumpSizes(jumpLoc+1);
% check for runs in which all of the steps are large, if this is 
% the case then choose the smallest of these to be a jump between runs
% use arbitrary cutoff of 3 seconds
if t(min(jumpSizes)) > 2
  trialJumpSize = min(jumpSizes);
end
% now find all the places where the trigger times exceed that (this
% should be a trial step
trialTrigs = [1 find(samplesBetweenTriggers>=trialJumpSize)+1];
nTrials = length(trialTrigs);
trigsPerTrial = unique(diff(trialTrigs));
if length(trigsPerTrial) > 1
  disp(sprintf('(motexGetLogs:motexGetPhotoDiodeTimes) %i trials do not have all the same length of photo triggers %s',nTrials,mlrnum2str(trigsPerTrial)));
end
% compute length in seconds of first trial
trialLen = t(falling(trialTrigs(2)))-t(rising(1));
titleStr = sprintf('Found %i trials of length %0.3fs (%s photo diode triggers per trial)',nTrials,trialLen,mlrnum2str(trigsPerTrial,'sigfigs=0'));
disp(sprintf('(motexGetPhotoDiodeTimes) %s',titleStr));
trialStartTime = t(edges.rising(videoTrigger.expandedRising(trialTrigs)));

% display
hFig = mlrSmartfig('motexGetLogs','reuse');clf;
subplot(3,1,1);
plot(t,photoDiode,'k-');
hold on
plot(t,trigTrace*(maxPhotoDiode-minPhotoDiode)+minPhotoDiode,'r-');
xlabel('time (sec)');
axis tight;
title(titleStr);

subplot(3,1,2);
startIndex = 1;
[~,endIndex] = min(abs(t-(trialStartTime(1)+5)));
plot(t(startIndex:endIndex),photoDiode(startIndex:endIndex),'k-');
hold on
plot(t(startIndex:endIndex),trigTrace(startIndex:endIndex)*(maxPhotoDiode-minPhotoDiode)+minPhotoDiode,'r-');
plot(t(startIndex:endIndex),cameraTrace(startIndex:endIndex)*(maxPhotoDiode-minPhotoDiode)+minPhotoDiode,'g-');
axis tight
legend('photo diode','frame detect','camera acquisition','Location','NorthWest');
xlabel('time (sec)');
ylabel('Analog amplitude (V)');
title(sprintf('Time up to first trial'));

subplot(3,1,3);
startIndex = 1;
[~,startIndex] = min(abs(t-(trialStartTime(1)-1)));
[~,endIndex] = min(abs(t-(trialStartTime(2)+1)));
plot(t(startIndex:endIndex),photoDiode(startIndex:endIndex),'k-');
hold on
plot(t(startIndex:endIndex),trigTrace(startIndex:endIndex)*(maxPhotoDiode-minPhotoDiode)+minPhotoDiode,'r-');
plot(t(startIndex:endIndex),cameraTrace(startIndex:endIndex)*(maxPhotoDiode-minPhotoDiode)+minPhotoDiode,'g-');
axis tight
xlabel('time (sec)');
title(sprintf('First trial: %i photo diode triggers',trialTrigs(2)-trialTrigs(1)));

% return computed info about the
photoDiodeTimes.nTrials = nTrials;
photoDiodeTimes.trialLen = trialLen;
photoDiodeTimes.trialStartTime = trialStartTime;
photoDiodeTimes.trialStartIndex = trialTrigs;
photoDiodeTimes.allPhotoDiodeStartTime = t(rising);
photoDiodeTimes.allPhotoDiodeEndTime = t(falling);

if ~askuser('Triggers ok')
  keyboard
end

close(hFig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% motexGetEdges.m
%
%      usage: motexGetEdges(timeseries,cutoff)
%         by: justin gardner
%       date: 12/08/03
%       e.g.: motexGetEdges(timeseries,4.5)
%    purpose: returns edge fall and rise times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function edges = motexGetEdges(timeseries,cutoff,dilation)

% check command line arguments
if (nargin == 2)
  dilation = 1;
elseif (nargin ~=3)
  help motexGetEdges;
  return
end

timeseries = timeseries(:)';

% find when timecourse exceeds cutoff value (if cutoff is negative
% assume that means less than cutoff. If cutoff is positive assume
% we are looking for values larger than cutoff).
if (cutoff < 0)
  cutofftimes = timeseries < cutoff;
else
  cutofftimes = timeseries > cutoff;
end

% no events, give up
if (isempty(find(cutofftimes)))
  edges.cutofftimes = cutofftimes;
  edges.rising = [];
  edges.falling = [];
  edges.n = 0;
  edges.dilated = edges.cutofftimes;
  return
end

% find rising edges
rising = [0 find(diff(find(cutofftimes)) > 1)]+1;
% make sure that last one is not problematic
if (rising(length(rising)) > length(cutofftimes))
  rising(length(rising)) = risingedgretimes(length(rising))-1;
end

% find falling edges
falling = [find(diff(find(cutofftimes)) > 1) length(find(cutofftimes))];

% match length
if (length(rising) ~= length(falling))
  falling = falling(1:length(rising));
end

% get times where the signal is over cutoff
findcutofftimes = find(cutofftimes);

% pack return structure
edges.cutofftimes = cutofftimes;
edges.rising = findcutofftimes(rising);
edges.falling = findcutofftimes(falling);

% dilate edges 
dilatedrise = edges.rising-dilation;
dilatedrise(dilatedrise <= 0) = 1;
dilatedfall = edges.falling+dilation;
dilatedfall(dilatedfall > length(timeseries)) = length(timeseries);

% set dilated edges to true
edges.dilated = cutofftimes;
for i = 1:length(dilatedrise);
  edges.dilated(dilatedrise(i):dilatedfall(i)) = 1;
end

edges.n = length(edges.rising);

% get edge length
edges.completeN = min(length(edges.rising),length(edges.falling));
edges.len = edges.falling(1:edges.completeN) - edges.rising(1:edges.completeN);
