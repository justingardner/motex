% motex2mrtools.m
%
%      usage: motex2mrtools(dataDir)
%         by: justin gardner
%       date: 07/19/19
%    purpose: Convert data formats from Benucci lab to mrTools
%       e.g.: motex2mrtools('M190621_MA')
%
function d = motex2mrtools(dataDir,varargin)

d = [];

% check arguments
if nargin < 1
  help motex2mrtools
  return
end

% process other arguments
[argNames argVals args] = getArgs(varargin,{'dataPath=~/data/motex'});

% get info about raw data
d = getMotexRawInfo(dataDir,'dataPath',dataPath);
if isempty(d),return,end

% sync to AI 
d = motexSyncToAI(d,args{:});

%%%%%%%%%%%%%%%%%%%%%%%
%    motexSyncToAI    %
%%%%%%%%%%%%%%%%%%%%%%%
function d = motexSyncToAI(d,varargin)

% function will find sync times in analog data files for
% each one of the images

% check arguments
getArgs(varargin,{'aiTraceTime=1','aiTracePhotoDiode=2','sessionNum=inf','runNum=inf'});

% for each session
for iSession = 1:d.nSessions
  % figure out how many runs to run over
  if isinf(runNum) nRun = d.nRuns(iSession); else nRun = runNum; end
  % for each run
  for iRun = 1:nRun
    % shortcut to runInfo
    runInfo = d.runInfo{iSession}{iRun};
    % check for log file
    if ~isfile(runInfo.logPath)
      disp(sprintf('(motex2mrtools:motexSyncToAI) Could not find log %s',runInfo.logPath));
      continue;
    end
    % load the log
    disppercent(-inf,sprintf('(motex2mrtools:motexSyncToAI) Loading log %s',runInfo.logPath));
    log = load(runInfo.logPath);
    disppercent(inf);
    if ~isfield(log,'AIdata')
      disp(sprintf('(motex2mrtools:motexSyncToAI) Log file %s does not contain AIdata',runInfo.logPath));
      continue;
    end

    % now figure out every time that the photoDiode trace went high
    d.photoDiodeTimes = motexGetPhotoDiodeTimes(log.AIdata(:,aiTraceTime),log.AIdata(:,aiTracePhotoDiode));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexGetPhotoDiodeTimes    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function photoDiodeTimes = motexGetPhotoDiodeTimes(t,photoDiode);

% default return value
photoDiodeTimes = [];

% make into arrays
t = t(:);photoDiode = photoDiode(:);

% get min and max
minPhotoDiode = min(photoDiode(:));
maxPhotoDiode = max(photoDiode(:));

% check that there were some reasonable events
if maxPhotoDiode < 1
  disp(sprintf('(motex2mrtools:motexGetPhotoDiodeTimes) The photoDiode trace does not have any values above 1'));
  return
end

% calculate a cutoff to find each video frame
cutoff = minPhotoDiode + (maxPhotoDiode-minPhotoDiode)*0.2;
edges = motexGetEdges(photoDiode,cutoff);

% get length of edges
len = median(edges.len(2:end-1));

% get the height of each edge by putting each edge
% into a matrix and taking median over the edge
for iLen = 1:len
  edgeHeight(iLen,:) = photoDiode(edges.rising+iLen-1);
end
edgeHeight = median(edgeHeight)';

% find the min and max height of edges
minEdgeHeight = min(edgeHeight);
maxEdgeHeight = max(edgeHeight);

% get a cutoff
cutoff = minEdgeHeight + (maxEdgeHeight-minEdgeHeight)*0.4;

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
% now find all the places where the trigger times exceed that (this
% should be a trial step
trialTrigs = [1 find(samplesBetweenTriggers>=trialJumpSize)+1];
nTrials = length(trialTrigs);
trigsPerTrial = unique(diff(trialTrigs));
if length(trigsPerTrial) > 1
  disp(sprintf('(motex2mrtools:motexGetPhotoDiodeTimes) %i trials do not have all the same length of photo triggers %s',nTrials,mlrnum2str(trigsPerTrial)));
end
% compute length in seconds of first trial
trialLen = t(falling(trialTrigs(2)))-t(rising(1));
titleStr = sprintf('Found %i trials of length %0.3fs (%s photo diode triggers per trial)',nTrials,trialLen,mlrnum2str(trigsPerTrial,'sigfigs=0'));
disp(sprintf('(motexGetPhotoDiodeTimes) %s',titleStr));
trialStartTime = t(edges.rising(videoTrigger.expandedRising(trialTrigs)));

% display
mlrSmartfig('motex2mrtools','reuse');clf;
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
axis tight
xlabel('time (sec)');
title(sprintf('Time up to first trial'));

subplot(3,1,3);
startIndex = 1;
[~,startIndex] = min(abs(t-(trialStartTime(1)-1)));
[~,endIndex] = min(abs(t-(trialStartTime(2)+1)));
plot(t(startIndex:endIndex),photoDiode(startIndex:endIndex),'k-');
hold on
plot(t(startIndex:endIndex),trigTrace(startIndex:endIndex)*(maxPhotoDiode-minPhotoDiode)+minPhotoDiode,'r-');
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
