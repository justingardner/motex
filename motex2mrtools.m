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
    log = load(runInfo.logPath);
    if ~isfield(log,'AIdata')
      disp(sprintf('(motex2mrtools:motexSyncToAI) Log file %s does not contain AIdata',runInfo.logPath));
      continue;
    end
    keyboard
    % now figure out every time that the photoDiode trace went high
    motexGetPhotoDiodeTimes(log.AIdata(:,aiTraceTime),log.AIdata(:,aiTracePhotoDiode));
      

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

% these should be going on and off at the video frame rate, so check that
% by calculating how long they usually take
edgeLength = unique(diff(edges.rising(2:end)))
keyboard
% now go find the events
photoDiodeEvents = find(photoDiode>cutoff);
keyboard
