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
getArgs(varargin,{'dataPath=~/data/motex'});

% get info about raw data
d = getMotexRawInfo(dataDir,'dataPath',dataPath);
if isempty(d),return,end

% sync to AI 
d = motexSyncToAI(d);

%%%%%%%%%%%%%%%%%%%%%%%
%    motexSyncToAI    %
%%%%%%%%%%%%%%%%%%%%%%%
function d = motexSyncToAI(d,varargin)

% function will find sync times in analog data files for
% each one of the images
getArgs(varargin,{'sessionNum=inf','runNum=inf'});

% for each session
for iSession = 1:d.nSessions
  % figure out how many runs to run over
  if isinf(runNum) nRun = d.nRuns(iSession); else nRun = runNum; end
  % for each run
  for iRun = 1:nRun
    % shortcut to runInfo
    runInfo = d.runInfo{iSession}{iRun};
    % load the log
    if ~isfile(runInfo.logPath)
      disp(sprintf('(motex2mrtools:motexSyncToAI) Could not find log %s',runInfo.logPath));
      return
    end
  end
end


keyboard
