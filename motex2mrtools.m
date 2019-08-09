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
getArgs(varargin,{'dataPath=~/data/motex'});

% get info about raw data
d = getMotexRawInfo(dataDir,'dataPath',dataPath);
if isempty(d),return,end

% get the log information
d = getMotexLogs(d,varargin{:});

keyboard