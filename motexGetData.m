% motexGetData.m
%
%        $Id:$ 
%      usage: motexGetData(dataName)
%         by: justin gardner
%       date: 06/26/19
%    purpose: rsync data over to google drive
%             will check both data and data-1 for correct server names for logs
%             Make sure you are connect via tunnel to 172.17.150.219/data (Labserver)
%             and also 172.17.150.7/tex for the data
%
%
function retval = motexGetData(dataName,varargin)

% check arguments
if nargin < 1
  help motexGetData;
  return
end

% get default args
getArgs(varargin,{'fromDataPath=/Volumes/tex/IMAGING','fromLogPaths',{'/Volumes/DATA-1/MOUSE/LOGS','/Volumes/DATA-1/MOUSE/LOGS'},'logDirs',{'MPEP_LOGS','VS_LOGS'},'toPath=~/data/motex','skipData=0'});

if ~skipData
  % look for data
  fromDataName = fullfile(fromDataPath,dataName);

  % check for data
  if ~isdir(fromDataName)
    disp(sprintf('(motexGetData) Could not find data: %s',fromDataName));
    return
  end
end

% check for log path
fromLogPath = '';
for iLogPath = 1:length(fromLogPaths)
  if isdir(fromLogPaths{iLogPath})
    fromLogPath = fromLogPaths{iLogPath};
  end
end

% look for logs
for iLog = 1:length(logDirs)
  % look for data
  fromLogName{iLog} = fullfile(fromLogPath,logDirs{iLog});
  % check for logs
  if ~isdir(fromLogName{iLog})
    disp(sprintf('(motexGetData) Could not find log directory: %s',fromLogName{iLog}));
    fromLogName{iLog} = '';
  end
end

if ~skipData
  % copy data
  command = sprintf('rsync -av --progress %s %s',fromDataName,fullfile(toPath));
  dispHeader(command);
  system(command);
end

% copy logs
for iLog = 1:length(logDirs)
  if ~isempty(fromLogName{iLog})
    command = sprintf('rsync -av --progress %s/%s* %s',fromLogName{iLog},dataName,fullfile(toPath,dataName,logDirs{iLog}));
    dispHeader(command);
    system(command);
  end
end




