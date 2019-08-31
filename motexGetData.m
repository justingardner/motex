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
%             To change paths you can set the following variables
%
%             fromDataPaths: A string (or cell array of strings, if there are multiple paths) which
%                            the camera data will be taken from
%              fromLogPaths: A string (or cell array of strings) for the path to the logs.
%                    toPath: A string specifying the path to where you want the data copied to.
%
%             e.g.
%             
%             motexGetData('M190802_RN','fromDataPaths=/Volumes/tex/IMAGING','fromLogPaths=/Volumes/DATA/MOUSE/LOGS','toPath=~/data/motex/raw');
%
%             If you want to download all availabel data (note that this uses rsync
%             so will only download data that has not been copied already):
%
%             motexGetData all
%
%             To download a single data set
%
%             motexGetData('M190802_RN');
%
%
function retval = motexGetData(dataName,varargin)

% check arguments
if nargin < 1
  help motexGetData;
  return
end

% get default args
getArgs(varargin,{'fromDataPaths',{'/Volumes/tex/IMAGING','/Volumes/tex/IMAGING/WIDEFIELD'},'fromLogPaths',{'/Volumes/DATA-1/MOUSE/LOGS','/Volumes/DATA/MOUSE/LOGS'},'logDirs',{'MPEP_LOGS','VS_LOGS'},'toPath=~/data/motex/raw','skipData=0'});

% see if we are being called to try to update all files
if strcmp(lower(dataName),'all')
  for iFromDataPath = 1:length(fromDataPaths)
    folders = dir(fromDataPaths{iFromDataPath});
    for iFolder = 1:length(folders)
      % check for folders that begin with M
      if folders(iFolder).name(1) == 'M'
	% call this function recursively
	motexGetData(folders(iFolder).name,varargin{:});
      end
    end
  end
  return;
end


% clear any funny characters
dataName = stripfilesep(strtrim(dataName));

if ~skipData
  % default to not found
  foundData = 0;

  % cycle through each possibe place data can be
  for iFromDataPath = 1:length(fromDataPaths)
    
    % look for data in each possible place
    fromDataName = fullfile(fromDataPaths{iFromDataPath},dataName);
    if isdir(fromDataName)
      % if found remember where it is found and break
      fromDataPath = fromDataPaths{iFromDataPath};
      foundData = 1;
      break;
    end
  end
  % could not find data
  if ~foundData
    % check for data
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




