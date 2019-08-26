% motexGetRawInfo.m
%
%        $Id:$ 
%      usage: d = motexGetRawInfo(dataDir)
%         by: justin gardner
%       date: 07/21/19
%    purpose: Gets the info about the raw file locations for a motex session and returns a structure
%       e.g.: d = motexGetRawInfo('M190621_MA')
%
function d = motexGetRawInfo(dataDir,varargin)

d = [];

% check arguments
if nargin < 1
  help motexGetRawInfo
  return
end

% process other arguments
getArgs(varargin,{'dataPath=~/data/motex/raw'});

% check to see if the data directory is there
d.dataDir = dataDir;
d.dataFullPath = fullfile(dataPath,dataDir);
if ~isdir(d.dataFullPath)
  disp(sprintf('(motexGetRawInfo) Could not find data directory: %s',dataDir));
  return
end

% check for data sessions
d.nSessions = 1;
while(isdir(fullfile(d.dataFullPath,num2str(d.nSessions))))
  d.nSessions = d.nSessions+1;
end
d.nSessions = d.nSessions-1;
disp(sprintf('(motexGetRawInfo) Found %i sessions in %s',d.nSessions,d.dataDir));

% check for runs within sessions
for iSession = 1:d.nSessions
  nRuns = 1;
  while(isdir(fullfile(d.dataFullPath,num2str(iSession),num2str(nRuns))))
    nRuns = nRuns+1;
  end
  d.nRuns(iSession) = nRuns-1;
  % now count the number of files in each run
  for iRun = 1:d.nRuns(iSession)
    % get the dataPath
    d.runInfo{iSession}{iRun}.dataPath = fullfile(d.dataFullPath,num2str(iSession),num2str(iRun));
    % get the number of mat files that live there
    d.runInfo{iSession}{iRun}.nFiles = length(dir(fullfile(d.runInfo{iSession}{iRun}.dataPath,'*.mat')));

    % load a single file to get image information
    filename = fullfile(d.runInfo{iSession}{iRun}.dataPath,sprintf('%s_%i.mat',dataDir,1));
    % load file
    [data hdr] = motexCamera2nifti(filename);
    if isempty(data)
      % display info for this file
      d.runInfo{iSession}{iRun}.imageInfoStr = '[Image file load error]';
      % no data
      d.runInfo{iSession}{iRun}.image = [];
      d.runInfo{iSession}{iRun}.hdr = [];
    else
      % display info for this file
      d.runInfo{iSession}{iRun}.imageInfoStr = sprintf('dims: %s pixdims: %s freq: %s Hz',mlrnum2str(hdr.dim,'sigfigs=0'),mlrnum2str(hdr.pixdim,'sigfigs=5'),mlrnum2str(1/hdr.pixdim(3)));
      % keep the mean image and header
      d.runInfo{iSession}{iRun}.image = median(data,3);
      d.runInfo{iSession}{iRun}.hdr = hdr;
    end

    % get the path for the AIdata
    d.runInfo{iSession}{iRun}.logPath = fullfile(d.dataFullPath,'VS_LOGS',sprintf('%s_%i_%i.mat',dataDir,iSession,iRun));
    if ~isfile(d.runInfo{iSession}{iRun}.logPath)
      disp(sprintf('(motexGetRawInfo) Could not find VS_LOGS file: %s',d.runInfo{iSession}{iRun}.logPath));
      d.runInfo{iSession}{iRun}.logPath = '';
    end

    % try to load the screen info
    d.runInfo{iSession}{iRun}.screenInfoPath = fullfile(d.dataFullPath,'VS_LOGS',dataDir,sprintf('%s_%i_%i.mat',dataDir,iSession,iRun));
    d.runInfo{iSession}{iRun}.screenInfo = [];
    if ~isfile(d.runInfo{iSession}{iRun}.screenInfoPath)
      disp(sprintf('(motexGetRawInfo) Could not find screenInfo: %s',d.runInfo{iSession}{iRun}.screenInfoPath));
    else
      warning off
      load(d.runInfo{iSession}{iRun}.screenInfoPath);
      if exist('myScreenInfo','var')
	d.runInfo{iSession}{iRun}.screenInfo = myScreenInfo;
      end
      warning on
    end

    % get the path for the 
    d.runInfo{iSession}{iRun}.protocolPath = fullfile(d.dataFullPath,'MPEP_LOGS',dataDir,num2str(iSession),num2str(iRun),'Protocol.mat');
    if ~isfile(d.runInfo{iSession}{iRun}.protocolPath)
      disp(sprintf('(motexGetRawInfo) Could not find protocol: %s',d.runInfo{iSession}{iRun}.protocolPath));
      d.runInfo{iSession}{iRun}.protocolPath = '';
    end
    
    % check for p file
    d.runInfo{iSession}{iRun}.pFilePath = fullfile(d.dataFullPath,'MPEP_LOGS',dataDir,num2str(iSession),num2str(iRun),sprintf('%s_%i_%i.p',dataDir,iSession,iRun));
    if ~isfile(d.runInfo{iSession}{iRun}.pFilePath)
      disp(sprintf('(motexGetRawInfo) Could not find p file: %s',d.runInfo{iSession}{iRun}.pFilePath));
      d.runInfo{iSession}{iRun}.pFilePath = '';
    end
  end
end
disp(sprintf('(motexGetRawInfo) Found %s runs in sessions',mlrnum2str(d.nRuns,'sigfigs=0')));

for iSession = 1:d.nSessions
  for iRun = 1:d.nRuns(iSession)
    % display information
    dispHeader(sprintf('Session: %i Run: %i',iSession,iRun));
    disp(sprintf('%i Files: %s',d.runInfo{iSession}{iRun}.nFiles,d.runInfo{iSession}{iRun}.imageInfoStr));
    if ~isempty(d.runInfo{iSession}{iRun}.logPath)
      disp(sprintf('Log: %s',d.runInfo{iSession}{iRun}.logPath));
    else
      disp(sprintf('Log: MISSING'));
    end
    if ~isempty(d.runInfo{iSession}{iRun}.protocolPath)
      disp(sprintf('Protocol: %s',d.runInfo{iSession}{iRun}.protocolPath));
    else
      disp(sprintf('Protocol: MISSING'));
    end
    if ~isempty(d.runInfo{iSession}{iRun}.screenInfo)
      disp(sprintf('Screen: %s',d.runInfo{iSession}{iRun}.screenInfo.MonitorType));
    end
  end
end


