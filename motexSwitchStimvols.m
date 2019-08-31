% motexSwitchStimvols.m
%
%      usage: motexSwitchStimvols(v)
%         by: justin gardner
%       date: 08/20/19
%    purpose: Switch the stimvols to compute a different set of event-related functions
%
function retval = motexSwitchStimvols(v)

% check arguments
if ~any(nargin == [1])
  help motexSwitchStimvols
  return
end

% check for legitimate view
if ~isview(v)
  disp(sprintf('(motexSwitchStimvols) Passed in variable v is not a view'));
  return
end

% set to raw
v = viewSet(v,'curGroup','Raw');

% load runInfo
etcDir = viewGet(v,'etcDir');

% get the rawInfo file
rawInfoFilename = fullfile(etcDir,'rawInfo.mat');
rawInfo = [];
if isfile(rawInfoFilename)
  rawInfo = load(rawInfoFilename);
end
if ~isfield(rawInfo,'d')
  disp(sprintf('(motexSwitchSimvols) Problem loading rawInfo file: %s',rawInfoFilename));
  return
end
rawInfo = rawInfo.d;

currentRunInfo = [];

iRun = 0;
% now go through each raw scan
nScans = viewGet(v,'nScans');
for iScan = 1:nScans
  % get the tseries filename to find out which session and run
  tseriesFilename = viewGet(v,'tseriesFile',iScan);
  % split the filename
  filenameSplit = strsplit(tseriesFilename,'_');
  if length(filenameSplit) < 3
    disp(sprintf('(motexSwitchStimvols) Filename %s does not have session and run',tseriesFilename));
    return
  end
  % get the session and run
  sessionNum = str2num(filenameSplit{end-2});
  runNum = str2num(filenameSplit{end-1});
  % make sure the session exists
  if (sessionNum < 1) || (sessionNum > rawInfo.nSessions)
    disp(sprintf('(motexSwitchStimvols) Session number %i is out of range',sessionNum));
    keyboard
  end
  % make sure the run exists
  if ~ismember(runNum,rawInfo.runNum{sessionNum})
    disp(sprintf('(motexSwitchStimvols) Session number %i is out of range',sessionNum));
    keyboard
  end
  % now load the runInfo if necessary
  if ~isequal(currentRunInfo,[sessionNum runNum])
    % if we have already collected a run, then save what we have for that
    if iRun > 0
      run(iRun).runNumScans = runNumScans;
    end
    iRun = iRun+1;
    % get runinfo filename
    runInfoFilename = fullfile(etcDir,sprintf('runInfo_%i_%i.mat',sessionNum,runNum));
    runInfo = [];
    if isfile(runInfoFilename)
      runInfo = load(runInfoFilename);
    end
    if ~isfield(runInfo,'runInfo')
      disp(sprintf('(motexSwitchStimvols) Problem loading runInfo file: %s',runInfoFilename));
      keyboard
    end
    runInfo = runInfo.runInfo;
    % remember which one is loaded
    currentRunInfo = [sessionNum runNum];
    % get a new field which sorts by texFamily_texGenType_texFolderName
    for iTrial = 1:length(runInfo.stimulusInfo.texFamily)
      runInfo.stimulusInfo.texAll{iTrial} = sprintf('%s_%s_%s',runInfo.stimulusInfo.texFamily{iTrial},runInfo.stimulusInfo.texGenType{iTrial},runInfo.stimulusInfo.texFolderName{iTrial});
    end

    % now ask user to choose reordering
    paramsInfo{1} = {'type',putOnTopOfList('inOrder',fieldnames(runInfo.stimNums))};
    paramsInfo{2} = {'sort',putOnTopOfList('texAll',union(setdiff(fieldnames(runInfo.stimvols),{'stimvol','labels','hdrlen'}),'texAll'))};
    params = mrParamsDialog(paramsInfo,sprintf('Stimvols for session: %i run: %i',sessionNum,runNum));
    if isempty(params),return,end
    % get the stimvols
    stimvols = motexGetStimvolFromStimnums(runInfo,runInfo.stimNums.(params.type));
    run(iRun).stimvols = stimvols.(params.sort);
    run(iRun).n = length(run(iRun).stimvols.stimvol);
    run(iRun).runInfo = runInfo;
    % keep track of which scan numbers we have
    runNumScans = 1;
  else
    runNumScans = runNumScans+1;
  end
  % save which run each scan is
  whichRun(iScan) = iRun;
  whichScan(iScan) = runNumScans;
end
run(iRun).runNumScans = runNumScans;

% make sure the scans match
for iRun = 1:length(run)
  if run(iRun).runNumScans ~= run(iRun).n
    disp(sprintf('(motexSwitchStimvols) Number of scans: %i does not match expected: %i',run(iRun).runNumScans,run(iRun).n));
    keyboard
  end
end

% make backupdir name
backupDir = sprintf('stimfile_%s_%s',datestr(now,'yyyymmdd'),datestr(now,'hhmmss'));
backupDir = fullfile(etcDir,backupDir);

% ask user if we should continue
if ~askuser(sprintf('(motexSwitchStimvols) Ok to backup stimfiles to %s and link new ones',getLastDir(backupDir)))
  return
end

% full path of backupDir
mkdir(backupDir);

disppercent(-inf,'(motexSwitchStimvols) Replacing stimvols');
for iScan = 1:nScans
  % get current stimfilename
  stimfileName = viewGet(v,'stimfileName',iScan);
  % move it to the backup directory
  if isfile(stimfileName{1})
    movefile(stimfileName{1},backupDir);
  end
  % make new variable
  stimvol = run(whichRun(iScan)).stimvols.stimvol{whichScan(iScan)};
  stimNames = run(whichRun(iScan)).stimvols.labels;
  runInfo = run(whichRun(iScan)).runInfo;
  % save the file
  save(stimfileName{1},'stimvol','stimNames','runInfo');
  % disppercent
  disppercent(iScan/nScans);
end
disppercent(inf);


