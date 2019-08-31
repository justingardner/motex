% motexGetStimvolFromStimnums.m
%
%      usage: stimvols = motexGetStimvolFromStimnums(runInfo,stimNums)
%         by: justin gardner
%       date: 08/20/19
%    purpose: Pass in runInfo and stimNums to get stimvols. used in motex2mrtools
%        
%
function stimvols = motexGetStimvolFromStimnums(runInfo,stimNums,varargin)

% check arguments
stimvols = [];
if nargin < 2
  help motexGetStimvolFromStimnums
  return
end
  
% get arguments
getArgs(varargin,{'makeStimvolsFor',{'texGenType','texFolderName','texFamily','texAll'}});

% make trial-by-trial stimvols
for iTrial = 1:runInfo.nMiniblocks
  % get stim nums for this trial
  stimNumsThisTrial = stimNums(iTrial,:);
  % get time points when the camera was on this trial
  cameraTimes = runInfo.acqMeanTime(find(runInfo.whichCameraOn==iTrial));
  % get time points when each stimulus came on
  imageStartIndex = runInfo.photoDiodeTimes.trialStartIndex(iTrial);
  imageEndIndex = imageStartIndex+runInfo.nImagesPerMiniblock-1;
  % if the trial ended early then we need to fix things
  % so check what the start of the next trial is
  if iTrial < runInfo.photoDiodeTimes.nTrials
    nextImageStartIndex = runInfo.photoDiodeTimes.trialStartIndex(iTrial+1);
  else
    nextImageStartIndex = length(runInfo.photoDiodeTimes.allPhotoDiodeStartTime)+1;
  end
  % see if we have gone off into the next trial
  if imageEndIndex >= nextImageStartIndex
    % set how many photo diode triggers we have
    nPhotoDiodeThisTrial = nextImageStartIndex-imageStartIndex;
    disp(sprintf('(motex2mrtools) Trial %i has %i photoDiodeTimes when it should have %i. Assuming that the trial ended early',iTrial,nPhotoDiodeThisTrial,imageEndIndex-imageStartIndex+1));
    imageEndIndex = nextImageStartIndex-1;
    % remove and stimNums that happen after the end
    stimIndexes = find(stimNumsThisTrial);
    missingStimIndexes = find(stimIndexes>nPhotoDiodeThisTrial);
    if ~isempty(missingStimIndexes)
      disp(sprintf('(motex2mrtools) Stimulus that happened at frame %s needs to be dropped',mlrnum2str(stimIndexes(missingStimIndexes),'sigfigs=0')));
      % set these to zero
      stimNumsThisTrial(stimIndexes(missingStimIndexes)) = 0;
    end
  end
  % now get all the image times we can
  imageTimes = runInfo.photoDiodeTimes.allPhotoDiodeStartTime(imageStartIndex:imageEndIndex);
  % now get times of each video frame where an image was shown in trial
  imageTimes = imageTimes(find(stimNumsThisTrial~=0));
  % convert these imageTimes into cameraFrames
  for iImage = 1:length(imageTimes)
    [timediff,cameraFrame(iImage)] = min(abs(imageTimes(iImage)-cameraTimes));
    % the difference in time between stimulus and camera should be less
    % than half a frame. Fudge a little since the camera frames might be slightly longer
    % than expected
    if timediff>(3*runInfo.cameraFrameLen/5)
      disp(sprintf('(motex2mrtools) Image %i in trial %i happened %fs from a camera frame which is more than half a camera frame length',iImage,iTrial,timediff));
    end
  end
  % get just the imageNums for the non gray frames
  imageNums = stimNumsThisTrial(stimNumsThisTrial~=0);
  
  % now we are ready to make stimvols. We make them sorted for
  % variables set in getArgs above - but making sure they exist
  makeStimvolsFor = intersect(makeStimvolsFor,fieldnames(runInfo.stimulusInfo));
  for iSortType = 1:length(makeStimvolsFor)
    thisStimvols = {};
    % make the labels for what each type is
    stimvols.(makeStimvolsFor{iSortType}).labels = unique(runInfo.stimulusInfo.(makeStimvolsFor{iSortType}));
    % find out label of each image in this sequence
    imageLabels = {runInfo.stimulusInfo.(makeStimvolsFor{iSortType}){imageNums}};
    % now go through each one of the labels
    for iLabel = 1:length(stimvols.(makeStimvolsFor{iSortType}).labels);
      % and find which one of the images match
      imageMatch = find(strcmp(imageLabels,stimvols.(makeStimvolsFor{iSortType}).labels{iLabel}));
      % and put the camera frames in here
      stimvols.(makeStimvolsFor{iSortType}).stimvol{iTrial}{iLabel} = cameraFrame(imageMatch);
    end
  end
end

