% motexSetRetinotopyStimImages.m
%
%      usage: motexSetRetinotopyStimImages()
%         by: justin gardner
%       date: 09/09/19
%    purpose: 
%
function retval = motexSetRetinotopyStimImages(varargin)

getArgs(varargin,{'retinotopyImageSequenceFilename=~/docs/2019/motex/retinotopy/retinotopyImageSequence.mat','stimImageWidth=90'});

% check for the retinotopy source data from Yuki
if ~isfile(retinotopyImageSequenceFilename)
  disp(sprintf('(motexSetRetinotopyStimImages) Could not find retinotopy image sequence: %s',retinotopyImageSequenceFilename));
  return
end

% load the retinotopy source data from Yuki
disppercent(-inf,'(motexSetRetinotopyStimImages) Loading original image sequence');
retinotopyImageSequence = load(retinotopyImageSequenceFilename);
disppercent(inf);

% easier names
X1 = retinotopyImageSequence.imageSequenceX1;
Y1 = retinotopyImageSequence.imageSequenceY1;

% scale down to a more useable size
[originalHeight originalWidth originalT] = size(X1);
[originalX originalY] = meshgrid(0:1/(originalWidth-1):1,0:1/(originalHeight-1):1);
stimImageHeight = round(originalHeight*stimImageWidth/originalWidth);
[newX newY] = meshgrid(0:1/(stimImageHeight-1):1,0:1/(stimImageWidth-1):1);

X1scale = nan(stimImageWidth,stimImageHeight,originalT);
disppercent(-inf,'(motexSetRetinotopyStimImages) Scaling images');
for iTime = 1:originalT
  X1scaled(:,:,iTime) = interp2(originalX,originalY,single(X1(:,:,iTime)),newX,newY,'nearest');
  disppercent(iTime/originalT);
end
disppercent(inf);

keyboard

