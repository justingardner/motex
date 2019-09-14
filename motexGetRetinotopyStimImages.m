% motexGetRetinotopyStimImages.m
%
%      usage: motexGetRetinotopyStimImages()
%         by: justin gardner
%       date: 09/09/19
%    purpose: 
%
function retval = motexGetRetinotopyStimImages(varargin)

getArgs(varargin,{'retPath=~/docs/2019/motex/retinotopy/','originalFile=retinotopyImageSequence.mat','scaledFile=retinotopyImageSequenceSmall.mat','stimImageWidth=90','forceRescale=0'});

% see if we need to load original and rescale.
if ~isfile(fullfile(retPath,scaledFile)) || forceRescale
  % get full filename
  originalFilename = fullfile(retPath,originalFile);
  if ~isfile(originalFilename)
    disp(sprintf('(motexGetRetinotopyStimImages) Could not find file: %s',originalFilename));
    return
  end
  % load the retinotopy source data from Yuki
  disppercent(-inf,sprintf('Loading original image sequence: %s',originalFilename));
  retinotopyImageSequence = load(originalFilename);
  disppercent(inf);

  % scale images
  retval.X1 = scaleImages(retinotopyImageSequence.imageSequenceX1,stimImageWidth);
  retval.X2 = scaleImages(retinotopyImageSequence.imageSequenceX2,stimImageWidth);
  retval.Y1 = scaleImages(retinotopyImageSequence.imageSequenceY1,stimImageWidth);
  retval.Y2 = scaleImages(retinotopyImageSequence.imageSequenceY2,stimImageWidth);

  % save it
  save(fullfile(fileparts(retPath,scaledFile)),'retval');
  
  return
end

% just load the already scaled file
load(fullfile(retPath,scaledFile));


%%%%%%%%%%%%%%%%%%%%%
%    scaleImages    %
%%%%%%%%%%%%%%%%%%%%%
function scaledImages = scaleImages(originalImages,stimImageWidth)

% scale down to a more useable size
[originalHeight originalWidth originalT] = size(originalImages);
[originalX originalY] = meshgrid(0:1/(originalWidth-1):1,0:1/(originalHeight-1):1);
stimImageHeight = round(originalHeight*stimImageWidth/originalWidth);
[newX newY] = meshgrid(0:1/(stimImageHeight-1):1,0:1/(stimImageWidth-1):1);

scaledImages = nan(stimImageWidth,stimImageHeight,originalT);
disppercent(-inf,'(motexGetRetinotopyStimImages) Scaling images');
for iTime = 1:originalT
  scaledImages(:,:,iTime) = interp2(originalX,originalY,single(originalImages(:,:,iTime)),newX,newY,'nearest');
  disppercent(iTime/originalT);
end
disppercent(inf);
