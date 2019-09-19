% motexGetRetinotopyStimImages.m
%
%      usage: motexGetRetinotopyStimImages()
%         by: justin gardner
%       date: 09/09/19
%    purpose: 
%
function retval = motexGetRetinotopyStimImages(varargin)

getArgs(varargin,{'retPath=~/docs/2019/motex/retinotopy/','originalFile=retinotopyImageSequence.mat','scaledFile=retinotopyImageSequenceSmall.mat','stimImageWidth=90','forceRescale=0','screenSize=[87.9 48.5]','screenDistance=20','screenRefreshRate=60'});

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
  retval.X1.im = scaleImages(retinotopyImageSequence.imageSequenceX1,stimImageWidth);
  retval.X2.im = scaleImages(retinotopyImageSequence.imageSequenceX2,stimImageWidth);
  retval.Y1.im = scaleImages(retinotopyImageSequence.imageSequenceY1,stimImageWidth);
  retval.Y2.im = scaleImages(retinotopyImageSequence.imageSequenceY2,stimImageWidth);

  % save it
  save(fullfile(retPath,scaledFile),'retval');

else
  % just load the already scaled file
  load(fullfile(retPath,scaledFile));
end

% create dimensions
retval.X1 = addDimensions(retval.X1,screenSize,screenDistance,screenRefreshRate);
retval.X2 = addDimensions(retval.X2,screenSize,screenDistance,screenRefreshRate);
retval.Y1 = addDimensions(retval.Y1,screenSize,screenDistance,screenRefreshRate);
retval.Y2 = addDimensions(retval.Y2,screenSize,screenDistance,screenRefreshRate);

% make into 0 for gray and 1 for any stimulus
retval.X1.im = makeContrastImage(retval.X1.im);
retval.X2.im = makeContrastImage(retval.X2.im);
retval.Y1.im = makeContrastImage(retval.Y1.im);
retval.Y2.im = makeContrastImage(retval.Y2.im);

%%%%%%%%%%%%%%%%%%%%%%%
%    addDimensions    %
%%%%%%%%%%%%%%%%%%%%%%%
function stimImage = addDimensions(stimImage,screenSize,screenDistance,screenRefreshRate)

% compute size in degrees of visual angle
screenSizeDeg = r2d(atan(screenSize/screenDistance));

% get x dimension in degs
xPixels = size(stimImage.im,2);
xDegs = 0:screenSizeDeg(1)/(xPixels-1):screenSizeDeg(1);
xDegs = xDegs - mean(xDegs);

% get y dimension in degs
yPixels = size(stimImage.im,1);
yDegs = 0:screenSizeDeg(2)/(yPixels-1):screenSizeDeg(2);
yDegs = yDegs - mean(yDegs);

% get x and y
[stimImage.x stimImage.y] = meshgrid(xDegs,yDegs);

% get time
stimImage.t = 0:1/screenRefreshRate:((size(stimImage.im,3)-1)*(1/screenRefreshRate));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeContrastImage    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function contrastImages = makeContrastImage(originalImages)

contrastImages = originalImages;
contrastImages(originalImages==128) = 0;
contrastImages(originalImages~=128) = 1;


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
