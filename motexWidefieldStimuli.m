% motexWidefieldStimuli.m
%
%        $Id:$ 
%      usage: motexWidefieldStimuli()
%         by: justin gardner
%       date: 06/20/19
%    purpose: Pulls images from Akshay's folders and converts them into
%             tif files alll with the same dimensions in an order that
%             can be read by Yuki's program. Essentially the program
%             displays images one at a time. There is a setting for frequency
%             which was set to 20 (2 Hz or 500ms per). The control system
%             sets which images to load from which directory (note there
%             is an off by one thing there since the folder we put them
%             in was run19 but started with run0, so you have to list the
%             number 20. On each line, you put the start and end numbers
%             of images. This program generates gray images for the
%             interstimulus intervals. Each try was set to have the stimulus
%             and then the scramble. Some optins:
%
%             flicker = Generate a grayscale image after each stim image
%                       control program should specify twice the frequency: 4
%             outputSize = 255. Sets size to scale image to in pixels
%             randomizeOrder = 1. Randomize order of textures so that no
%                       texture family repeats
%
function retval = motexWidefieldStimuli(varargin)

getArgs(varargin,{'sourceDir=/Volumes/GoogleDrive/My Drive/docs/2019/NSF CRCNS/Stimuli','texFolder=tex_eq','noiseFolder=noise_eq','destDir=~/Desktop/expt','texFamily',{'glass','fronds','spikes','beans','crystals','rocks','scales','clouds'},'texGenTypes',{'pool4','PS'},'destType=tif','nImageRepeat=5','destFilenameStem=im','grayValue=128','readmeFilename=readme.txt','stimulusInfoFilename=stimulusInfo','trialTime=10','outputSize=255','flicker=1','iti=5','imageFrequency=2.0','randomizeOrder=0','generateMethod=2','preStimDur=1','postStimDur=1.5','stimDur=2.5','miniblock=1'});

% get files from both tex and noise
searchFolders = {texFolder, noiseFolder};

% get all texture classes
if isequal(lower(texFamily),'all')
  % display what we are doing
  dispHeader('(motexWidefieldStimuli) Finding all texture types');
  texFamily = getAllFieldNum(sourceDir,searchFolders,2);
end

% get all algorithms
if isequal(lower(texGenTypes),'all')
  % display what we are doing
  dispHeader('(motexWidefieldStimuli) Finding all texture generation algorithms');
  texGenTypes = getAllFieldNum(sourceDir,searchFolders,3);
end

% get all the matching textures
matchTex = getMatchTextures(sourceDir,searchFolders,texGenTypes,texFamily);

if generateMethod==1
  % generate the sequence the old way
  generateImageSequence1(sourceDir,searchFolders,destDir,matchTex,texFamily,texGenTypes,destType,nImageRepeat,destFilenameStem,grayValue,readmeFilename,stimulusInfoFilename,trialTime,outputSize,flicker,iti,imageFrequency,randomizeOrder);
else
  % generate the new way (independent trials for each stimulus
  generateImageSequence2(sourceDir,searchFolders,destDir,matchTex,texFamily,texGenTypes,destType,destFilenameStem,grayValue,readmeFilename,stimulusInfoFilename,outputSize,flicker,randomizeOrder,imageFrequency,preStimDur,stimDur,postStimDur,miniblock);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    gereateImageSequence2    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function generateImageSequence2(sourceDir,searchFolders,destDir,matchTex,texFamily,texGenTypes,destType,destFilenameStem,grayValue,readmeFilename,stimulusInfoFilename,outputSize,flicker,randomizeOrder,imageFrequency,preStimDur,stimDur,postStimDur,miniblock);

% confirm
if ~askuser('(motexWidefiledStimuli) Ok to copy'), return,end

% make the destination dir
makeDestDir(searchFolders,destDir);

% readme string
readmeText = '';

% set random number see to always generate the same sequence
e.randomSeed = 10;
rng(e.randomSeed);

% get the stimulus order
stimulusOrder = getRandomizeOrder(matchTex,randomizeOrder);

% set the folder number 
for iFolder = 1:length(searchFolders)
  folderNums(iFolder,1:length(stimulusOrder)) = iFolder;
end

% now make it the length of how many folders (i.e. one for tex and phase-scramble)
stimulusOrder = repmat(stimulusOrder,length(searchFolders),1);

% now make them into single array
stimulusOrder = stimulusOrder(:);
folderNums = folderNums(:);

% randomize for miniblocks
if miniblock>1
  % get randomized order
  e.randOrder = randperm(length(stimulusOrder));
  stimulusOrder = stimulusOrder(e.randOrder);
  folderNums = folderNums(e.randOrder);
end

% write out image sequence
imageNum = 1;iTrial = 0;
for iFile= 1:length(stimulusOrder)

  % get filenumber and folderNum
  fileNum = stimulusOrder(iFile);
  folderNum = folderNums(iFile);

  % keep information about image
  iTrial = iTrial + 1;
  e.texFamily{iTrial} = matchTex.family{1}{fileNum};
  e.texGenType{iTrial} = matchTex.genType{1}{fileNum};
  e.texSampleNum(iTrial) = matchTex.sampleNum{1}(fileNum);
  e.texFolderNum(iTrial) = folderNum;
  e.texFolderName{iTrial}= searchFolders{folderNum};
  e.miniBlockNumber(iTrial) = floor((iTrial-1)/miniblock)+1;
  e.miniBlockOrder(iTrial) = mod(iTrial-1,miniblock)+1;
  
  % get source filename
  sourceFilename = fullfile(matchTex.paths{folderNum}{fileNum},matchTex.filenames{folderNum}{fileNum});

  % get source filetype
  sourceType = getext(sourceFilename);

  % read image
  sourceImage = imread(sourceFilename,sourceType);

  % rescale
  imageSize = size(sourceImage);
  if (length(imageSize)>2) && (imageSize(3) > 1)
    sourceImage = mean(sourceImage,3);
  end
  [sourceMeshX sourceMeshY] = meshgrid((0:1/(imageSize(1)-1):1)-0.5,(0:1/(imageSize(2)-1):1)-0.5);
  [outputMeshX outputMeshY] = meshgrid((0:1/(outputSize-1):1)-0.5,(0:1/(outputSize-1):1)-0.5);
  sourceImage = interp2(sourceMeshX,sourceMeshY,double(sourceImage),outputMeshX,outputMeshY);
  sourceImage = uint8(sourceImage);

  % keep image
  e.texImage{iTrial} = sourceImage;
  e.texFilename{iTrial} = sourceFilename;

  % display image
  mlrSmartfig('motexWidefieldStimuli','reuse');clf;
  imagesc(sourceImage);colormap(gray);title(sourceFilename);drawnow
  disp(sprintf('(motexWidefileStimuli) ImageNum: %i Converting from %s: %s',imageNum,getLastDir(matchTex.paths{folderNum}{fileNum}),getLastDir(sourceFilename)));

  % make gray image
  grayImage = sourceImage;
  grayImage(:) = grayValue;

  % get number of images for each section
  nPreStimImages = round(preStimDur*imageFrequency);
  if flicker
    nStimImages = round(stimDur*imageFrequency/2);
  else
    nStimImages = round(stimDur*imageFrequency);
  end
  nPostStimImages = round(postStimDur*imageFrequency);

  if e.miniBlockOrder(end) == 1
    % write out the gray image for the prestimdur
    for iImageRepeat = 1:nPreStimImages
      % create destFilenamea
      destFilename = sprintf('%s%04i.%s',destFilenameStem,imageNum,destType);
      destFilename = fullfile(destDir,destFilename);
      % write out what was done
      readmeText = sprintf('%s%04i\tgray\tgray\n',readmeText,imageNum);
      % update image num
      imageNum = imageNum+1;
      % write image
      imwrite(grayImage,destFilename,destType,'Compression','lzw');
    end
  end

  % write out the image
  for iImageRepeat = 1:nStimImages
    % create destFilename
    destFilename = sprintf('%s%04i.%s',destFilenameStem,imageNum,destType);
    destFilename = fullfile(destDir,destFilename);
    % write out what was done
    readmeText = sprintf('%s%04i\t%s\t%s\n',readmeText,imageNum,getLastDir(matchTex.paths{folderNum}{fileNum}),getLastDir(matchTex.filenames{folderNum}{fileNum}));
    % update image num
    imageNum = imageNum+1;
    % write image
    imwrite(sourceImage,destFilename,destType,'Compression','lzw');
    % if flicker, then put in gray image, after image stimulus image
    if flicker
      % create destFilename
      destFilename = sprintf('%s%04i.%s',destFilenameStem,imageNum,destType);
      destFilename = fullfile(destDir,destFilename);
      % write out what was done
     readmeText = sprintf('%s%04i\tgray\tgray\n',readmeText,imageNum);
      % update image num
      imageNum = imageNum+1;
      % write image
      imwrite(grayImage,destFilename,destType,'Compression','lzw');
    end
  end

  % write out the gray image for the poststimdur
  for iImageRepeat = 1:nPostStimImages
    % create destFilename
    destFilename = sprintf('%s%04i.%s',destFilenameStem,imageNum,destType);
    destFilename = fullfile(destDir,destFilename);
    % write out what was done
    readmeText = sprintf('%s%04i\tgray\tgray\n',readmeText,imageNum);
    % update image num
    imageNum = imageNum+1;
    % write image
    imwrite(grayImage,destFilename,destType,'Compression','lzw');
  end
end

% write out readme file
f = fopen(fullfile(destDir,readmeFilename),'w');
fprintf(f,readmeText);
fclose(f);

% save stimulus information
save(fullfile(destDir,stimulusInfoFilename),'e');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    gereateImageSequence1    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function generateImageSequence1(sourceDir,searchFolders,destDir,matchTex,texFamily,texGenTypes,destType,nImageRepeat,destFilenameStem,grayValue,readmeFilename,stimulusInfoFilename,trialTime,outputSize,flicker,iti,imageFrequency,randomizeOrder)

% this is the first way we did it while I was in Japan

% calculate how long sequence will take
nFiles = length(matchTex.filenames{1});
nTexFamily = length(texFamily);
nTexGenTypes = length(texGenTypes);
nSamples = nFiles/(nTexFamily*nTexGenTypes);
timePerTextureClass = (trialTime+iti)*nTexFamily*nTexGenTypes;
timePerRun = timePerTextureClass*nSamples;
if flicker,imageFrequency = imageFrequency*2; end
nFrames = nTexFamily*nTexGenTypes*nSamples * (trialTime * imageFrequency);
disp(sprintf('(motexWidefieldStimuli) %i textures classes,  %i algorithms, %0.1f samples total: %i',nTexFamily,nTexGenTypes,nSamples,nFiles));
disp(sprintf('(motexWidefieldStimuli) %0.1fs (%s) per texture class (iti = %0.1fs)',timePerTextureClass,mlrDispElapsedTime(timePerTextureClass),iti));
disp(sprintf('(motexWidefieldStimuli) %0.1fs (%s) per run (iti = %0.1fs)',timePerRun,mlrDispElapsedTime(timePerRun),iti));
disp(sprintf('(motexWidefieldStimuli) %i frames',nFrames));


% confirm
if ~askuser('(motexWidefiledStimuli) Ok to copy'), return,end

% make the destination dir
makeDestDir(searchFolders,destDir);

% readme string
readmeText = '';

% get the stimulus order
stimulusOrder = getRandomizeOrder(matchTex,randomizeOrder);

% write out image sequence
imageNum = 1;iTrial = 0;
for iFile= 1:length(matchTex.filenames{1})
  % get filenumber
  fileNum = stimulusOrder(iFile);
  % keep information about image
  iTrial = iTrial + 1;
  e.texFamily{iTrial} = matchTex.family{1}{fileNum};
  e.texGenType{iTrial} = matchTex.genType{1}{fileNum};
  e.texSampleNum(iTrial) = matchTex.sampleNum{1}(fileNum);
  % create images
  for iFolder = 1:length(searchFolders)
    % get source filename
    sourceFilename = fullfile(matchTex.paths{iFolder}{fileNum},matchTex.filenames{iFolder}{fileNum});
    % get source filetype
    sourceType = getext(sourceFilename);
    % read image
    sourceImage = imread(sourceFilename,sourceType);
    % rescale
    imageSize = size(sourceImage);
    if (length(imageSize)>2) && (imageSize(3) > 1)
      sourceImage = mean(sourceImage,3);
    end
    [sourceMeshX sourceMeshY] = meshgrid((0:1/(imageSize(1)-1):1)-0.5,(0:1/(imageSize(2)-1):1)-0.5);
    [outputMeshX outputMeshY] = meshgrid((0:1/(outputSize-1):1)-0.5,(0:1/(outputSize-1):1)-0.5);
    sourceImage = interp2(sourceMeshX,sourceMeshY,double(sourceImage),outputMeshX,outputMeshY);
    sourceImage = uint8(sourceImage);
    % keep image
    if ~isempty(findstr(searchFolders{iFolder},'tex'))
      e.texImage{iTrial} = sourceImage;
      e.texFilename{iTrial} = sourceFilename;
    else
      e.controlImage{iTrial} = sourceImage;
      e.controlFilename{iTrial} = sourceFilename;
    end
    % display image
    imagesc(sourceImage);colormap(gray);title(sourceFilename);drawnow
    disp(sprintf('(motexWidefileStimuli) ImageNum: %i Converting from %s: %s',imageNum,getLastDir(matchTex.paths{iFolder}{fileNum}),getLastDir(sourceFilename)));
    % make gray image
    grayImage = sourceImage;
    grayImage(:) = grayValue;
    % get number of gray images
    nGrayImages = nImageRepeat;
    % if flickering, then double
    if flicker, nGrayImages = nGrayImages*2; end
    % write out the gray image
    for iImageRepeat = 1:nGrayImages
      % create destFilename
      destFilename = sprintf('%s%04i',destFilenameStem,imageNum);
      destFilename = setext(destFilename,destType);
      destFilename = fullfile(destDir,destFilename);
      % write out what was done
      readmeText = sprintf('%s%04i\tgray\tgray\n',readmeText,imageNum);
      % update image num
      imageNum = imageNum+1;
      % write image
      imwrite(grayImage,destFilename,destType,'Compression','lzw');
    end
    % write out the image
    for iImageRepeat = 1:nImageRepeat
      % create destFilename
      destFilename = sprintf('%s%04i',destFilenameStem,imageNum);
      destFilename = setext(destFilename,destType);
      destFilename = fullfile(destDir,destFilename);
      % write out what was done
      readmeText = sprintf('%s%04i\t%s\t%s\n',readmeText,imageNum,getLastDir(matchTex.paths{iFolder}{fileNum}),getLastDir(matchTex.filenames{iFolder}{fileNum}));
      % update image num
      imageNum = imageNum+1;
      % write image
      imwrite(sourceImage,destFilename,destType,'Compression','lzw');
      % if flicker, then put in gray image, after image stimulus image
      if flicker
	% create destFilename
	destFilename = sprintf('%s%04i',destFilenameStem,imageNum);
	destFilename = setext(destFilename,destType);
	destFilename = fullfile(destDir,destFilename);
	% write out what was done
	readmeText = sprintf('%s%04i\tgray\tgray\n',readmeText,imageNum);
	% update image num
	imageNum = imageNum+1;
	% write image
	imwrite(grayImage,destFilename,destType,'Compression','lzw');
      end
    end
  end
end

% write out readme file
f = fopen(fullfile(destDir,readmeFilename),'w');
fprintf(f,readmeText);
fclose(f);

% save stimulus information
save(fullfile(destDir,stimulusInfoFilename),'e');


%%%%%%%%%%%%%%%%%%%%%%%%%
%    getAllTexFamily    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function texFamily = getAllTexFamily(sourceDir,searchFolders)

% clear texFamily
texFamily = [];
% display what we are doing
dispHeader('(motexWidefieldStimuli) Finding all texture types');
for iFolder = 1:length(searchFolders)
  % init texFamily for this folder
  thisTexFamily = {};
  allFiles = dir(fullfile(sourceDir,searchFolders{iFolder}));
  for iFile = 1:length(allFiles)
    % make it sure looks like an image file with an extension
    if ~isempty(getext(allFiles(iFile).name))
      % grab the name of the texture class
      className = strsplit(stripext(allFiles(iFile).name),'_');
      className = className{end-1};
      % add it to list of textures
      thisTexFamily = union(thisTexFamily,className);
    end
  end
  % now keep all of the textures found
  if iFolder == 1
    texFamily = thisTexFamily;
  else
    % or compare to make sure that we have matching textures in all directoris
    if (length(thisTexFamily) ~= length(texFamily)) || any(~ismember(thisTexFamily,texFamily))
      disp(sprintf('(motexWidefieldStimuli) Folder %s has mismateched texture families',searchFolders{iFolder}))
      % display them
      texFamily
      thisTexFamily
      keyboard
    end
  end
end

% display what we found
disp(sprintf('(motexWidefieldStimuli) Found %i texture famililes',length(texFamily)));
disp(texFamily);


%%%%%%%%%%%%%%%%%%%%%%%%
%    getAllFieldNum    %
%%%%%%%%%%%%%%%%%%%%%%%%
function allMatch = getAllFieldNum(sourceDir,searchFolders,fieldNum)

% function gets all of the files with the k from last field in file name
% e.g. filenames are like noise_1x1_pool1_beans_smp1.png, so k = 2 returns
% the texFamily and k = 3 the algorithm


% clear 
allMatch = [];
for iFolder = 1:length(searchFolders)
  % init thisMatch for this folder
  thisMatch = {};
  % list of files
  allFiles = dir(fullfile(sourceDir,searchFolders{iFolder}));
  for iFile = 1:length(allFiles)
    % make it sure looks like an image file with an extension
    if ~isempty(getext(allFiles(iFile).name))
      % grab the name of the texture class
      thisName = strsplit(stripext(allFiles(iFile).name),'_');
      thisName = thisName{end-fieldNum+1};
      % add it to list
      thisMatch = union(thisMatch,thisName);
    end
  end
  % now keep all that we found
  if iFolder == 1
    allMatch = thisMatch;
  else
    % or compare to make sure that we have matching textures in all directoris
    if (length(thisMatch) ~= length(allMatch)) || any(~ismember(thisMatch,allMatch))
      disp(sprintf('(motexWidefieldStimuli) Folder %s has mismatch',searchFolders{iFolder}))
      % display them
      allMatch
      thisMatch
      keyboard
    end
  end
end

% display what we found
disp(sprintf('(motexWidefieldStimuli) Found %i matches',length(allMatch)));
disp(allMatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getMatchTextures    %
%%%%%%%%%%%%%%%%%%%%%%%%%%

function matchTex = getMatchTextures(sourceDir,searchFolders,texGenTypes,texFamily)

% get the names of all textures to convert
for iFolder = 1:length(searchFolders)
  % init list of files to copy
  matchTex.filenames{iFolder} = {};
  matchTex.paths{iFolder} = {};
  matchTex.genType{iFolder} = {};
  matchTex.family{iFolder} = {};
  matchTex.sample{iFolder} = {};
  matchTex.sampleNum{iFolder} = [];

  dispHeader(sprintf('Searching: %s',searchFolders{iFolder}));
  for iTexGenType = 1:length(texGenTypes)
    dispHeader(sprintf('Texure type: %s',texGenTypes{iTexGenType}));
    for iTexFamily = 1:length(texFamily)
      dispHeader(sprintf('Texure class: %s',texFamily{iTexFamily}));
      % create search name
      searchName = fullfile(sourceDir,searchFolders{iFolder},sprintf('*_%s_%s_smp*.*',texGenTypes{iTexGenType},texFamily{iTexFamily}));
      % look for match
      thisMatchFilenames = dir(searchName);
      % no match
      if isempty(thisMatchFilenames)
	disp(sprintf('(motexWidefieldStimuli) Could not find any files: %s',searchName));
	return
      end
      % add to list of files to convert
      for i = 1:length(thisMatchFilenames)
	% add to lists
	matchTex.paths{iFolder} = {matchTex.paths{iFolder}{:} fullfile(sourceDir,searchFolders{iFolder})};
	matchTex.filenames{iFolder} = {matchTex.filenames{iFolder}{:} thisMatchFilenames(i).name};
	matchTex.genType{iFolder} = {matchTex.genType{iFolder}{:} texGenTypes{iTexGenType}};
	matchTex.family{iFolder} = {matchTex.family{iFolder}{:} texFamily{iTexFamily}};
	% get the sample number
	sampleString = stripext(thisMatchFilenames(i).name);
	while ~isempty(sampleString)
	  [sampleNum sampleString] = strtok(sampleString,'_');
	end
	numloc = regexp(sampleNum,'[0-9]');
	matchTex.sample{iFolder} = {matchTex.sample{iFolder}{:} sampleNum};
	if ~isempty(numloc)
	  sampleNum = str2num(sampleNum(numloc(1):end));
	  matchTex.sampleNum{iFolder}(end+1) = sampleNum;
	else
	  matchTex.sampleNum{iFolder}(end+1) = nan;
	end
	% display filename
	disp(sprintf('%s',fullfile(getLastDir(matchTex.paths{iFolder}{end}),matchTex.filenames{iFolder}{end})));
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%
%    makeDestDir    %
%%%%%%%%%%%%%%%%%%%%%
function makeDestDir(searchFolders,destDir)

% check old directory
for iFolder = 1:length(searchFolders)
  if isdir(destDir) && ~isempty(ls(destDir))
    % remove existing files
    dispHeader(destDir);
    disp(sprintf('(motexWidefieldStimuli) Removing old files'));
    if 1%askuser(sprintf('(motexWidefiledStimuli) Ok to remove old files from %s',destDir))
      delete(fullfile(destDir,'*'));
    end
  else
    % otherwise make the directory
    if ~isdir(destDir)
      mkdir(destDir);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getRandomizeOrder    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulusOrder = getRandomizeOrder(matchTex,randomizeOrder)

nFiles = length(matchTex.filenames{1});

% get the randomizeOrder
if randomizeOrder
  % get a random stimulus order
  stimulusOrder = randperm(nFiles);
  % just in case there is no possible ordering
  tryNumber = 1;maxTryNumber = 10000;
  % see if any texture family repeats
  while any(strcmp({matchTex.family{1}{stimulusOrder(2:end)}},{matchTex.family{1}{stimulusOrder(1:end-1)}})) && (tryNumber < maxTryNumber)
    % get a new stimulus ordering
    stimulusOrder = randperm(length(matchTex.filenames{1}));
    tryNumber = tryNumber+1;
  end
  if tryNumber > maxTryNumber
    disp(sprintf('(motexWidefieldStimuli) Could not find a random stimulus ordering which does not repeat classes after %i tries',tryNumber));
  else
    disp(sprintf('(motexWidefieldStimuli) Found random stimulus ordering which does not repeat classes after %i tries',tryNumber));
  end
else
  % do in order
  stimulusOrder = 1:nFiles;
end

