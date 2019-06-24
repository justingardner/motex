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

getArgs(varargin,{'sourceDir=/Volumes/GoogleDrive/My Drive/docs/2019/NSF CRCNS/Stimuli','texFolder=tex_eq','noiseFolder=noise_eq','destDir=~/Desktop/expt','texFamily',{'glass','fronds','spikes','beans','crystals','rocks','scales','clouds'},'texGenTypes',{'pool4','PS'},'destType=tif','nImageRepeat=5','destFilenameStem=im','grayValue=128','readmeFilename=readme.txt','stimulusInfoFilename=stimulusInfo','trialTime=10','outputSize=255','flicker=0','iti=5','imageFrequency=2.0','randomizeOrder=1'});

% get files from both tex and noise
searchFolders = {texFolder, noiseFolder};

% get the names of all textures to convert
for iFolder = 1:length(searchFolders)
  % init list of files to copy
  matchFilenames{iFolder} = {};
  matchPaths{iFolder} = {};
  matchTexGenType{iFolder} = {};
  matchTexFamily{iFolder} = {};
  matchTexSample{iFolder} = {};
  matchTexSampleNum{iFolder} = [];

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
	matchPaths{iFolder} = {matchPaths{iFolder}{:} fullfile(sourceDir,searchFolders{iFolder})};
	matchFilenames{iFolder} = {matchFilenames{iFolder}{:} thisMatchFilenames(i).name};
	matchTexGenType{iFolder} = {matchTexGenType{iFolder}{:} texGenTypes{iTexGenType}};
	matchTexFamily{iFolder} = {matchTexFamily{iFolder}{:} texFamily{iTexFamily}};
	% get the sample number
	sampleString = stripext(thisMatchFilenames(i).name);
	while ~isempty(sampleString)
	  [sampleNum sampleString] = strtok(sampleString,'_');
	end
	numloc = regexp(sampleNum,'[0-9]');
	matchTexSample{iFolder} = {matchTexSample{iFolder}{:} sampleNum};
	if ~isempty(numloc)
	  sampleNum = str2num(sampleNum(numloc(1):end));
	  matchTexSampleNum{iFolder}(end+1) = sampleNum;
	else
	  matchTexSampleNum{iFolder}(end+1) = nan;
	end
	% display filename
	disp(sprintf('%s',fullfile(getLastDir(matchPaths{iFolder}{end}),matchFilenames{iFolder}{end})));
      end
    end
  end
end

% calculate how long sequence will take
nFiles = length(matchFilenames{1});
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

% readme string
readmeText = '';

% get the randomizeOrder
if randomizeOrder
  % get a random stimulus order
  stimulusOrder = randperm(nFiles);
  % just in case there is no possible ordering
  tryNumber = 1;maxTryNumber = 10000;
  % see if any texture family repeats
  while any(strcmp({matchTexFamily{1}{stimulusOrder(2:end)}},{matchTexFamily{1}{stimulusOrder(1:end-1)}})) && (tryNumber < maxTryNumber)
    % get a new stimulus ordering
    stimulusOrder = randperm(length(matchFilenames{1}));
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

% write out image sequence
imageNum = 1;iTrial = 0;
for iFile= 1:length(matchFilenames{1})
  % get filenumber
  fileNum = stimulusOrder(iFile);
  % keep information about image
  iTrial = iTrial + 1;
  e.texFamily{iTrial} = matchTexFamily{1}{fileNum};
  e.texGenType{iTrial} = matchTexGenType{1}{fileNum};
  e.texSampleNum(iTrial) = matchTexSampleNum{1}(fileNum);
  % create images
  for iFolder = 1:length(searchFolders)
    % get source filename
    sourceFilename = fullfile(matchPaths{iFolder}{fileNum},matchFilenames{iFolder}{fileNum});
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
    disp(sprintf('(motexWidefileStimuli) ImageNum: %i Converting from %s: %s',imageNum,getLastDir(matchPaths{iFolder}{fileNum}),getLastDir(sourceFilename)));
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
      readmeText = sprintf('%s%04i\t%s\t%s\n',readmeText,imageNum,getLastDir(matchPaths{iFolder}{fileNum}),getLastDir(matchFilenames{iFolder}{fileNum}));
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

