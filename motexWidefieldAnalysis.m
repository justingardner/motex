% motexWidefieldAnalysis.m
%
%        $Id:$ 
%      usage: motexWidefieldAnalysis()
%         by: justin gardner
%       date: 06/23/19
%    purpose: 
%
function retval = motexWidefieldAnalysis(varargin)

getArgs(varargin,{'dataPath=/Volumes/GoogleDrive/My Drive/docs/2019/NSF CRCNS/data','dataDir=M190621_MA','processedDataPath=VDAQtensor/1/1','stimulusInfoPath=ANALYZED'});

% load processed data
d = motexLoadProcessedData(fullfile(dataPath,dataDir,processedDataPath));

% load stimulus info
d = motexLoadStimulusInfo(d,fullfile(dataPath,dataDir,stimulusInfoPath));;

% load anatomy image
d = motexLoadAnatomyImage(d,fullfile(dataPath,dataDir));

% analysize tseries
d = motexTSeriesAnalysis(d)

keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexTSeriesAnalysis    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = motexTSeriesAnalysis(d)

% create a mask
[d.mask.x d.mask.y] = meshgrid(80:100,40:80);
d.mask.i = sub2ind(d.dataSize(1:2),d.mask.x,d.mask.y);

% reshape data
data = reshape(d.data,d.nTrials,prod(d.dataSize(1:2)),d.dataSize(3));

% get average time series over mask
d.tSeries = squeeze(mean(data(:,d.mask.i,:),2));

% now find all of the data associated with good images
strongEffectMembers = ismember(d.texFamily,{'glass','fronds','spikes','beans'});
weakEffectMembers = ismember(d.texFamily,{'rocks','scales','clouds','crystals'});

% compute mean and standard error
meanTSeries = 100*mean(d.tSeries,1);
steTSeries = 100*std(d.tSeries,1)/sqrt(d.nTrials);

% compute for strong effect texture families
strongEffect = 100*mean(d.tSeries(strongEffectMembers,:),1);
strongEffectSTE = 100*std(d.tSeries(strongEffectMembers,:),1)/sqrt(sum(strongEffectMembers));

% compute for weak effect texture families
weakEffect = 100*mean(d.tSeries(weakEffectMembers,:),1);
weakEffectSTE = 100*std(d.tSeries(weakEffectMembers,:),1)/sqrt(sum(weakEffectMembers));

% compute time base
sampleRate = 10;
t = (1:length(meanTSeries))/sampleRate;

% display figure
mlrSmartfig('motexWidefieldAnalysis1','reuse');clf
subplot(2,1,1);
myerrorbar(t,meanTSeries,'yError',steTSeries,'yErrorBarType=fill','Color=k');
hold on;
vline([2.5 5 7.5]);
xlabel('Time (s)');
ylabel('delatF/f (%)');
%drawPublishAxis;
drawnow;

subplot(2,1,2);
myerrorbar(t,strongEffect,'yError',strongEffectSTE,'yErrorBarType=fill','Color=r');
xlabel('Time (s)');
ylabel('delatF/f (%)');
hold on
myerrorbar(t,weakEffect,'yError',weakEffectSTE,'yErrorBarType=fill','Color=b');
hold on;vline([2.5 5 7.5]);
mylegend({'Strong effect families','Weak effect families'},{'r','b'});

makeEqualYaxis(2,1,[1 2]);
drawPublishAxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexLoadAnatomyImage    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = motexLoadAnatomyImage(d,filename)

% look for a bmp file
thisdir = dir(fullfile(filename,'*.bmp'));

if length(thisdir) ~= 1
  disp(sprintf('(motexWidefieldAnalysis:motexLoadAnatomyImage) Found %i bmp files inf %s',length(thisdir),filename));
  if length(thisdir) < 1,return,end
end

% load anatomy image
d.anatomyImage = imread(fullfile(filename,thisdir(1).name),'bmp');

% scale image to size of data
d.anatomyImageScaled = motexScaleImage(d.anatomyImage,d.dataSize(1:2));


%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexScaleImage    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function im = motexScaleImage(im,outputSize)

imageSize = size(im);
[sourceMeshX sourceMeshY] = meshgrid((0:1/(imageSize(2)-1):1)-0.5,(0:1/(imageSize(1)-1):1)-0.5);
[outputMeshX outputMeshY] = meshgrid((0:1/(outputSize(2)-1):1)-0.5,(0:1/(outputSize(1)-1):1)-0.5);
im = interp2(sourceMeshX,sourceMeshY,double(im),outputMeshX,outputMeshY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexLoadStimulusInfo    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = motexLoadStimulusInfo(d,filename)

% check for file
filename = fullfile(filename,'stimulusInfo.mat');
if ~isfile(filename)
  disp(sprintf('(motexWidefieldAnalysis:motexLoadStimulusInfo) Could not find file: %s',filename));
  return
end

% otherwise load
d.stimulusInfo = load(filename);
if isfield(d.stimulusInfo,'e')
  d.stimulusInfo = d.stimulusInfo.e;
  d.stimulusInfo.nTypes = length(d.stimulusInfo.texFamily);
else
  disp(sprintf('(motexWidefieldAnalysis:motexLoadStimulusInfo) Could not find e filed within stimulusInfo file: %s',filename));
  return
end

% figure out what stimulus each trial belonds to
if isfield(d,'data')
  d.stimNum = mod(1:d.nTrials,d.stimulusInfo.nTypes);
  d.stimNum(d.stimNum == 0) = d.stimulusInfo.nTypes;
end

% get some info
d.texFamily = {d.stimulusInfo.texFamily{d.stimNum}};
d.texGenType = {d.stimulusInfo.texGenType{d.stimNum}};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    motexLoadProcessedData(filename)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = motexLoadProcessedData(filename);

d = [];

% check valid directory
if ~isdir(filename)
  disp(sprintf('(motexWidefiledAnalysis:motexLoadProcessedData) %s not found',filename));
  return
end

% look for mat files
d.filename = filename;
d.dir = dir(fullfile(filename,'*.mat'));

% load the files
disppercent(-inf,'(motexWidefiledAnalysis:motexLoadProcessedData) Loading processed data');
for iFile = 1:length(d.dir)
  % load matlab file
  loadedData = load(fullfile(filename,d.dir(iFile).name));
  % get dimensions of data
  if isfield(loadedData,'trialdata')
    d.data{iFile} = loadedData.trialdata;
    d.dataSize(iFile,:,:,:) = size(loadedData.trialdata);
  else
    disp(sprintf('(motexWidefiledAnalysis:motexLoadProcessedData) Could not find field trialdata in file %s',d.dir(iFile).name));
    keyboard
  end
  disppercent(iFile/length(d.dir));
end
disppercent(inf);

% check sizes
for iDim = 1:3
  uniqueDimsize =unique(d.dataSize(:,iDim));
  if length(uniqueDimsize) > 1
    disp(sprintf('(motexWidefiledAnalysis:motexLoadProcessedData) Dimension %i of data has %i values (%s). Cropping all to the value %i',iDim,length(uniqueDimsize),num2str(uniqueDimsize(:)'),min(uniqueDimsize)));
  end
  dataSize(iDim) = min(uniqueDimsize);
end

% crop data
disppercent(-inf,sprintf('(motexWidefiledAnalysis:motexLoadProcessedData) Cropping data to %i %i %i',dataSize(1),dataSize(2),dataSize(3)));
d.dataSize = dataSize;
d.nTrials = length(d.data);
for iTrial = 1:d.nTrials
  data(iTrial,1:dataSize(1),1:dataSize(2),1:dataSize(3)) = d.data{iTrial}(1:dataSize(1),1:dataSize(2),1:dataSize(3));
  disppercent(iTrial/d.nTrials);
end
disppercent(inf);
d.data = data;

