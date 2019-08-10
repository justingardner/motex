% motexCamera2nifti.m
%
%      usage: [data hdr] = motexCamera2nifti(cameraFilename)
%         by: justin gardner
%       date: 07/19/19
%    purpose: loads a camera file and converts using mlrImage to
%             nifti compatabile data format - actually a data / hdr that 
%             mlrImage functions can manipulate and save as nifti
%
function [data hdr] = motexCamera2nifti(cameraFilename,varargin)

data = [];hdr = [];

% check arguments
if nargin < 1
  help motexCamera2nifti
  return
end

% check arguments
getArgs(varargin,{'autoCorrect=0','spoof3d=0'});

% check for file
if ~isfile(cameraFilename)
  disp(sprintf('(motexCamera2nifti) Could not find file %s',cameraFilename));
  return
end

% Use Yuki's function to read PCO file
[nRows, nCols, timeStamps, data, startTime,nanFrames]= loadPCOFile(cameraFilename, autoCorrect);

% see if we can find the CameraInfo file in the directory above
mmPerCameraPix = 1;
cameraInfoFilename = fullfile(fileparts(cameraFilename),'..','CameraInfo.mat');
if isfile(cameraInfoFilename)
  load(cameraInfoFilename)
end
if ~exist('CameraInfo','var')
  disp(sprintf('(motexCamera2nifti) Could not find CameraInfo for %s so mm dimension unknwon',getLastDir(cameraFilename)));
  CameraInfo = [];
else
  if isfield(CameraInfo,'MmPerCameraPix')
    mmPerCameraPix = CameraInfo.MmPerCameraPix;
  end
end

% get time per image
deltaT = median(diff(timeStamps));

% make a header
hdr.nDim = length(size(data));
hdr.pixdim = [mmPerCameraPix mmPerCameraPix deltaT];
hdr.dim = size(data);
hdr.qform = diag([mmPerCameraPix mmPerCameraPix 0 deltaT]);
hdr.sform = hdr.qform;
hdr.filename = cameraFilename;
hdr.startTime = startTime;
hdr.nanFrames = nanFrames;
hdr.cameraInfo = CameraInfo;
hdr.ext = 'PCO';
[tf hdr] = mlrImageIsHeader(hdr);

% if we need to spoof as a 3d volume 
if spoof3d
  % spoof a 3rd dimension
  data = reshape(data,hdr.dim(1),hdr.dim(2),1,hdr.dim(3));
  hdr.dim = [hdr.dim(1) hdr.dim(2) 1 hdr.dim(3)];
  hdr.nDim = 4;
  hdr.pixdim = [hdr.pixdim(1) hdr.pixdim(2) 1 hdr.pixdim(3)];
end

