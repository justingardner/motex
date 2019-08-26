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
[nRows, nCols, timeStamps, data, startTime,nanFrames]= motexLoadPCOFile(cameraFilename, autoCorrect);

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
hdr.qform = diag([mmPerCameraPix mmPerCameraPix 1 deltaT]);
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

function [nRows, nCols, timeStamps, imgData, startTime,nanFrames]= motexLoadPCOFile(filename,autoCorrect)
%% LOADPCOFILE load data saved by PCO.edge camera.
% Slightly better performance than tools.LoadPCO function.
% @author yuki
% last update 2018-10-05

    if nargin<2, autoCorrect=1; end % auto-correct by default
    
    fid=fopen(filename,'r');
    if fid ~=-1
        A = fread(fid,3, '*uint32'); % image data size
        startTime=fread(fid,1,'*double'); % UNIX time stamp
        imgData=fread(fid,prod(A),'*uint16'); % vectorized image data   
        timeStamps=fread(fid,A(3),'*double'); % time stamps
        fclose(fid);          
        nRows=A(1); % image height
        nCols=A(2); % image width
        nFrames=A(3); % number of frames
        % reshape to a 3-D matrix
        imgData=reshape(imgData,[nRows nCols nFrames]);
        % frames dropped by intensity correction (2018-07-12@yuki)
        nanFrames=[];
        % apply some corrections if requested
        if autoCorrect
            % drop some frames
            medInterval = median(diff(timeStamps));
            % if there is long gap between first two frames, drop the 1st
            if diff(timeStamps(1:2)) > 2.0 *medInterval
                fprintf('Correcting first frame.\n');
                timeStamps = timeStamps(2:end)-timeStamps(2);
                imgData=imgData(:,:,2:end);                
            end
            % if there is zero gap between last two frames, drop the last
            if diff(timeStamps(end-1:end)) < 0.5 *medInterval
                fprintf('Correcting last frame.\n');
                imgData=imgData(:,:,1:end-1);
                timeStamps=timeStamps(1:end-1);
            end
            
            % drop low intensity frames (2018-07-12@yuki)
            if autoCorrect==2     
                fprintf('Checking intensity\n');
                cropped=imgData(1:10:end,1:10:end,:);
                medIntensity=median(cropped(:)); % median intensity 
                disp(medIntensity)
                figure; plot(squeeze(mean(mean(cropped))))
                for i=1:size(cropped,3)
                    frameR = cropped(:,:,i); % get cropped frame     
                    if mean(frameR(:)) <medIntensity*0.5 % less than half the median intensity value
                        imgData(:,:,i)=NaN*imgData(:,:,i); % set to NaN                      
                        nanFrames=[nanFrames;i]; %#ok<AGROW>                        
                        fprintf('dropping frame %d: LOW INTENSITY!\n',i);
                    end
                end
            end
            
        end % end of autocorrect
    else
        error('Unable to open %s!\n',filename);
    end
 % end of function
%% LOG
%{
*2018-10-05@yuki: improved memory usage.
*2018-07-12@yuki: low intensity autocorrection added.
*2015-01-29@yuki: creation based on [tools.LoadPCO]
%}
 %end of file