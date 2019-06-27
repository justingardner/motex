function CameraInfo = GetCameraInfo(animal,iseries,rewriteflag)
% GetCameraInfo prompts the user for basic camera info
%
% CameraInfo = GetCameraInfo(animal,iseries)
%
% CameraInfo = GetCameraInfo(animal,iseries,1) forces a rewrite
%
% Example: 
% CameraInfo = GetCameraInfo('catz057',1)
%
% 2006-01 Matteo Carandini

if nargin < 3
    rewriteflag = 0;
end

global DIRS

if isempty(DIRS) || ~isfield(DIRS,'camera')
    error('Should have declared a global DIRS with the camera files');
end
        
SeriesDir = fullfile(DIRS.camera,animal,num2str(iseries));

if ~exist(SeriesDir,'dir')
    error('Could not find series directory');
end

CameraFile = fullfile(SeriesDir,'CameraInfo.mat');

if exist(CameraFile,'file') && ~rewriteflag
    
    load(CameraFile);

else
    
   prompt = {'Camera model','Top lens focal length (mm)','Bottom lens focal length (mm)','Hardware binning factor'};
   name = sprintf('Camera Details for %s series %d',animal,iseries);
   numlines = 1;
   defaultanswer = {'PCO edge','50','50','1'};
 
   answer=inputdlg(prompt,name,numlines,defaultanswer);
 
   CameraInfo.Model             = answer{1};
   CameraInfo.TopFocalLength    = str2double(answer{2});
   CameraInfo.BotFocalLength    = str2double(answer{3});
   CameraInfo.HardwareBinning   = str2double(answer{4});
   
   switch CameraInfo.Model
       case 'Dalsa M160'
           CameraInfo.MmPerCameraPix = 0.014 * CameraInfo.HardwareBinning;
       case 'Photonfocus'
           prompt = {'Zoom gain'};
           name = sprintf('What is the zoom gain for this photonfocus?'); %ask AP for details
           numlines = 1;
           defaultanswer = {'1'};
           answer = inputdlg(prompt,name,numlines,defaultanswer);
           zoom_gain = str2double(answer{1});
           
           CameraInfo.MmPerCameraPix = 0.0105 * CameraInfo.HardwareBinning*zoom_gain;
       case 'PCO edge'
           CameraInfo.MmPerCameraPix = 0.0063 * CameraInfo.HardwareBinning; %AP Oct 11
       otherwise
           error('Do not know this camera type');
   end
   
   save(CameraFile,'CameraInfo');
end
