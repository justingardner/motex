function func_DefaultDirs (flag, anim, seri, serv)
% DefaultDirs sets the default directories MA 20150123
% Directory of Files on Labserver
global DIRS

if nargin <4, serv = []; end
if serv ==1,  serv = []; end %onlt server 2 is tagged with #2

%% directory for pco and et data
%if serv==6
%    Server  = ['\\NCB-Labserver' num2str(serv) '\data\'];
Server  = fullfile(DIRS.labserver{serv});
%else
%    Server  = ['\\Labserver' num2str(serv) '\data\'];
%end
switch flag
    case 'YC'
        DIRS.camera    = fullfile(Server,'MOUSE','IMAGING','YC');
    case 'GCAMP'
        DIRS.camera    = fullfile(Server,'MOUSE','IMAGING','GCAMP');
        DIRS.ET        = fullfile(Server,'MOUSE','EyeTracking','GCAMP');
    case 'OBOX'
        DIRS.camera    = fullfile(Server,'MOUSE','IMAGING','OBOX');
        DIRS.ET        = fullfile(Server,'MOUSE','EyeTracking','OBOX');
end
DIRS.processeddata  = fullfile(Server,'MOUSE','PROCESSED DATA');
DIRS.behavior       = fullfile(Server,'OHARA_DATA');

%% Log Directories
DIRS.MPEP   = fullfile(DIRS.labserver{1},'MOUSE','LOGS','MPEP_LOGS'); %'\\Labserver\data\MOUSE\LOGS\MPEP_LOGS\'
DIRS.data   = fullfile(DIRS.labserver{1},'MOUSE','LOGS','MPEP_LOGS');
DIRS.BVS    = fullfile(DIRS.labserver{1},'OHARA_DATA','IMAGING');
DIRS.VS     = fullfile(DIRS.labserver{1},'MOUSE','LOGS','VS_LOGS',anim,[anim,'_',num2str(seri)]);

%% Directory for Segmentation result
DIRS.segres = fullfile(Server,'MOUSE','Segmentation','animalseg');

%% Save (and load) directory for analyzed mat files and figures
DIRS.SaveDir = fullfile(DIRS.camera,anim,'ANALYZED',num2str(seri));
if ~exist(DIRS.SaveDir, 'dir')
    disp ('No Save Directory. mkdir ...'); 
    mkdir(DIRS.SaveDir), 
end


%% LOG of this m file
%{
Mohammad 20170707 changed the DIRS.data to DIRS.MPEP and DIRS.ScreenInfo to DIRS.VS
%}

