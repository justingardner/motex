% Visual Area Segmentation
%{
Background science: 
1- Kalatsky and Stryker 2003 Neuron; 
2- Marshel et al 2011 Neuron
3- Garret et al 2014 J Neurosci
4- Zhuang et al Waters 2017 eLife

How to use:
1- input the experiment information in the dialog box
2- func_DefaultDirs sets the directories based on the input above. You may need to adjust this on your own.
3- GetCameraInfo gets the camera information used for experiments. If you use a different camera, you need to add it to GetCameraInfo
4- ProtocolLoad loads the experimental protocol used for experiment, which is saved to mpep logs (DIRS.MPEP)
5- load([DIRS.VS '_' num2str(iexp)]) loads the ScreenInfo calss object. Direct loading will not work. It is used to define visual_field.
6- Create the list of trials filenames to be loaded by loadWFtrials
7- loadWFtrials is a complicated function for many different experiment types. Use it with the given params
8- params = ('singletrials',false,'repFlag',false,'resFac',1,'LoCutFreq',0.05,'HiCutFreq',2,'spatfiltwid',3,'frame0list',1:2,'moco_flag', true);
    singletrials: keep single trial data
    repFlag: do retinotopy for every repeat
    resFac: resize factor for the tensor [0 1]
    LoCutFreq and HiCutFreq: the cutoff frequencies for temporal filtering
    spatfiltwid: the width of gaussian filter for spatial smoothing
    frame0list: the list of the frames for frame-0 correction
    moco_flag: motion correction
9- All the neural and experimental data are stored in VDAQ that is a global variable
10- Use processRetVDAQ function to apply all basic pre-processing (frame trimming, frame0-blank corrections, spatiotemporal filtering)
11- Use TensorFrequency to compute the complex maps that contain phase and amplitude maps
12- Save the complex maps into specified directory for later use
13- func_getRetinoMaps plots the azimuth and elevation maps and the visual field sign maps and the thresholded sign maps
14- Create & save sessregim to be used for registering future experiments to this retinotopy blod vessel image
15- Check the power spectra for more insight into the quality of this experiment using func_plot_VDAQFourietAmpSpect
16- bpViewComplexMaps plots the individual phase maps (forward and reverse), amplitude maps and their combination


- Mohammd Abdolrahmani 20190626: Created an independent toolbox for retinotopy, to be shared.
- For questions, contact:  mohammad@benuccilab.net; m.abd.ac@gmail.com

History of the toolbox:
1- Some of the files used in this toolbox are inherited from Andrea Benucci from Carandini lab. For instance bpViewComplexMaps, bpMakeJointMap, ...
2- Mohammad Abdolrahmani, in Benucci lab wrote the main_Retinotopy and its dependencies.
3- Mike Morais organized Mohammad's previous WF analysis code (e.g. mainNeural, mainWidefield,...) into the loadWFtrials, with additional functioalities
3- Mohammad made major changes to loadWFtrials, e.g. changed the order of processing steps, added plotting, removed redundant parts
%}

%clearvars; clear global

%% 01 // Declarations and setting required directories
global DIRS VDAQ

DIRS.labserver{1} = '/Volumes/DATA';
DIRS.labserver{2} = '/Volumes/DATA-1';
DIRS.labserver{3} = '/Volumes/DATA-1';
DIRS.labserver{4} = '/Volumes/DATA-1';
DIRS.labserver{5} = '/Volumes/DATA-2';
DIRS.labserver{6} = '/Volumes/DATA-1';

% MODE 1 -- single trial retinotopy
params = struct('singletrials',false,'repFlag',false,'resFac',1,'LoCutFreq',0.05,'HiCutFreq',2,'spatfiltwid',3,'frame0list',1:2,'moco_flag', true);
params.guisave = true;
guiSaveDir = '\\Labserver5\data\MOUSE\Segmentation\guidata';
params.regsave = false;
regSaveDir = '\\Labserver5\data\MOUSE\Segmentation\animalseg';

% prompt user for session information (without hard coding it!)
prompt = {'Indicator:','Animal/date:','ID: (only if saving reg files)','Series number:','Experiments: (enter as "# spacebar #")','Server:'};
def    = {'GCAMP','M090621_MA','18273','2','1 2','6'};
answer = inputdlg(prompt,'Enter session info...',1,def);
[indicator,animal,idnum,iseries,Expts,Serv] = deal(answer{:});
iseries = str2double(iseries);
Expts   = cellfun(@str2double,regexp(Expts,'[0-9]*','match'));
Serv    = str2double(Serv);

func_DefaultDirs(indicator, animal, iseries, Serv);
CameraInfo = GetCameraInfo(animal, iseries);

%% 02 // load neural data
for idx  = 1:length(Expts)
    
    % prompt user for which VDAQ to load, or to make a new one
    canddir = fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(Expts(1,idx)));
    candfiles = dir(fullfile(canddir,'VDAQ*.mat')); candfiles = {candfiles.name}; candfiles = cat(2,{'[MAKE NEW VDAQ]'},candfiles);
    fileselect = listdlg('PromptString','Select VDAQ loading method:','SelectionMode','single','ListString',candfiles,'Name',canddir);
    if fileselect == 1
        loadfile = [];
        savefile = inputdlg('Define a savename for the new VDAQ file:','',1,{'VDAQ.mat'}); savefile = savefile{1};
    else
        loadfile = fullfile(canddir,candfiles{fileselect});
    end
    
    % get the experiment number
    iexp  = Expts(1,idx);
    fprintf('starting expt %d of %d...\n',idx,length(Expts));
    % load associated data (params)
    p = ProtocolLoad(animal,iseries,iexp);
    if isempty(p.blankstims), p.blankstims = p.nstim; end
    load([DIRS.VS '_' num2str(iexp)]);  % ScreenInfo calss object
    % get data (neural images) filenames
    fileDir = fullfile(DIRS.camera,animal,num2str(iseries),num2str(iexp));
    triallist = dir(fullfile(fileDir,'*.mat')); triallist = {triallist.name};
    triallist = sort_nat(triallist);
    
    %% loop loading trials data
    if isempty(loadfile)
        %% load and populate VDAQ global variable (if tensor doesn't already exist) initialize VDAQ
        fprintf('loading trials and building VDAQ...\n');
        VDAQ = struct('animal',animal,'iseries',iseries,'iexp',iexp,'Frame0List',...
            params.frame0list,'ResizeFactor',params.resFac,'durs',p.pfiledurs(1),...
            'nstim',p.nstim,'fileList',{triallist},'MeanIntensities',[],...
            'MmPerCameraPix',CameraInfo.MmPerCameraPix*CameraInfo.BotFocalLength/CameraInfo.TopFocalLength/params.resFac);
        
        % load and average tensors (this new block of code replaces Main_WideField        ;
        if ~params.singletrials % averaged across trials
            for k = 1:p.nstim
                fprintf('    stimulus %d of %d\n',k,p.nstim);
                W = loadWFtrials( fileDir, triallist(p.seqnums(k,:)),[], 'alltrials_flag',true, ...
                    'downsamp',params.resFac, 'PreFrm',params.frame0list, 'trialavg_flag',true, ...
                    'LoCutFreq',params.LoCutFreq, 'HiCutFreq',params.HiCutFreq, 'trimendframe',0,...
                    'moco_flag',false, 'filter_flag',false );
                VDAQ.tensor{1,k} = W.xavg;
                VDAQ.MeanIntensities(k,:) = W.bg';
            end
        else % single trials
            for k = 1:p.nstim
                for kk = 1:size(p.seqnums,2)
                    fprintf('    stimulus %d of %d\n',k,p.nstim);
                    W = loadWFtrials( fileDir, triallist(p.seqnums(k,kk)),[], 'alltrials_flag',true, ...
                        'downsamp',params.resFac, 'PreFrm',params.frame0list, 'trialavg_flag',false, ...
                        'LoCutFreq',params.LoCutFreq, 'HiCutFreq',params.HiCutFreq, 'trimendframe',0,...
                        'moco_flag',false, 'filter_flag',false );
                    VDAQ.tensor(kk,k) = W.x;
                    VDAQ.MeanIntensities(k,kk) = W.bg;
                end
            end
        end
        
        % make all tensor units the same size
        fprintf('VDAQ structure postprocessing...\n');
        minsz = cellfun(@size,VDAQ.tensor,repmat({3},size(VDAQ.tensor)));
        S = struct('type','()','subs',cellfun(@transpose,cellfun(@squeeze,mat2cell(cat(3,repmat({':'},size(minsz)),...
            repmat({':'},size(minsz)),repmat({1:min(minsz)},size(minsz))),ones(size(minsz,1),1),ones(size(minsz,2),1),3),'uni',0),'uni',0));
        VDAQ.tensor = cellfun(@subsref,VDAQ.tensor,num2cell(S),'uni',0);
        VDAQ.nsummedframes = [];
        VDAQ.nrepeats = size(p.seqnums,2);
        VDAQ.meanvalue = nanmean(VDAQ.MeanIntensities(:));
        VDAQ.durs    = p.pfiledurs(1);
        VDAQ.FrameRate = median(reshape(cellfun(@size,VDAQ.tensor,repmat({3},size(VDAQ.tensor)))/p.pfiledurs(1),[],1));
        
        % build VDAQ global to match with the rest of the code
        [nr,nc,nt]   = size(VDAQ.tensor{1});
        VDAQ.durs    = nt/VDAQ.FrameRate;
        VDAQ.nframes = nt;
        VDAQ.ny      = nr;
        VDAQ.nx      = nc;
        VDAQ.tt      = linspace(0, VDAQ.durs, nt);
        VDAQ.PCOdata_transposed = true;
        VDAQ.BuildDate = datestr(now);
        p.pars(1,:)  = ones(size(p.pars(1,:))).*(VDAQ.durs*10);
        
        % deltaF/F normalization
        meanresps = zeros(size(VDAQ.tensor));
        for istim = 1:p.nstim
            for iten = 1:size(VDAQ.tensor,1)
                meanresps(iten,istim) = mean(VDAQ.tensor{iten,istim}(:));
            end
        end
        VDAQ.meanvalue = mean(meanresps(:));
        for istim = 1:p.nstim
            for iten = 1:size(VDAQ.tensor,1)
                VDAQ.tensor{iten,istim} = (VDAQ.tensor{iten,istim}/VDAQ.meanvalue)-1;
            end
        end
        
        % saving VDAQ
        if ~params.singletrials % only save trial-averaged data
            fprintf('saving VDAQ...\n');
            mkdir(fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(iexp)));
            % save the VDAQ file
            save(fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(iexp),savefile),'VDAQ','p','-v7.3')
            % save a registration file for behavior (if flagged)
            if idx == 1 && params.regsave
                matches = dir(fullfile(regSaveDir,sprintf('R%s_*',idnum))); matches = {matches.name};
                matches = regexp(matches,'[A-Z0-9]*\_([0-9]*)','tokens');
                if ~isempty(matches)
                    ridx = max(cellfun(@str2double,cellfun(@cell2mat,vertcat(matches{:}),'uni',0)));
                else
                    ridx = 0;
                end
                sessregim = nanmean(VDAQ.tensor{1},3);
                save(fullfile(regSaveDir,sprintf('R%s_%d.mat',idnum,ridx+1)),'sessregim');
            end
        else
            fprintf('cannot save this version VDAQ -- too big...\n');
        end
    else
        fprintf('VDAQ save already detected! automatically loading that for you... overwrite manually if you want to reanalyze\n')
        load(loadfile);
    end
    
    %% Peform all basic pre-processing (frame trimming, frame0-blank corrections, spatiotemporal filtering)
    [p,stimfreqs] = processRetVDAQ(p,params);
    
    %% Compute complex maps
    VDAQ_full = VDAQ;
    for k = 1:size(VDAQ_full.tensor,1)
        VDAQ.tensor = VDAQ_full.tensor(k,:);
        [AbsMaps{k,idx}, AngleMaps{k,idx}, CmplxMaps{k,idx}] = TensorFrequency([], stimfreqs ); % Get Azm and elv maps
    end
    
    %% Calculate visual field (screen x and y in degrees of visual angle)
    monitor_distance = myScreenInfo.Dist;
    if strcmp(p.xfile,'stimMarshel.x')
        if  p.pars(13,1)==1 % ori 90 or 0 (if 1 or 2)
            screen_fractionX = (p.pars(3,1)-p.pars(2,1))/360;
            L = myScreenInfo.Xmax*myScreenInfo.PixelSize*screen_fractionX; % in cm
        elseif p.pars(13,1)==2
            screen_fractionY = (p.pars(3,1)-p.pars(2,1))/360;
            L = myScreenInfo.Ymax*myScreenInfo.PixelSize*screen_fractionY; % in cm
        end
    elseif strcmp(p.xfile,'stimKalatsky.x')
        if  p.pars(end,1)==1 % ori 90 or 0 (if 1 or 2)
            screen_fractionX = (p.pars(3,1)-p.pars(2,1))/360;
            L = myScreenInfo.Xmax*myScreenInfo.PixelSize*screen_fractionX; % in cm
        elseif p.pars(end,1)==2
            screen_fractionY = (p.pars(3,1)-p.pars(2,1))/360;
            L = myScreenInfo.Ymax*myScreenInfo.PixelSize*screen_fractionY; % in cm
        end
    end
    visual_field(idx) = 2*atan(L/(2*monitor_distance))*180/pi;% in visual angle    
    
end

% Apply eye position (meridians) with respect to the screen
EyeY = 0;
EyeX = 0;

%% Save the file for segmentation analysis
fprintf('saving AnalyzedRet...\n');
save([DIRS.SaveDir 'AnalyzedRet'],'CmplxMaps','EyeX', 'EyeY','visual_field')
if params.guisave
    save(fullfile(guiSaveDir,sprintf('AnalyzedRet_%s-%s-%d',indicator,animal(1:end-3),iseries)),'CmplxMaps','EyeX', 'EyeY','visual_field')
end
fprintf('... complete ...\n');


%% Get all maps(Azm, Alt, VFS, VFS_thr)
[map_hor, map_vert, VFS, VFS_thr] = func_getRetinoMaps(CmplxMaps, visual_field, params.resFac);

%% (for registration) is not made, manually create & save sessregim to be used for registering future experiments to this image
dirPCO = fullfile(DIRS.camera, animal, num2str(iseries), num2str(iexp));
pcoFiles  = dir( fullfile(dirPCO,'*.mat' ));
pcoFiles  = struct2cell(pcoFiles);
pcoFilesS = sort_nat(pcoFiles(1,:)); % sort files in natural order
pcoFilesS = pcoFilesS (1:numel(pcoFilesS)); % now ignore the last trials that are corrupted
filename = fullfile(dirPCO, pcoFilesS{1});
[nRows, nCols, timeStamps, rawData, startTime] = loadPCOFile(filename);
rawData = (fliplr(imrotate(mean(rawData,3),-90)));
%sessregim = mean(rawData,3);
sessregim = rawData;

%Dir = '\\Labserver5\data\MOUSE\Segmentation\animalseg\';
%save([Dir 'R' idnum '_1'],'sessregim')
figure; imshow(sessregim(:,:,1),[]); axis image
hold on; imcontour(imresize(VFS_thr,size(sessregim),'nearest'),[1,1],'r')

%% check power spectra
for idx  = 1:length(Expts)
    iexp  = Expts(1,idx);
    canddir = fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(Expts(1,idx)));
    load(fullfile(canddir,'VDAQ.mat'))
    p = ProtocolLoad(animal,iseries,iexp);
    if idx==1, ROI = []; end
    [p,stimfreqs] = processRetVDAQ(p,params);
    [TimeSer{idx}, Frq{idx}, Amp{idx}, ROI] = func_plot_VDAQFourietAmpSpect(VDAQ,ROI);
end
% func_saveFig([],[],200,[],[idnum '_' animal '_' num2str(iseries) '_FourierAmp'],0,0,1)

%% view complex maps
tmp{1,1}=CmplxMaps{1}{1};
tmp{2,1}=CmplxMaps{1}{2};
tmp{1,2}=CmplxMaps{2}{1};
tmp{2,2}=CmplxMaps{2}{2};
bpViewComplexMaps(tmp) 
hold on; imcontour(VFS_thr,[1 1],'k')
clear tmp
% func_saveFig([],[],200,[],[idnum '_' animal '_' num2str(iseries) '_PhaseAmp'],0,0,1)



%% Log
%{
- Mohammad 2015: wrote the main files
- Mohammd 20190626: Put together as an independent toolbox for retinotopy.
%}
