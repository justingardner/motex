% script to load preprocessed data for VS-type experiments, including
% suite2p-segmented ROIs, eye-tracking, behavioral output and their
% timestamps. compatible for general VS experiments including texture
% presentation.
% author:ryo
% last updated@190627

% preparation
global DIRS   %necessary for func_DefaultDirs.m
Serv = 6;
addpath('C:\Users\ryoao\Documents\SVN\LAB_CODE\MATLAB\MOUSE\2P\2PAnalysisToolbox2\npy-matlab'); % function to load .npy in matlab
addpath('C:\Users\ryoao\Documents\SVN\LAB_CODE\MATLAB\MOUSE\WF\BehaviorVS\bvsEvent'); % for func_wheelZeroCross.m

%% set variables
% animal = 'M190405_18226';
% iser = '2'
% animal = 'M190421_18228';
% iser = '1'
% animal = 'M190604_18272';
% iser = '1';
animal = 'M190625_18273';
iser = '1';
func_DefaultDirs('GCAMP', animal, iser, Serv, false); % necessary for ProtocolLoad.m
do_ET = false; % no eyetracking

filepath = fullfile('//NCB-LABSERVER6/data/MOUSE/IMAGING/GCAMP', animal, 'Preprocessed', sprintf('%s_%s.h5', animal, iser));

%% fetch eye tracking info
if do_ET
    % DeepLabCut model for eyetracking. If DLC_ver in .h5 file is different,
    % re-do the analysis
    %dlc_ver = 'DeepCut_resnet50_mouse_ET_MA190314shuffle1_80000';
    dlc_ver = 'DeepCut_resnet50_mouse_ET_MA190613shuffle1_550000';
    info = h5info(filepath);
    attr_names = {info.Attributes.Name};
    if any(cellfun(@(x) strcmp(x, 'dlc_ver'), attr_names)) && strcmp(dlc_ver, h5readatt(filepath, '/','dlc_ver'))
        error('eyetracking already done!');
    end
    
    % update dlc_ver (attribute)
    h5writeatt(filepath, '/', 'dlc_ver', dlc_ver);
end

%% load: common across iexp
% expSet
expSet = h5read(filepath,sprintf('/expSet'));

% expType
expType = {};
for iexp = expSet'
    expType = [expType; h5read(filepath, sprintf('/%d/expType',iexp))];
end

% number of trials
ntri = h5read(filepath, sprintf('/ntri'));

% suite2p output
twoP_F = readNPY(fullfile('\\NCB-LABSERVER6\data\MOUSE\IMAGING\GCAMP', animal, '\2P\', [animal '_' num2str(iser)], '\suite2p\plane0\F.npy'));
twoP_Fneu = readNPY(fullfile('\\NCB-LABSERVER6\data\MOUSE\IMAGING\GCAMP', animal, '\2P\', [animal '_' num2str(iser)], '\suite2p\plane0\Fneu.npy')); % neuropil
twoP_spks = readNPY(fullfile('\\NCB-LABSERVER6\data\MOUSE\IMAGING\GCAMP', animal, '\2P\', [animal '_' num2str(iser)], '\suite2p\plane0\spks.npy')); % deconvolved spikes

% load twoP_iscell and save in h5
twoP_iscell = readNPY(fullfile('\\NCB-LABSERVER6\data\MOUSE\IMAGING\GCAMP', animal, '\2P\', [animal '_' num2str(iser)], '\suite2p\plane0\iscell.npy')); % 1st col: manual selection, 2nd col: reliabiity
info = h5info(filepath);
dataset_names = {info.Datasets.Name};
if ~any(cellfun(@(x) strcmp(x, 'twoP_iscell'), dataset_names))
    h5create(filepath, '/twoP_iscell', size(twoP_iscell));
end
h5write(filepath, '/twoP_iscell',twoP_iscell);

twoP_n_frame = h5read(filepath,sprintf('/twoP_n_frame'));

% driver cell index
if any(cellfun(@(x) strcmp(x, 'dc_idx'), dataset_names))
    dc_idx = h5read(filepath, sprintf('/dc_idx'));
end

%% load: loop for each iexp
for iexp = 1:numel(expSet)
    exp_this = expSet(iexp);
    
    % wheel timestamp
    load(fullfile('//Labserver/data/MOUSE/LOGS/VS_LOGS', sprintf('%s_%s_%d.mat', animal, iser, exp_this))); % load AIdata
    wheeltimestamp = AIdata(:,[1,11]); % channel 11 (CTRL 1) corresponds to wheel
    mosaicTTL = AIdata(:,7); % on for MOSAIC stim period
    mosaicTiming = [wheeltimestamp(find(diff(mosaicTTL) > 2)+1,1), wheeltimestamp(find(diff(mosaicTTL) < -2)+1,1)];
    %    wheelvel = ([0; diff(smooth(wheeltimestamp(:,2),330))]/median(diff(wheeltimestamp(:,1))));
    wheelvel = ([0; diff(wheeltimestamp(:,2))]/median(diff(wheeltimestamp(:,1))));
    
    % wheel onset detection
    wheel_dsmp = wheeltimestamp(1:1000:end,:); % downsampled wheel to 10 Hz
    wheelvel_dsmp = [0; diff(wheel_dsmp(:,2))]./median(diff(wheel_dsmp(:,1)));
    [signchT,signch,mvtsize] = func_wheelZeroCross([wheel_dsmp(:,2),wheel_dsmp(:,1)],[wheelvel_dsmp, wheel_dsmp(:,1)],5,10); % threshold was determined empirically for detection of non-trained wheel movement
    %    [signchT,signch,mvtsize] = func_wheelZeroCross([],[wheelvel_dsmp, wheel_dsmp(:,1)],10,10); % threshold was determined empirically for detection of non-trained wheel movement
    title(sprintf('wheel position and velcicty, exp: %d', exp_this));
    
    % photo diode
    photodiode_pulse_rise = AIdata(find(diff(AIdata(:,2))>2)+1,1); % on during stim presentation but flicker due to 60 Hz refresh of monitor
    photodiode_onset = photodiode_pulse_rise([1; find(diff(photodiode_pulse_rise) > 2*0.016)+1]); % > twice longer than monitor refresh is considered as stimulus onset
    photodiode_pulse_fall = AIdata(find(diff(AIdata(:,2))<-2)+1,1);
    photodiode_offset = [photodiode_pulse_fall(find(diff(photodiode_pulse_fall) > 2*0.016)); photodiode_pulse_fall(end)];
    
    %% ET Deeplabcut
    if do_ET
        ET_timestamp = h5read(filepath, sprintf('/%d/ET_timestamp',exp_this)); % ET frame onset
        
        % need sanity check that csv file have same #frame with original data!!
        mp4_path = '\\Ncb-labserver6\data\MOUSE\IMAGING\GCaMP\';
        [pos, area, sacT, not_fitted_frms, blink] = deal(cell(1));
        for itr = 1:ntri(iexp)
            eye_file = [DIRS.ET animal '\' num2str(iser) '\' num2str(exp_this) '\' animal '_' num2str(itr) '.eye'];
            csv_file = fullfile(mp4_path, animal, 'ET', num2str(iser), [animal '_' num2str(exp_this) '_' num2str(itr) dlc_ver '.csv']);
            [pos{itr}, area{itr}, sacT{itr}, not_fitted_frms{itr}] =  extractEyeInfo(eye_file, csv_file, false, false,false);
            area{itr} = medfilt1(area{itr},10); % smoothing, on the assumpution that area changes slowly
        end
                
        % semi-manual annotation for blink
        fig_handle = figure('position',[156 121 1023 819]);
        for itr = 1:ntri(iexp)
            blink_thistrial = [];
            not_fitted_thistrial = not_fitted_frms{itr};
            for iFr = find(diff(not_fitted_thistrial)==-1)'+1 % scan for all fall point
                if iFr <= size(not_fitted_thistrial,1)-2 && any(not_fitted_thistrial(iFr+1:iFr+2)) % replace 0 with 1 if rise again in 1 or 2 frames
                    not_fitted_thistrial(iFr:iFr+1) = true;
                end
            end
            
            idx_blink_onset = find(diff(not_fitted_thistrial)==1)+1;
            blink_offsets = [find(diff(not_fitted_thistrial)==-1); size(not_fitted_thistrial,1)]; % falls or last frame
            for iFr = idx_blink_onset'
                if iFr <= size(not_fitted_thistrial,1)-2 && all(not_fitted_thistrial(iFr:iFr+2)) % if not-fitted for at least 3 consective frames
                    idx_blink_offset = blink_offsets(find(blink_offsets>iFr, 1, 'first'));
                    eye_file = [DIRS.ET animal '\' num2str(iser) '\' num2str(exp_this) '\' animal '_' num2str(itr) '.eye'];
                    [data,timestamps,meta]=func_ReadFleaData(eye_file); % load eye file
                    clf(fig_handle);
                    for i = 1:length(iFr:min([iFr+8, idx_blink_offset]))
                        imagesc(data(:,:,iFr+i-1),'parent', subplot(3,3,i));
                        title(sprintf('trial: %d, frame: %d', itr, iFr+i-1));
                    end
                    
                    isBlink = menu('                Is this blink?            ', 'blink', 'not blink') == 1; % true if blink
                    if isBlink
                        blink_thistrial = [blink_thistrial; [ET_timestamp{itr}(iFr), ET_timestamp{itr}(idx_blink_offset)+median(diff(ET_timestamp{itr}))]];
                    end
                end
            end
            blink{itr} = blink_thistrial;
        end
        close(fig_handle); drawnow;
    end
    %% two photon
    % tp timestamp
    twoP_frame_onset = h5read(filepath, sprintf('/%d/twoP_frame_onset',exp_this));
    
    
    % load Protocol
    if any(cellfun(@(x) strcmp(expType(iexp),x), {'ori','texture'})) % if protocol exists
        do_make_folder = false; % no subfolder "ANALYZED" made
        [protocol, successflag] = ProtocolLoad(animal, str2double(iser), exp_this, do_make_folder);
    else
        protocol.seqnums = 1:ntri(iexp); % non-"ori" type experiments are treated as 1 stim x ntri(iexp) repeats.
    end
    
    % sort for each trial
    % twoP_F frame index for each trial onset
    exp_onset = sum(twoP_n_frame(1:iexp-1)); % idx for each exp onset (NOTE: 0-index)
    if any(cellfun(@(x) strcmp(expType(iexp),x), {'ori','texture'}))% for non-continuous recording by ND Sequence Acqusition
        trial_onset = [1; 1+ find(diff(twoP_frame_onset) > 1)]; % gap > 1 sec
    else
        trial_onset = [];
        for itr = 1:ntri(iexp) % for continuous recording, set 1st 2p frame after 1st ET frame as trial onset
            trial_onset(itr) = find(twoP_frame_onset > ET_timestamp{itr}(1), 1,'first');
        end
    end
    min_n_frame = min(diff(trial_onset));
    frm_interval = median(diff(twoP_frame_onset));
    
    % check frame for mosaic onset/offset
    mosaicTTL_per_twoPFrame = zeros(size(twoP_frame_onset));
    for iTTL = 1:size(mosaicTiming,1)
        mosaicOnset_idx = find(mosaicTiming(iTTL,1) - frm_interval < twoP_frame_onset, 1, 'first'); % frame during which mosaic onset happened
        mosaicOffset_idx = find(mosaicTiming(iTTL,2) - frm_interval < twoP_frame_onset, 1, 'first'); % frame during which mosaic offset happened
        mosaicTTL_per_twoPFrame(mosaicOnset_idx:mosaicOffset_idx) = 1;
    end
    
    % twoP data tranform
    [twoP_F_reshaped, twoP_Fneu_reshaped, twoP_spks_reshaped] = deal(NaN(min_n_frame, size(protocol.seqnums,1), size(protocol.seqnums,2), size(twoP_F,1))); % initialization
    [pupilArea_reshaped, wheel_reshaped, blink_reshaped, wheel_onset_reshaped, saccade_reshaped, mosaicTTL_reshaped, twoP_frame_onset_reshaped] = deal(NaN(min_n_frame, size(protocol.seqnums,1), size(protocol.seqnums,2)));
    for itr = 1: ntri(iexp) % loop for stim x rep combination
        disp(sprintf('processing %d th of %d trials...',itr, ntri(iexp)));
        strt_idx = trial_onset(itr); % onset index for trial (stim x rep)
        [istm, irep] = ind2sub(size(protocol.seqnums),find(protocol.seqnums(:)==itr));
        if strt_idx + min_n_frame -1 <= size(twoP_frame_onset,1)
            n_twoPframe_thistr = min_n_frame;
            twoP_frame_onset_thistr = twoP_frame_onset(strt_idx:strt_idx + min_n_frame -1);
        else
            % hack: some 2p data finished in the middle of last repeat
            %             twoP_frame_onset_thistr = [twoP_frame_onset(strt_idx:end);...
            %                 (1:strt_idx + min_n_frame - 1 - size(twoP_frame_onset,1))'*frm_interval + twoP_frame_onset(end)];
            n_twoPframe_thistr = size(twoP_frame_onset,1)-strt_idx+1;
            twoP_frame_onset_thistr = twoP_frame_onset(strt_idx:end);
        end
        twoP_F_reshaped(1:size(twoP_frame_onset_thistr,1),istm,irep,:) = permute(twoP_F(:,(exp_onset+strt_idx):(exp_onset+strt_idx + size(twoP_frame_onset_thistr,1) -1)),[2,3,4,1]);
        twoP_Fneu_reshaped(1:size(twoP_frame_onset_thistr,1),istm,irep,:) = permute(twoP_Fneu(:,(exp_onset+strt_idx):(exp_onset+strt_idx + size(twoP_frame_onset_thistr,1) -1)),[2,3,4,1]);
        twoP_spks_reshaped(1:size(twoP_frame_onset_thistr,1),istm,irep,:) = permute(twoP_spks(:,(exp_onset+strt_idx):(exp_onset+strt_idx + size(twoP_frame_onset_thistr,1) -1)),[2,3,4,1]);
        twoP_frame_onset_reshaped(1:size(twoP_frame_onset_thistr,1),istm,irep) = twoP_frame_onset_thistr;
        
        if do_ET
            pupilArea_interp = interp1(ET_timestamp{itr}(1:min([size(ET_timestamp{itr},1), size(area{itr},2)])), area{itr}(1:min([size(ET_timestamp{itr},1), size(area{itr},2)]))',... % HACK...
                twoP_frame_onset_thistr);
        end
        % wheel velocity, wheel onset and saccade
        [wheelvel_per_frame, wheelonset_per_frame, saccade_per_frame] = deal([]);
        for iFrm = 1:size(twoP_frame_onset_thistr,1) %twoP_frame_onset(strt_idx:strt_idx + min_n_frame -1)'
            frm_start = find(wheeltimestamp(:,1) > twoP_frame_onset_thistr(iFrm),1,'first');
            frm_end = find(wheeltimestamp(:,1) < twoP_frame_onset_thistr(iFrm)+frm_interval,1,'last');
            wheelvel_per_frame(iFrm) = mean(abs(wheelvel(frm_start:frm_end)));
            wheelonset_per_frame(iFrm) = any(all([signchT> twoP_frame_onset_thistr(iFrm),...
                signchT < twoP_frame_onset_thistr(iFrm)+frm_interval],2));
            if do_ET
                saccade_per_frame(iFrm) = any(all([ET_timestamp{itr}(sacT{itr})> twoP_frame_onset_thistr(iFrm),...
                    ET_timestamp{itr}(sacT{itr}) < twoP_frame_onset_thistr(iFrm)+frm_interval],2));
            end
        end
        
        wheel_reshaped(1:size(twoP_frame_onset_thistr,1),istm,irep) = wheelvel_per_frame;
        wheel_onset_reshaped(1:size(twoP_frame_onset_thistr,1),istm,irep) = wheelonset_per_frame;
        mosaicTTL_reshaped(1:size(twoP_frame_onset_thistr,1),istm,irep) = mosaicTTL_per_twoPFrame(strt_idx:strt_idx+n_twoPframe_thistr -1);
        if do_ET
            pupilArea_reshaped (1:size(twoP_frame_onset_thistr,1),istm,irep) = pupilArea_interp;
            saccade_reshaped(1:size(twoP_frame_onset_thistr,1),istm,irep) = saccade_per_frame;
            % blink
            for iEv = 1: size(blink{itr},1) % skipped if blink{itr} is empty
                frm_start = find(blink{itr}(iEv,1)-frm_interval < twoP_frame_onset_thistr,1,'first');
                frm_end = find(blink{itr}(iEv,2)+frm_interval > twoP_frame_onset_thistr,1,'last');
                % blink_reshaped: boolean, 1 from brinking start until end
                blink_reshaped(frm_start:frm_end,istm,irep) = 1;
            end
            % plot area for sanity check
            figure; imagesc(pupilArea_reshaped(:,:)');
            title(sprintf('pupil area in exp:%d', exp_this));
            xlabel('frame');
            ylabel('stim x repeat');
        end
        
    end
    
    
    %% append to h5 file
    info = h5info(filepath, sprintf('/%d',exp_this));
    dataset_names = {info.Datasets.Name};
    % twoP_F
    if ~any(cellfun(@(x) strcmp(x, 'twoP_F'), dataset_names))
        h5create(filepath, sprintf('/%d/twoP_F',exp_this), size(twoP_F_reshaped));
    end
    h5write(filepath, sprintf('/%d/twoP_F',exp_this), twoP_F_reshaped);
    % twoP_Fneu
    if ~any(cellfun(@(x) strcmp(x, 'twoP_Fneu'), dataset_names))
        h5create(filepath, sprintf('/%d/twoP_Fneu',exp_this), size(twoP_Fneu_reshaped));
    end
    h5write(filepath, sprintf('/%d/twoP_Fneu',exp_this), twoP_Fneu_reshaped);
    % twoP_spks
    if ~any(cellfun(@(x) strcmp(x, 'twoP_spks'), dataset_names))
        h5create(filepath, sprintf('/%d/twoP_spks',exp_this), size(twoP_spks_reshaped));
    end
    h5write(filepath, sprintf('/%d/twoP_spks',exp_this), twoP_spks_reshaped);
    % twoP_frame_onset
    if ~any(cellfun(@(x) strcmp(x, 'twoP_frame_onset_reshaped'), dataset_names))
        h5create(filepath, sprintf('/%d/twoP_frame_onset_reshaped',exp_this), size(twoP_frame_onset_reshaped));
    end
    h5write(filepath, sprintf('/%d/twoP_frame_onset_reshaped',exp_this), twoP_frame_onset_reshaped);
    % wheel velocity
    if ~any(cellfun(@(x) strcmp(x, 'wheel_velocity'), dataset_names))
        h5create(filepath, sprintf('/%d/wheel_velocity',exp_this), size(wheel_reshaped));
    end
    h5write(filepath, sprintf('/%d/wheel_velocity',exp_this), wheel_reshaped);
    % wheel onset
    if ~any(cellfun(@(x) strcmp(x, 'wheel_onset'), dataset_names))
        h5create(filepath, sprintf('/%d/wheel_onset',exp_this), size(wheel_onset_reshaped));
    end
    h5write(filepath, sprintf('/%d/wheel_onset',exp_this), wheel_onset_reshaped);
    % mosaicTTL
    if ~any(cellfun(@(x) strcmp(x, 'mosaicTTL'), dataset_names))
        h5create(filepath, sprintf('/%d/mosaicTTL',exp_this), size(mosaicTTL_reshaped));
    end
    h5write(filepath, sprintf('/%d/mosaicTTL',exp_this), mosaicTTL_reshaped);
    % photodiode_onset
    if ~any(cellfun(@(x) strcmp(x, 'photodiode_onset'), dataset_names))
        h5create(filepath, sprintf('/%d/photodiode_onset',exp_this), size(photodiode_onset));
    end
    h5write(filepath, sprintf('/%d/photodiode_onset',exp_this), photodiode_onset);
    % photodiode_offset
    if ~any(cellfun(@(x) strcmp(x, 'photodiode_offset'), dataset_names))
        h5create(filepath, sprintf('/%d/photodiode_offset',exp_this), size(photodiode_offset));
    end
    h5write(filepath, sprintf('/%d/photodiode_offset',exp_this), photodiode_offset);
    % seqnums
    if ~any(cellfun(@(x) strcmp(x, 'seqnums'), dataset_names))
        h5create(filepath, sprintf('/%d/seqnums',exp_this), size(protocol.seqnums));
    end
    h5write(filepath, sprintf('/%d/seqnums',exp_this), protocol.seqnums);
    % EYETRACKING
    if do_ET
        if any(cellfun(@(x) strcmp(x, 'pupilArea'), dataset_names)) || any(cellfun(@(x) strcmp(x, 'saccade'), dataset_names))
            str = input('pupil area area or saccade (or both) already saved! overwrite? [y/n]','s');
            if strcmpi(str,'y') % overwrite
                h5write(filepath, sprintf('/%d/pupilArea',exp_this), pupilArea_reshaped);
                h5write(filepath, sprintf('/%d/saccade',exp_this), saccade_reshaped);
                h5write(filepath, sprintf('/%d/blink',exp_this), blink_reshaped);
            end
        else
            
            % pupilArea
            if ~any(cellfun(@(x) strcmp(x, 'pupilArea'), dataset_names))
                h5create(filepath, sprintf('/%d/pupilArea',exp_this), size(pupilArea_reshaped));
            end
            h5write(filepath, sprintf('/%d/pupilArea',exp_this), pupilArea_reshaped);
            % saccade
            if ~any(cellfun(@(x) strcmp(x, 'saccade'), dataset_names))
                h5create(filepath, sprintf('/%d/saccade',exp_this), size(saccade_reshaped));
            end
            h5write(filepath, sprintf('/%d/saccade',exp_this), saccade_reshaped);
            % blink
            if ~any(cellfun(@(x) strcmp(x, 'blink'), dataset_names))
                h5create(filepath, sprintf('/%d/blink',exp_this), size(blink_reshaped));
            end
            h5write(filepath, sprintf('/%d/blink',exp_this), blink_reshaped);
       end
    end
    
   
    
    
end

