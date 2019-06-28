% example code to analyze texture data (in collaboration with Dr.
% Gardner@Standford)
% author: Ryo Aoki
% last updated@190627

%%

animal = 'M190625_18273';
iser = '1';
h5path = fullfile('//NCB-LABSERVER6/data/MOUSE/IMAGING/GCAMP', animal, 'Preprocessed', sprintf('%s_%s.h5', animal, iser));

frameRate = 15.24; % Hz

%% load
% expSet
%expSet = h5read(h5path,sprintf('/expSet'));
exp_this = h5read(h5path,sprintf('/expSet'));
twoP_isCell = h5read(h5path, sprintf('/twoP_iscell'));
idx_validCell = find(logical(twoP_isCell(:,1)));

twoP_F = h5read(h5path, sprintf('/%d/twoP_F',exp_this)); % timebin x stim x rep x cell
twoP_Fneu = h5read(h5path, sprintf('/%d/twoP_Fneu',exp_this));
twoP_frame_onset = h5read(h5path, sprintf('/%d/twoP_frame_onset_reshaped',exp_this));
photodiode_onset= h5read(h5path, sprintf('/%d/photodiode_onset',exp_this));
photodiode_offset= h5read(h5path, sprintf('/%d/photodiode_offset',exp_this));
seqnums= h5read(h5path, sprintf('/%d/seqnums',exp_this));


%% compute deltaF/F
baseline_window = [0, 2.3]; % 0 - 2.3 sec is used for baseline
baseline_idx = floor(baseline_window(1)*frameRate)+1:floor(baseline_window(2)*frameRate);
twoP_F_corrected = twoP_F - twoP_Fneu * 0.7; % neuropil-correction
twoP_dF_F = bsxfun(@(x,y) (x-y)/y, twoP_F_corrected, mean(mean(mean(twoP_F_corrected(baseline_idx,:,:,:),1),2),3));

% sanity check
twoP_F_reshaped = permute(twoP_dF_F,[4, 1, 2, 3]);
figure; imagesc(twoP_F_reshaped(idx_validCell,:)); 
set(gca,'clim',[-0.2 1]); 
xlabel('timebin x stim x rep'); ylabel('cell');

%% some example cells
for iCell = 1:10
    figure;
    dF_F_thiscell = twoP_dF_F(:,:,:,idx_validCell(iCell)); % timebin x stim x rep

    subplot(221);
    errorbar((1:size(twoP_dF_F,1))/frameRate, mean(dF_F_thiscell(:,:),2), std(dF_F_thiscell(:,:),[],2)./sqrt(size(dF_F_thiscell(:,:),2)));
    title(sprintf('cell: %d (%d), mean response', idx_validCell(iCell), iCell));
    axis tight;
    xlabel('time [sec]');
    ylabel('dF/F');
    
    subplot(223);
    imagesc(dF_F_thiscell(:,:,1)');
    set(gca,'clim',[-0.2 3],'xtick',([0:2.5:9.9].*frameRate),'xticklabel',[0:2.5:9.9],'tickdir','out');
    colorbar;
    xlabel('time [s]');
    ylabel('stim');
    title('repeat 1');

    subplot(224);
    imagesc(dF_F_thiscell(:,:,2)');
    set(gca,'clim',[-0.2 3],'xtick',([0:2.5:9.9].*frameRate),'xticklabel',[0:2.5:9.9],'tickdir','out');
    colorbar;
    xlabel('time [s]');
    ylabel('stim');
    title('repeat 2');  
    
end

%% timing sanity check
i = 0; figure; hold on;
for iStim = 1: size(twoP_F,2)
    for iRep = 1:size(twoP_F,3)
        i = i+1;
        plot(twoP_frame_onset(:,iStim,iRep)-twoP_frame_onset(1,iStim,iRep), repmat(i,[1,size(twoP_F,1)]),'.k');
        plot(photodiode_onset(seqnums(iStim,iRep))-twoP_frame_onset(1,iStim,iRep),i,'ro');
        plot(photodiode_offset(seqnums(iStim,iRep))-twoP_frame_onset(1,iStim,iRep),i,'bo');
    end
end
xlabel('time [s]');
ylabel('stim x rep');









