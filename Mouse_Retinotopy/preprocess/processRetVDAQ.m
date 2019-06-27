function [p,stimfreqs] = processRetVDAQ(p,params)
% this function performs all necessary processing steps to the saved VDAQ
% structures. those files are only df/f normalized.
% steps:
%   1. trimming invalid frames from the end of files
%   2. spatial filtering
%   3. temporal filtering
%   4. blank correction
%   5. frame0 correction


if isempty(params)
    params = struct('LoCutFreq',0.05,'HiCutFreq',2,'spatfiltwid',3);
end

global VDAQ
% trim invalid frames from tensor
    figure; meansig = cell(size(VDAQ.tensor));
    for istim = 1:p.nstim
        for iten = 1:size(VDAQ.tensor,1)
            col = {'ko-'}; if istim == p.blankstims, col = {'ro-'}; end
            meansig{iten,istim} = squeeze(mean(mean(VDAQ.tensor{iten,istim})));
            plot(meansig{iten,istim}, col{1}); hold on
        end
    end
    title('Click to choose an appropriate Threshold (y) for trimming')
    [~, tshld2] = ginput(1); close;
    keepidx = cellfun(@gt,meansig,repmat({tshld2},size(meansig)),'uni',0);
    [~,tmpidx] = min(reshape(cellfun(@sum,keepidx),[],1));
    keepidx = repmat(keepidx(tmpidx),size(keepidx));
    minsz = cellfun(@size,VDAQ.tensor,repmat({3},size(VDAQ.tensor)));
    S = struct('type','()','subs',cellfun(@transpose,cellfun(@squeeze,mat2cell(cat(3,repmat({':'},size(minsz)),repmat({':'},size(minsz)),keepidx),ones(size(keepidx,1),1),ones(size(keepidx,2),1),3),'uni',0),'uni',0));
    VDAQ.tensor = cellfun(@subsref,VDAQ.tensor,num2cell(S),'uni',0);
    [nr,nc,nt]   = size(VDAQ.tensor{1});
    VDAQ.durs    = nt/VDAQ.FrameRate;
    VDAQ.nframes = nt;
    VDAQ.ny      = nr;
    VDAQ.nx      = nc;
    VDAQ.tt      = linspace(0, VDAQ.durs, nt);
    p.pars(1,:)  = ones(size(p.pars(1,:))).*(VDAQ.durs*10);
% prepare p-file for filtering procedures
if isempty(p.blankstims), p.blankstims = p.nstim; end
p.pfilefreqs = p.pars(4,:)/100;
stimfreqs = p.pfilefreqs;
stimfreqs(p.blankstims) = [];

VDAQ_full = VDAQ;
for k = 1:size(VDAQ_full.tensor,1)
    VDAQ.tensor = VDAQ_full.tensor(k,:);
    % Temporal and Spatial Filtering
    func_TempFiltering(params.LoCutFreq, params.HiCutFreq);
    func_SpatFiltering(params.spatfiltwid);
    % Blank and Frame0 corrections
    %blanklist = p.blankstims;
    %func_CorrectBlank(blanklist, 1);
    func_CorrectFrm0(params.frame0list);

    VDAQ_full.tensor(k,:) = VDAQ.tensor;
end
VDAQ = VDAQ_full;



%% Log
%{
- Mohammad 2015; wrote the file
- mjm organized 2016
%}


