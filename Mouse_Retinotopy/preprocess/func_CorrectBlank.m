function func_CorrectBlank(blanklist, avgB)
%Blank and average blank corrections
% avgB: input 1 for average blank and 0 for blank correction
% if avgB=1, the blank tensors are first averaged over time, and then across trials
% if avgB=0, the blank tensors are avergaed across trials.

global VDAQ
% avgB = 0;

% -------- BLANK CORRECTION ----------
for istim = 1:VDAQ.nstim
    nt = size(VDAQ.tensor{istim},3);
    Correction = zeros(size(VDAQ.tensor{istim}));
    
    for iblank = 1:length(blanklist)
        iblnk = blanklist(iblank);
        if avgB
            Tmp = repmat(squeeze(mean(VDAQ.tensor{iblnk},3)),[1 1 nt]);
            fprintf('Avg Blank Correction ...\n');
        else
            Tmp = VDAQ.tensor{iblnk};
            fprintf('Blank Correction ...\n');
        end
        Correction = Correction + Tmp/length(blanklist);
    end
    clear Tmp;
    
    % Do the actual correction
    %for istim = 1:size(VDAQ.tensor,2)
    VDAQ.tensor{istim} = VDAQ.tensor{istim} - Correction;
    clear TheCorrection1
end

%% Log
%{
- Mohammad 2015; wrote the file
%}