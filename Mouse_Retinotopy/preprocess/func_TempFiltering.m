function tensor = func_TempFiltering(LoCutFreq, HiCutFreq, tensor, FrameRate)
%% Temporal Filtering in Hz (in Freq domain)
global VDAQ

disp(['Cutoff Fqs: [' num2str(LoCutFreq) ', ' num2str(HiCutFreq) '], Temporal filtering ...']);
if nargin < 3
    for istim = 1:VDAQ.nstim
        TmpFilt_Tensor = TensorBPTime(  LoCutFreq, HiCutFreq, VDAQ.FrameRate, VDAQ.tensor{istim}, [] );
        VDAQ.tensor{istim} = TmpFilt_Tensor;
    end
else
    if nargin < 4,      FrameRate = 30;    end
    if ~iscell(tensor), tensor = {tensor}; end
    for istim = 1:numel(tensor)
        tensor{istim} = TensorBPTime(  LoCutFreq, HiCutFreq, FrameRate, tensor{istim}, [] );
    end
end


%% Log
%{
- Mohammad 2015; wrote the file
%}