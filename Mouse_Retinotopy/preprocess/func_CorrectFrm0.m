function func_CorrectFrm0(f0List)
%Blank and average blank corrections
% input 1 for average blank and 0 for blank correction
global VDAQ
nt = size(VDAQ.tensor{1},3);

for istim = 1:VDAQ.nstim
    Tmp = mean(VDAQ.tensor{istim}(:,:,f0List),3);
    VDAQ.Frame0List = Tmp;
    TheCorrection = repmat(Tmp,[1 1 nt]);
    VDAQ.tensor{istim} = VDAQ.tensor{istim} - TheCorrection;
    fprintf('Frm0 Correction on stim #%d ...\n', istim);
end
clear Tmp;


%% Log
%{
- Mohammad 2015; wrote the file
%}