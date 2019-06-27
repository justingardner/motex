function [AbsMaps, AngleMaps, FourierMaps] = TensorFrequency(StimList, freqlist, CoherenceFlag)
% TensorFrequency performs power and phase analyses of VSD data 
% 
% [AbsMaps, AngleMaps] = TensorFrequency(StimList, freqlist) operates on
% stimuli in StimList of VDAQ.tensor and returns amplitude AbsMaps and phase
% AngleMaps of the responses at frequencies freqlist. The phase AngleMaps is
% the angle of the Fourier component (from -pi to pi)
%
% [AbsMaps, AngleMaps, FourierMaps] = TensorFrequency(StimList, freqlist)
% also returns the complex maps
%
% [AbsMaps, AngleMaps, FourierMaps] = TensorFrequency(StimList, freqlist,1) 
% returns the signal coherence in AbsMaps, which is AbsMaps./noisePow,
% where noisePow is the sum of all the frequency components.
%
% Use 'PlotAmpPhaseMaps' to plot the results.
%
% 2004-04 AB MC and A.Wade
% 2004-05 MC + AB removed calls to fft
% 2004-06 AB introduced VDAQ structure and pi/2 correction
% 2005-04 MC removed the pi/2 "correction" which was a problem
% 2005-12 MC made the computation of CoherenceMaps conditional on nargout
% 2005-12 MC handled special case in which frequency is 0
% 2006-09 MC from TensorCoherence

if nargin == 2
    CoherenceFlag = 0;
end

global VDAQ
nstim  = VDAQ.nstim-1;
% nstim = length(VDAQ.tensor);

if isempty(StimList), StimList = 1:nstim; end

nframes = size(VDAQ.tensor{1},3);
durs    = VDAQ.durs; % might be safer to use expt.stimdurs
dur = durs(1);

AbsMaps     = {};
AngleMaps   = {};
FourierMaps = {};

[nx, ny, ~] = size(VDAQ.tensor{1});
DiffFreqs = unique(freqlist);

for iDiffFreq = 1:length(DiffFreqs)
    ThisFreq = DiffFreqs(iDiffFreq);
    ComplexMaps = TensorFreq( {VDAQ.tensor{StimList}}, ThisFreq, dur );
    for iStimList = 1:length(StimList)
        if ThisFreq==0
            AbsMaps{iStimList,iDiffFreq}   = ComplexMaps{iStimList};
            AngleMaps{iStimList,iDiffFreq} = zeros(nx,ny);
            FourierMaps{iStimList,iDiffFreq} = ComplexMaps{iStimList};
        else
            AbsMaps{iStimList,iDiffFreq}   = abs(ComplexMaps{iStimList});
            AngleMaps{iStimList,iDiffFreq} = angle(ComplexMaps{iStimList}); % AB had added + pi/2;
            FourierMaps{iStimList,iDiffFreq} = ComplexMaps{iStimList};
        end
    end
end

if CoherenceFlag
    StdMap = {};
    for iStimList = 1:length(StimList)
        StdMap{iStimList} = std( VDAQ.tensor{iStimList}, [], 3 );
        for iDiffFreq = 1:length(DiffFreqs)
            AbsMaps{iStimList,iDiffFreq}  = AbsMaps{iStimList,iDiffFreq} ./ StdMap{iStimList};
        end
    end
end

