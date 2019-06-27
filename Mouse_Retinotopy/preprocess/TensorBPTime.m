function tensor = TensorBPTime(  LoCutFreq, HiCutFreq, FrameRate, tensor, graphicsflag )
% TensorBPTime bandpass-filters in time using FFT method
%
% TensorBPTime( LoCutFreq, HiCutFreq ) takes as input VDAQ and
% lets you specify the low and high cutoffs LoCutFreq, HiCutFreq in
% Hz.
%
% Notice that the filtering is done in the frequency domain, with the
% sharpest possible cutoffs (equivalent to convolving with sinc functions).
% 
% It figures out what they correspond in terms of frames from the field
% VDAQ.FrameRate
%
% Set LoCutFreq to 0   if you don't want to high-pass filter
% Set HiCutFreq to Inf if you don't want to low-pass filter
%
% tensor = TensorBPTime( LoCutFreq, HiCutFreq, FrameRate, tensor)
% lets you specify the units FrameRate (DEFAULT: [], which means
% VDAQ.FrameRate) and a tensor that can be in the following formats:
% tensor{nstim}(nx ny frames)
% tensor(nx ny nframes)
%
% tensor = TensorBPTime( ..., 'showgraphics') shows the power spectra
% before and after filtering
%
% See also: TensorBandPassTime

% THERE IS CURRENTLY NO WINDOWING

% 2006-03 Matteo Carandini, based on TensorBandPassTime

global VDAQ

if nargin < 5, graphicsflag = [];   end
if nargin < 4, tensor = [];   end
if nargin < 3, FrameRate = []; end
if nargin < 2, error('Need to specify at least two inputs'); end

if LoCutFreq>HiCutFreq, warning('You got the low cutoff and the high cutoff switched around'); end

if isempty(tensor)
    goVDAQ = 1;
else
    goVDAQ = 0;
end

if isempty(FrameRate)
    FrameRate = VDAQ.FrameRate;
end

if HiCutFreq<=LoCutFreq
    error('Highpass and lowpass specs do not make sense - not doing highpass');
end

%------------------------------------

if goVDAQ
    nstim = length(VDAQ.tensor);
    [nx, ny, nt] = size(VDAQ.tensor{1});
else
    if iscell(tensor)
        % format is tensor{nstim}(nx ny frames)
        nstim = numel(tensor);
        [nx, ny, nt] = size(tensor{1});
    elseif ndims(tensor)==3
        [nx, ny, nt] = size(tensor);
    end
end

ff = freq(nt, nt/FrameRate);

% Code that might be useful if one day we want to use Gaussian filters:
% ftFilter = zeros(1,1,VDAQ.nframes);
% ftFilter(:) = 1 - gaussian( sigmaf, ff );
% ftFilter = repmat(ftFilter,[VDAQ.nx, VDAQ.ny, 1]);
% VDAQ.tensor{istim} = real(ifft( ftFilter.*ftTensor, [], 3 ));

%------------------- do the work ----------------------------------

%fprintf(1, 'Filtering in time');

if goVDAQ

    for istim = 1:nstim
        %fprintf(1,'.');
        ftTensor = fft( VDAQ.tensor{istim}, [], 3 );
        if LoCutFreq~=0, ftTensor(:,:,abs(ff) <= LoCutFreq) = 0; end
        ftTensor(:,:,abs(ff) >= HiCutFreq) = 0;
        VDAQ.tensor{istim} = real(ifft( ftTensor, [], 3 ));
    end

else
    if iscell(tensor)
        % format is tensor{nstim}(nx ny frames)
        for istim = 1:nstim
            %fprintf(1,'.');
            ftTensor = fft( tensor{istim}, [], 3 );
            if LoCutFreq~=0, ftTensor(:,:,abs(ff) <= LoCutFreq) = 0; end
            ftTensor(:,:,abs(ff) >= HiCutFreq) = 0;
            tensor{istim} = real(ifft( ftTensor, [], 3 ));
        end

    elseif ndims(tensor)==3

        %fprintf(1,'...');
        ftTensor = fft( tensor, [], 3 );
        if LoCutFreq~=0, ftTensor(:,:,abs(ff) <= LoCutFreq) = 0; end
        ftTensor(:,:,abs(ff) >= HiCutFreq) = 0;
        tensor = real(ifft( ftTensor, [], 3 ));

    else

        fprintf(1,'Do not quite understand the format of the input\n')
        return
    end
end

%fprintf(1, 'done\n');

return



