function [TimeSer, Frq, Amp, ROI] = func_plot_VDAQFourietAmpSpect(VDAQ,ROI)
% Compute fourier amplitue spectrum in ROIs using PowerTimeTrace function
clc
display('Choose ROI for each stim, and press enter ...')

if isempty(ROI)
    FgH = figure; clf;
    clf; drawnow
    istim =1;
    Tmp{istim} = squeeze(std(VDAQ.tensor{istim},[],3));
    mymin = prctile(Tmp{istim}(:), 2.5);
    mymax = prctile(Tmp{istim}(:), 97.5);
    h_im = imagesc(fliplr(rot90( Tmp{istim},2)), [mymin mymax]); %,[min(min(Tmp{istim})) max(max(Tmp{istim}))]
    axis image; axis xy; colormap bone;
    e = imrect(gca,[10 10 100 100]);
    pause;
    ROI = createMask(e,h_im);
    close(FgH)
end
clear Tmp e R C h_im mymax mymin;

% estimate amplitude
TimeSer = []; Frq = []; Amp = []; Tmp = {};
for istim = 1:size(VDAQ.tensor, 2)
    [R, C] = find(ROI);
    TimeSer(istim,:)     = squeeze(mean(mean(VDAQ.tensor{istim}(R(1):R(end),C(1):C(end),:))));
    [Frq, Amp(istim,:)] = func_PowerTimeTrace(TimeSer(istim,:), VDAQ.durs); 
end

% Plot amplitude spectrum
figure('color','w')
subplot(2,1,1)
for istim=1:size(VDAQ.tensor, 2)
    col = {'k'}; if istim == size(VDAQ.tensor, 2), col = {'b'}; end
    plot(VDAQ.tt, TimeSer(istim,:), col{1}); hold on    
end
xlim([0 VDAQ.durs(1)]); xlabel('time (s)'); ylabel('resp')
set(gca,'tickdir','out'); set(gcf,'color','w'); box off
pbaspect([3,1,1])

subplot(2,1,2)
idxFrq = find(Frq > 0.02);
for istim=1:size(VDAQ.tensor, 2)
    col = {'k'}; if istim == size(VDAQ.tensor, 2), col = {'b'}; end
    plot(Frq(idxFrq), Amp(istim,idxFrq), col{1}); hold on    
end
xlim([0 2]); xlabel('Freq (Hz)'); ylabel('Amplitude')
set(gca,'tickdir','out'); set(gcf,'color','w'); box off
pbaspect([3,1,1])


function [Freq, Amp] = func_PowerTimeTrace(tTrace, dur)
% POWERTIMETRACE computes the power of a generic time trace
% [ff Py] = PowerTimeTrace(tTrace, dur)
% tTrace is a [1,number of time points] vector, 'dur' is a scalar for duration in seconds
% 'ff' is the frequency axis and Py the amplitude spectrum

% get time pars
nt = length(tTrace);

%  Compute amplitude with fft
Fs = nt/dur;
x = 1+tTrace;

% Alternatives for actual power - don't have the toolboxes
% [Py(istim,:),f(istim,:)] = pwelch(avgs(istim,:),nt,[],[],VDAQ.FrameRate);
% [Py(istim,:),f(istim,:)] = periodogram(avgs(istim,:),nt,[],VDAQ.FrameRate);

xdft = fft(x); xdft = xdft(1:round(nt/2)+1);
psdx = (1/(Fs*nt)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
Amp = psdx;
Freq = 0:Fs/length(x):Fs/2;
return


%% Log
%{
- Mohammad 2015; wrote the file
%}



