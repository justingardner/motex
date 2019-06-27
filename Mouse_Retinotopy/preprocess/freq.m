function f = freq(nsamples, duration)
%FREQ frequencies corresponding to the elements of an fft
%
% 	freq(nsamples, duration) returns a vector of frequencies
%	[ 0 1/duration 2/duration ... -1/duration ]
%
%
% 1995 Matteo Carandini
% 2004 VM separated even and odd case
%
% part of the Matteobox toolbox

if mod(nsamples,2)
   %odd (the original freq by MC did always the odd case)
   f = [0:floor(nsamples/2) , -floor((nsamples-1)/2):1:-1]/duration;
else
   %even
   f = [0:(floor(nsamples/2)-1) , -floor((nsamples-1)/2)-1:1:-1]/duration;
end



return


tt = 1:11;
rr = ones(size(tt));

RR = abs(fftshift(fft(rr)))


fftshift(myfreq(10,10))

return

%% Example of usage

% Example:

w = 4; % Hz frequency of the stimulus

fs = 1000;
dur = 1; 

tt = (1:1000)/1000;

signal = 5*sin( w*2*pi*tt ) + normrnd(0,4,1, 1000);

figure; plot(tt, signal);

ftsignal = fft(signal);

ff = freq(dur*fs,dur);

figure; stem( ff, abs(ftsignal) );

signalamplitude = abs(ftsignal(ff==w)) / (1000/2);
