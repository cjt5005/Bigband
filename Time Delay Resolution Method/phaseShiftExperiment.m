clear all
close all
clc

%% Setup Parameters

effectiveBW = 900e6; %Nyquist Region for Fs = p*channelFs
p = 20; %Undersampling Factor
nfft = 128;

effectiveFs = effectiveBW*2;
channelFs = effectiveFs/p;

wrapTime = 1/effectiveFs;


%% Time Domain

t1 = (0:nfft*p-1)*1/effectiveFs;
f = 700e6;
timeDelay = wrapTime;

t2 = t1+timeDelay;

s1 = exp(2i*pi*f*t1);
s2 = exp(2i*pi*f*t2);

% plot(s1)
% hold on
% plot(s2)

s1d = s1(1:p:end);
s2d = s2(1:p:end);
figure()
plot(imag(s1d))
hold on
plot(imag(s2d),'r')

phaseEstimate = angle(s2d/s1d);
f = phaseEstimate/(2*pi*timeDelay);

%% Frequency Domain

s1d2 = imag(s1d);
s2d2 = imag(s2d);
w = 0:channelFs/nfft:(channelFs-channelFs/nfft);
s1f = fft(s1d2,nfft);
s2f = fft(s2d2,nfft);

figure()
plot(w(1:nfft/2),abs(s1f(1:nfft/2)))

% Find peak (for f = 411e6, peak is in bin 30)

peak = find(abs(s1f(1:nfft/2))>.9*max(abs(s1f)));

bhat1 = s1f(peak);
bhat2 = s2f(peak);

phaseEst = abs(angle(bhat2/bhat1))

freqEst = phaseEst/(2*pi*timeDelay)