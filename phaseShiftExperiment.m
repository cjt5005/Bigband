clear all
close all
clc

%% Setup Parameters

effectiveBW = 900e6; %Nyquist Region for Fs = p*channelFs
p = 4; %Undersampling Factor
nfft = 128;

effectiveFs = effectiveBW*2;
channelFs = effectiveFs/p;

wrapTime = 1/effectiveBW;


%% Experiment with Phase Shifts

t1 = (0:nfft*p-1)*1/effectiveFs;
f = 411e6;
timeDelay = wrapTime;

t2 = t1+timeDelay;

s1 = exp(2i*pi*f*t1);
s2 = exp(2i*pi*f*t2);

plot(s1)
hold on
plot(s2)

s1d = s1(1:p:end);
s2d = s2(1:p:end);
figure()
plot(imag(s1d))
hold on
plot(imag(s2d),'r')