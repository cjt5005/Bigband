close all
clear all






%% Receiver Parameters

effectiveBW = 900e6; %Nyquist Region for Fs = p*channelFs
p = 4; %Undersampling Factor
nfft = 128;

effectiveFs = effectiveBW*2;
channelFs = effectiveFs/p;

wrapTime = 1/effectiveBW;

timeDelays = [0,.33,1]'*wrapTime;


%% Experiment with Phase Shifts

t1 = (0:nfft*p-1)*1/effectiveFs;
f = 411e6;
timeDelay = wrapTime;

t2 = t1+timeDelay;

s1 = sin(2*pi*f*t1);
s2 = sin(2*pi*f*t2);

plot(s1)
hold on
plot(s2)

s1d = s1(1:p:end);
s2d = s2(1:p:end);
figure()
plot(s1d)
hold on
plot(s2d)

%% Signal Model
t = (0:nfft*p-1)*1/effectiveFs; %Time vector for non-aliased signal
f = [122e6,287e6,411e6]; %signal frequencies
A = [1;1;1];

% Figure out time vector for different channels after application of time
% delay

t = repmat(t,3,1);
t = t+repmat(timeDelays,1,size(t,2));

% Generate signal for different channels

for i = 1:3
    temp = repmat(A,1,size(t,2)).*sin(2*pi*f'*t(i,:));
    s(i,:) = sum(temp,1);
end
plot(s')

%Undersample

s = s(:,1:p:end);
figure()
plot(s')
for i = 1:3
    sF(i,:) = fft(s(i,:),nfft);
end

%% Estimate Frequency

bucketThreshold = 2000;
buckets = find(abs(sF(1,1:nfft/2)).^2>1000);

b1 = sF(1,:); %Non-delayed spectrum
b2 = sF(3,:); %Delayed spectrum

for i = 1:length(buckets)
    deltaPhase(i) = wrapTo2Pi((angle(b2(buckets(i))/b1(buckets(i)))));
    freqEst(i) = deltaPhase(i)/(2*pi*timeDelays(3));
end

freqEst/1e6
freq = 0:channelFs/nfft:(channelFs-channelFs/nfft);
figure()
for i = 1:3
%     subplot(121)
    magSpect = abs(sF(i,:)).^2;
    hold on
    plot(freq(1:nfft/2)/1e6,magSpect(1:nfft/2))
%     subplot(122)
%     hold on
%     angSpect = unwrap(angle(sF(i,:)));
%     plot(freq(1:nfft/2)/1e6,angSpect(1:nfft/2))
end