close all
clear all

%% Receiver Parameters

effectiveBW = 900e6; %Nyquist Region for Fs = p*channelFs
p = 18; %Undersampling Factor
nfft = 128;

effectiveFs = effectiveBW*2;
channelFs = effectiveFs/p;

wrapTime = 1/effectiveFs;

timeDelays = [0,.33,1]'*wrapTime;

%% Signal Model
t = (0:nfft*p-1)*1/effectiveFs; %Time vector for non-aliased signal
f = [137.1e6,239e6,871.375e6]; %signal frequencies
A = [1;1;1];

% Figure out time vector for different channels after application of time
% delay

t = repmat(t,3,1);
t = t+repmat(timeDelays,1,size(t,2));

% Generate signal for different channels (i.e. at different delays)
mask=[zeros(3,30*p),ones(3,p*(nfft-30))];
mask = repmat(blackman(p*nfft)',3,1);
mask = ones(3,p*nfft);

s = zeros(3,p*nfft);
for i = 1:3
    temp = repmat(A,1,size(t,2)).*sin(2*pi*f'*t(i,:)).*mask;
    s(i,:) = sum(temp,1);
end


%Downsample
sd = s(:,1:p:end);
sF = zeros(3,nfft);
for i = 1:3
    sF(i,:) = fft(sd(i,:),nfft);
end

%% Estimate Frequency

[~, buckets] = findpeaks(abs(sF(1,1:nfft/2))); % Find peaks in buckets
% buckets = find(abs(sF(1,1:nfft/2)).^2 > 0);

b1 = sF(1,:); %Non-delayed spectrum
b2 = sF(3,:); %Delayed spectrum

values = b1(buckets);
%For each occupied bucket, estimate the true frequency by comparing phase
%difference between the delayed spectrum and non-delayed spectrum
freqEst = zeros(1,length(buckets));
for i = 1:length(buckets)
    ind = buckets(i);
    phaseEst = abs(angle(b2(ind)/b1(ind)));
    freqEst(i) = phaseEst/(2*pi*timeDelays(3));
end

%% Construct xhat

freq = 0:effectiveFs/(p*nfft):(effectiveFs-effectiveFs/(p*nfft));
xf = zeros(size(freq));
for i = 1:length(buckets)
   [v,k] = min(abs(freqEst(i)-freq));
   xf(k) = values(i);
   freq(k) = inf;
end

sEst = ifft(xf,nfft*p);

subplot(311)
plot(2*p*imag(sEst))
subplot(312)
plot(s(1,:))
subplot(313)
plot(abs(imag(sEst)-s(1,:)))

figure()
freq = 0:effectiveFs/(p*nfft):(effectiveFs-effectiveFs/(p*nfft));
plot(freq/1e6,abs(fft(sEst,p*nfft)).^2)

