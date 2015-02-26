close all
clear all

%% Receiver Parameters

effectiveBW = 900e6; %Nyquist Region for Fs = p*channelFs
p = 9; %Undersampling Factor
nfft = 128;

effectiveFs = effectiveBW*2;
channelFs = effectiveFs/p;

wrapTime = 1/effectiveFs;

timeDelays = [0;.33;1]*wrapTime;

%% Signal Model

f = [85e6;140e6;450e6]; %signal frequencies
A = [1;1;1];
phi = rand(3,1)*2*pi;

% f=340e6;
% A=1;
% phi = rand()*2*pi;
%Generate signal at multiple delays
[s,t] = sigGen(A,f,phi,effectiveFs,timeDelays,p*nfft);

%Resample
if p<1
    sd = upsample(s,1/p);
elseif p>1
    sd = downsample(s,p);
end

%% Frequency Bucketization

% Bucketize the signal (i.e. Compute the periodogram of the resampled
% signal)
sF = fft(sd,nfft);
freq = 0:channelFs/(nfft):(channelFs-channelFs/(nfft));
pgram = abs(sF(1:nfft/2,1)).^2;
plot(freq(1:nfft/2)/1e6,pgram)
% Find the occupied buckets
[~, buckets] = findpeaks(abs(sF(1:nfft/2,1)));


%% Collision Detection

b1 = sF(buckets,1);
b2 = sF(buckets,2);
b3 = sF(buckets,3);

magDiff = abs(b3)-abs(b1);




% Save the value of the buckets
values = sF(buckets,1);

%For each occupied bucket, estimate the true frequency by comparing phase
%difference between the delayed spectrum and non-delayed spectrum
phaseEst = abs(angle(sF(buckets,3)./sF(buckets,1)));
freqEst = phaseEst/(2*pi*timeDelays(3));

%% Find Collisions

%% Construct xhat

freq = 0:effectiveFs/(p*nfft):(effectiveFs-effectiveFs/(p*nfft));
xf = zeros(p*nfft,1);

for i = 1:length(buckets)
   [v,k] = min(abs(freqEst(i)-freq));
   xf(k) = values(i);
   freq(k) = inf;
end

sEst = ifft(xf,nfft*p);

% subplot(311)
% plot(2*p*imag(sEst))
% subplot(312)
% plot(s(:,1))
% subplot(313)
% plot(abs(imag(sEst)-s(:,1)))

figure()
freq = 0:effectiveFs/(p*nfft):(effectiveFs-effectiveFs/(p*nfft));
subplot(121)
plot(freq/1e6,abs(fft(sEst,p*nfft)).^2)
subplot(122)
plot(pgram)

