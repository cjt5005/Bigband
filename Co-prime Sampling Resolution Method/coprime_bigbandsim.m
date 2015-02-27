close all
clear all

%% Receiver Parameters

effectiveBW = 500e6; %Nyquist Region for Fs = p*channelFs
p = [9;10]; % Row M contains the Channel M Undersampling Factor
nfft = 128;

effectiveFs = effectiveBW*2;
channelFs = effectiveFs./p;

timeDelays = [0;1/effectiveFs];


%% Signal Model

f = [185e6;210e6;185e6]; %signal frequencies
A = [1;1;1];
phi = rand(3,1)*2*pi;

% f = 200e6;
% A=1;
% phi = rand()*2*pi;

%Form the (N x K x 2) signal matrix where:
% The number of samples is N
% The number of channels is K
% The third dimension contains the non-delayed and delayed version of the
% signal

s = zeros(nfft,length(p),2);
for i = 1:length(p)
    [s(:,i,:),t] = sigGen(A,f,phi,effectiveFs/p(i),timeDelays,nfft);
end

% figure()
% for k=1:length(p)
%     hold on
%     subplot(length(p),1,k)
%     plot(s(:,k,1))
%     hold on
%     plot(s(:,k,2))
% end

%% Frequency Bucketization

% Bucketize the signal (i.e. Compute the periodogram of the resampled
% signal)
sF = fft(s,nfft);
for i = 1:length(p);
    freq(:,i) = 0:channelFs(i)/nfft:(channelFs(i)-channelFs(i)/nfft);
end

pgram = abs(sF(1:nfft/2,:,:)).^2;
plotFreq = freq(1:nfft/2,:)/1e6;

figure()
for k=1:length(p)
    hold on
    subplot(length(p),1,k)
    plot(plotFreq(:,k),pgram(:,k,1))
    hold on
    plot(plotFreq(:,k),pgram(:,k,2))
end

%% Estimate True Frequencies

N = [0:.5:floor(max(p)/2)]';

% Find the occupied buckets
[~, buckets1] = findpeaks([0;pgram(:,1,1);0]);
[~, buckets2] = findpeaks([0;pgram(:,2,1);0]);

%
fpos1 = plotFreq(buckets1,1);
fpos2 = plotFreq(buckets2,2);

Fs1 = channelFs(1)/1e6;
Fs2 = channelFs(2)/1e6;
f1h = fpos1;
f2h = fpos2;

N = 0:floor(p(1)/2);
M = 0:floor(p(2)/2);

f1 = repmat(f1h,1,length(N)) + repmat(N*Fs1,length(buckets1),1);
f2 = repmat(-f1h,1,length(N)) + repmat(N*Fs1,length(buckets1),1);
f3 = repmat(f2h,1,length(M)) + repmat(M*Fs2,length(buckets1),1);
f4 = repmat(-f2h,1,length(M)) + repmat(M*Fs2,length(buckets1),1);

% Comparison 1 - (f1,f3)
val1 = inf;
ind1 = 0;

for k = 1:size(f1,1)
    for i =1:length(f1)
        for j = 1:length(f3)
            if abs(f1(k,i)-f3(k,j))^2 < val1
                ind1 = i;
                val1 = abs(f1(k,i)-f3(k,j)).^2;
            end
        end
    end
    % Comparison 2 - (f1,f4)
    val2 = inf;
    ind2 = 0;
    for i =1:length(f1)
        for j = 1:length(f4)
            if abs(f1(k,i)-f4(k,j))^2 < val2
                ind2 = i;
                val2 = abs(f1(k,i)-f4(k,j)).^2;
            end
        end
    end
    % Comparison 3 - (f2,f3)
    val3 = inf;
    ind3 = 0;
    for i =1:length(f2)
        for j = 1:length(f3)
            if abs(f2(k,i)-f3(k,j))^2 < val3
                ind3 = i;
                val3 = abs(f2(k,i)-f3(k,j)).^2;
            end
        end
    end
    % Comparison 4 - (f2,f4)
    val4 = inf;
    ind4 = 0;
    for i =1:length(f2)
        for j = 1:length(f4)
            if abs(f2(k,i)-f4(k,j))^2 < val4
                ind4 = i;
                val4 = abs(f2(k,i)-f4(k,j)).^2;
            end
        end
    end

    [val1,ind1;val2,ind2;val3,ind3;val4,ind4]
    ind = [ind1;ind2;ind3;ind4];
    val = [val1;val2;val3;val4];
    [~,i] = min(val);

    if i<=2
        freqEst(k) = f1(ind(i));
    else
        freqEst(k) = f2(ind(i))
    end
end

freqEst
