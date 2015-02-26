function [ pairs ] = findFreqPairs(f,Fs,p)
%This function returns the possible frequency pairs (f,f') for a given
%input frequency f, sampling rate Fs, and undersampling factor p.

f = 60;
Fs = 200;
p=9;

N = 1:floor(p/2);

falias = abs(f-N*Fs);
falias2 = abs(-f-N*Fs);
flist = [f falias falias2];
pairs = nchoosek(flist,2);

end

