function [s] = sigGen(A,f,Fs,Td,N)
%% Description
%This function generates N samples of a monotone sinusoidal signal with parameters:
% A: Amplitude
% f: frequency in Hz
% Fs: Sampling frequency in Hz
% Td: Time delay in seconds
%%
T = 1/f;
phi = Td/T*2*pi;
h = dsp.SineWave(A,f,phi,'SampleRate',Fs,'SamplesPerFrame',N);
s = step(h);

end

