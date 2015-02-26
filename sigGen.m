function [s,t] = sigGen(A,f,phi,Fs,Td,N)
%% Description
%This function generates N samples of a M-tone sinusoidal signal at K time
%delays.
%
%Inputs:
% A   (Mx1): Amplitude
% f   (Mx1): Frequency in Hz
% phi (Mx1): Phase Offset in Radians
% Fs  (1x1): Sampling frequency in Hz
% Td  (Kx1): Time delay in seconds
%
%Outputs:
% s (KxN): Output Signal
% t (1xN): Time vector

%% Main

M = size(A,1);
K = size(Td,1);

t = repmat((0:1:N-1)*1/Fs,K,1)+repmat(Td,1,N);
A = repmat(A,1,N);
phi = repmat(phi,1,N);

s = zeros(M,N,K);
for i=1:K
    s(:,:,i) = A.*sin(2*pi*f*t(i,:)+phi);
end
s = squeeze(sum(s,1));

end

