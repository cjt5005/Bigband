function [ output_args ] = resolveCollision(pairs,b,td)
%RESOLVECOLLISION Summary of this function goes here
%   Detailed explanation goes here

for k = 1:length(pairs)
    f1 = pairs(k,1);
    f2 = pairs(k,2);
    A = [exp(2i*pi*f1*td) exp(2i*pi*f2*td)];
    x = pinv(A)*b;
    mse(k) = mean(abs(A*x-b).^2);
end

end

