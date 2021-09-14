close all; 
clc;

%---------------Define Declaration---------------%
%Num of bits = 1024 (given)
nBits = 1024;

%Generating random binary digits
Data = randi([0 1],1,nBits);
%Converting to +1 and -1
Signal = 2 .* (Data - 0.5);

%Generating noise with mean = 0 and variance = 1
Noise = randn(1,nBits);

%Fixing SNR value to 10 dB
SNR = 10;
%SNR = 10log(S/N), where S = Signal Power, N = Noise Power
%Signal has unit power -> S = 1
%Therefore, using above equation, N = 0.1
noisePower = 0.1^(SNR/10);
%Adjusting Noise based on the SNR value
Noise = sqrt(noisePower) .* Noise;

%Calculating the received signal
Received = Signal + Noise;

%Fixing a threshold value for the received signal
Threshold = 0;
%Calculating the output based on the threshold level
Output = zeros(1,nBits);
for k = 1 : nBits
    if Received(k) > Threshold
        Output(k) = 1;
    else
        Output(k) = 0;
    end
end

%Computing the total number of errors
Error = 0;
for k = 1 : nBits
    if Output(k) ~= Data(k) 
        Error = Error + 1;
    end
end

%Calculating the bit error rate
% BER = (Total errors)/(Number of bits)
BER = double(Error)/double(nBits);
disp(BER)