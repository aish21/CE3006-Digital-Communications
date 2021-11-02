close all;
clear all;
clear workspace;
clc;

% Phase 3: OOK Hamming Encoding

nBits = 1024;
codeword_len = 7; 
msgword_len = 4;

%initialisations for CBC encoding
genpoly = cyclpoly(codeword_len,msgword_len); %n=7, k=4
parmat = cyclgen(codeword_len,genpoly);
trt = syndtable(parmat);

encBits = nBits*codeword_len/msgword_len;                   %Number of bits for encoded signal
carrierFrequency = 10000;
carrierSignal = carrierFrequency * 16;
dataRate = 1000;
samplingPeriod = carrierSignal / dataRate;
[lowB, lowA] = butter(6,0.2);
amp = 1;

%time samples for signals
t = 0: 1/carrierSignal : nBits/dataRate;                        
tEnc = 0: 1/carrierSignal : encBits/dataRate;

SNR_dB = 0:5:50;
SNR = (10.^(SNR_dB/10));
modifySNR_dB = 5;

%initialisations for error counting
errorRateOOK = zeros(length(SNR));
errorRateOOKHam = zeros(length(SNR));
errorRateOOKCBC = zeros(length(SNR));

carrier = amp .* cos(2*pi*carrierFrequency*t);
carrierEnc = amp .* cos(2*pi*carrierFrequency*tEnc);

signalLength = carrierSignal*nBits/dataRate + 1;
signalLengthEnc = carrierSignal*encBits/dataRate + 1;

numRuns = 20;

for i = 1 : length(SNR)
    avgErrorOOK = 0;
    avgErrorOOKHam = 0;
    avgErrorOOKCBC = 0;
    
    for j = 1 : numRuns
        data = round(rand(1,nBits));
        dataStream = zeros(1, signalLength);
        
        hammingSignal = encode(data,codeword_len,msgword_len,'hamming/binary');
        dataStreamHam = zeros(1, signalLengthEnc);
        
        cbcSignal = encode(data,codeword_len,msgword_len,'cyclic/binary',genpoly);
        dataStreamCBC = zeros(1, signalLengthEnc);
        
        for k = 1: signalLength - 1
            dataStream(k) = data(ceil(k*dataRate/carrierSignal));
        end
        dataStream(signalLength) = dataStream(signalLength - 1);
        OOKSignal = carrier .* dataStream;

        for k = 1: signalLengthEnc - 1
            dataStreamHam(k) = hammingSignal(ceil(k*dataRate/carrierSignal));
            dataStreamCBC(k) = cbcSignal(ceil(k*dataRate/carrierSignal));
        end
     
        dataStreamHam(signalLengthEnc) = dataStreamHam(signalLengthEnc - 1);
        OOKSignalHam = carrierEnc .* dataStreamHam;

        dataStreamCBC(signalLengthEnc) = dataStreamCBC(signalLengthEnc - 1);
        OOKSignalCBC = carrierEnc .* dataStreamCBC;
        
        
        % Noise Generation
        OOKSignalPower = (norm(OOKSignal)^2)/signalLength;  
        OOKNoise = OOKSignalPower ./SNR(i);
        OOKNoise = sqrt(OOKNoise/2) .*randn(1,signalLength);
        
        OOKSignalPowerHam = (norm(OOKSignalHam)^2)/signalLengthEnc;  
		    OOKNoiseHam = OOKSignalPowerHam ./SNR(i);
		    OOKNoiseHam = sqrt(OOKNoiseHam/2) .*randn(1,signalLengthEnc);
        
        OOKSignalPowerCBC = (norm(OOKSignalCBC)^2)/signalLengthEnc;  
        OOKNoiseCBC = OOKSignalPowerCBC ./SNR(i);
        OOKNoiseCBC = sqrt(OOKNoiseCBC/2) .*randn(1,signalLengthEnc);
        
        % Transmit signal
		    OOKTransmit = OOKSignal + OOKNoise;
        OOKTransmitHam = OOKSignalHam + OOKNoiseHam;
        OOKTransmitCBC = OOKSignalCBC + OOKNoiseCBC;
        
        % Use non-coherent detection - square law detector
        sqLawOOK = OOKTransmit .* OOKTransmit;
        sqLawOOKHam = OOKTransmitHam .* OOKTransmitHam;
        sqLawOOKCBC = OOKTransmitCBC .* OOKTransmitCBC;
        
        % Pass this through the LP filter
        LPfilterOOK = filtfilt(lowB, lowA, sqLawOOK);
        LPfilterOOKHam = filtfilt(lowB, lowA, sqLawOOKHam);
        LPfilterOOKCBC = filtfilt(lowB, lowA, sqLawOOKCBC);

        % Sample the signal
        OOKsample = sample(LPfilterOOK, samplingPeriod, nBits);
        finalResultOOK = decision_device(OOKsample,nBits, amp/2);
        
        OOKsampleHam = sample(LPfilterOOKHam, samplingPeriod, encBits);
        finalResultOOKHam = decision_device(OOKsampleHam,encBits, amp/2);
        finalResultOOKDecodedHam = decode(finalResultOOKHam, codeword_len, msgword_len, 'hamming/binary');

        OOKsampleCBC = sample(LPfilterOOKCBC, samplingPeriod, encBits);
        finalResultOOKCBC = decision_device(OOKsampleCBC,encBits, amp/2);
        finalResultOOKDecodedCBC = decode(finalResultOOKCBC,codeword_len,msgword_len,'cyclic/binary',genpoly,trt);
        
       
        % Calculate Error
        errorOOK = 0;
        errorOOKHam = 0;
        errorOOKCBC = 0;
        for k = 1: nBits - 1
            if(finalResultOOK(k) ~= data(k))
                errorOOK = errorOOK + 1;
            end
            if(finalResultOOKDecodedHam(k) ~= data(k))
                errorOOKHam = errorOOKHam + 1;
            end
            if(finalResultOOKDecodedCBC(k) ~= data(k))
                errorOOKCBC = errorOOKCBC + 1;
            end
        end 
        avgErrorOOK = errorOOK + avgErrorOOK;
        avgErrorOOKHam = errorOOKHam + avgErrorOOKHam;
        avgErrorOOKCBC = errorOOKCBC + avgErrorOOKCBC;
    end
    
    % Variable Plots
    if (modifySNR_dB == SNR_dB(i))
        plotSignal = data;
        plotOOK = OOKSignal;
        plotOOKTransmit = OOKTransmit;
        plotLPfilterOOK = LPfilterOOK;
        plotfinalResultOOK = finalResultOOK;
        
        %HAMMING
        plotSignalHam = hammingSignal;
        plotOOKHam = OOKSignalHam;
        plotOOKTransmitHam = OOKTransmitHam;
        plotLPfilterOOKHam = LPfilterOOKHam;
        plotfinalResultOOKHam = finalResultOOKDecodedHam;

        %CBC
        plotSignalCBC = cbcSignal;
        plotOOKCBC = OOKSignalCBC;
        plotOOKTransmitCBC = OOKTransmitCBC;
        plotLPfilterOOKCBC = LPfilterOOKCBC;
        plotfinalResultOOKCBC = finalResultOOKDecodedCBC;
    end
    
    errorRateOOK(i) = (avgErrorOOK / numRuns)/nBits;
    errorRateOOKHam(i) = (avgErrorOOKHam/numRuns)/nBits;
    errorRateOOKCBC(i) = (avgErrorOOKCBC/numRuns)/nBits;
end

figure(1);
semilogy (SNR_dB,errorRateOOK,'r');
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
title("OOK - Encoding Comparison")
hold on

semilogy (SNR_dB,errorRateOOKHam,'b');
hold on
semilogy (SNR_dB,errorRateOOKCBC,'g');

leg = legend("OOK",'OOK-Hamming','OOK-CBC');
set(leg,'location','southwest');
hold off


figure(2);
subplot(511);plot(plotSignal);title('Data Generated: ');
subplot(512);plot(plotOOK,'k');title('Step - Modulation (OOK): ');
subplot(513);plot(plotOOKTransmit, 'k');title('Signal Recieved: ');
subplot(514);plot(plotLPfilterOOK, 'k');title('Step - Demodulation (OOK): ');
subplot(515);plot(plotfinalResultOOK);title('Decoded: ');

figure(3);
subplot(611);plot(plotSignal);title('Data Generated: ');
subplot(612);plot(plotSignalHam);title('Encoded Data - Hamming: ');
subplot(613);plot(plotOOKHam,'k');title('Step - Modulation (OOK): ');
subplot(614);plot(plotOOKTransmitHam, 'k');title('Signal Recieved: ');
subplot(615);plot(plotLPfilterOOKHam, 'k');title('Step - Demodulation (OOK): ');
subplot(616);plot(plotfinalResultOOKHam);title('Decoded: ');

figure(4);
subplot(611);plot(plotSignal);title('Data Generated: ');
subplot(612);plot(plotSignalCBC);title('Encoded Data - CBC: ');
subplot(613);plot(plotOOKCBC,'k');title('Step - Modulation (OOK): ');
subplot(614);plot(plotOOKTransmitCBC, 'k');title('Signal Recieved: ');
subplot(615);plot(plotLPfilterOOKCBC, 'k');title('Step - Demodulation (OOK): ');
subplot(616);plot(plotfinalResultOOKCBC);title('Decoded: ');


% Sampling and Decision Device Simulation
function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end

function binary_out = decision_device(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) > threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end
