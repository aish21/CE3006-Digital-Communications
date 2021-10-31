close all;
clc;

% Phase 3: OOK Hamming Encoding

nBits = 1024;
codeword_len = 7; 
msgword_len = 4;
genpoly = cyclpoly(codeword_len,msgword_len); %n=7, k=4
parmat = cyclgen(codeword_len,genpoly);
trt = syndtable(parmat);
encBits = nBits*codeword_len/msgword_len; %need to find values for cycle block code
carrierFrequency = 10000;
carrierSignal = carrierFrequency * 16;
dataRate = 1000;
samplingPeriod = carrierSignal / dataRate;
[lowB, lowA] = butter(6,0.2);
amp = 1;
t = 0: 1/carrierSignal : nBits/dataRate;
tEnc = 0: 1/carrierSignal : encBits/dataRate;
SNR_dB = 0:5:50;
SNR = (10.^(SNR_dB/10));
modifySNR_dB = 5;
errorRateOOK = zeros(length(SNR));
errorRateOOKenc = zeros(length(SNR));
carrier = amp .* cos(2*pi*carrierFrequency*t);
carrierEnc = amp .* cos(2*pi*carrierFrequency*tEnc);
signalLength = carrierSignal*nBits/dataRate + 1;
signalLengthEnc = carrierSignal*encBits/dataRate + 1;
numRuns = 20;

for i = 1 : length(SNR)
    avgErrorOOK = 0;
    avgErrorOOKEnc = 0;
    
    for j = 1 : numRuns
        data = round(rand(1,nBits));
        %hammingSignal = encode(data,7,4,'hamming/binary');
        cbcSignal = encode(data,codeword_len,msgword_len,'cyclic/binary',genpoly);
        dataStreamEnc = zeros(1, signalLengthEnc);
        
        for k = 1: signalLengthEnc - 1
            %dataStreamEnc(k) = hammingSignal(ceil(k*dataRate/carrierSignal));
            dataStreamEnc(k) = cbcSignal(ceil(k*dataRate/carrierSignal));
        end
        
        dataStreamEnc(signalLengthEnc) = dataStreamEnc(signalLengthEnc - 1);
        OOKSignalEnc = carrierEnc .* dataStreamEnc;
        
        dataStream = zeros(1, signalLength);
        for k = 1: signalLength - 1
            %dataStream(k) = hammingSignal(ceil(k*dataRate/carrierSignal));
            dataStreamEnc(k) = cbcSignal(ceil(k*dataRate/carrierSignal));
        end
        
        dataStream(signalLength) = dataStream(signalLength - 1);
        OOKSignal = carrier .* dataStream;
        
        % Noise Generation
        OOKSignalPower = (norm(OOKSignal)^2)/signalLength;  
		OOKNoise = OOKSignalPower ./SNR(i);
		OOKNoise = sqrt(OOKNoise/2) .*randn(1,signalLength);
        
        OOKSignalPowerEnc = (norm(OOKSignalEnc)^2)/signalLengthEnc;  
		OOKNoiseEnc = OOKSignalPowerEnc ./SNR(i);
		OOKNoiseEnc = sqrt(OOKNoiseEnc/2) .*randn(1,signalLengthEnc);
        
        % Transmit signal
		OOKTransmit = OOKSignal + OOKNoise;
        OOKTransmitEnc = OOKSignalEnc + OOKNoiseEnc;
        
        % Use non-coherent detection - square law detector
        sqLawOOK = OOKTransmit .* OOKTransmit;
        sqLawOOKEnc = OOKTransmitEnc .* OOKTransmitEnc;
        
        % Pass this through the LP filter
        LPfilterOOK = filtfilt(lowB, lowA, sqLawOOK);
        LPfilterOOKEnc = filtfilt(lowB, lowA, sqLawOOKEnc);

        % Sample the signal
        OOKsample = sample(LPfilterOOK, samplingPeriod, nBits);
        finalResultOOK = decision_device(OOKsample,nBits, amp/2);
        
        OOKsampleEnc = sample(LPfilterOOKEnc, samplingPeriod, encBits);
        finalResultOOKEnc = decision_device(OOKsampleEnc,encBits, amp/2);
        %finalResultOOKDecoded = decode(finalResultOOKEnc, 7, 4, 'hamming/binary');
        finalResultOOKDecoded = decode(finalResultOOKEnc,codeword_len,msgword_len,'cyclic/binary',genpoly,trt);
        % Calculate Error
        errorOOK = 0;
        errorOOKEnc = 0;
        for k = 1: nBits - 1
            if(finalResultOOK(k) ~= data(k))
                errorOOK = errorOOK + 1;
            end
        end 
        avgErrorOOK = errorOOK + avgErrorOOK;
        
        for k = 1: nBits - 1
            if(finalResultOOKDecoded(k) ~= data(k))
                errorOOKEnc = errorOOKEnc + 1;
            end
        end 
        avgErrorOOKEnc = errorOOKEnc + avgErrorOOKEnc;
    end
    
    % Variable Plots
    if (modifySNR_dB == SNR_dB(i))
        plotSignal = data;
        plotOOK = OOKSignal;
        plotOOKTransmit = OOKTransmit;
        plotLPfilterOOK = LPfilterOOK;
        plotfinalResultOOK = finalResultOOK;
    end
    errorRateOOK(i) = (avgErrorOOK / numRuns)/nBits;
    
    if (modifySNR_dB == SNR_dB(i))
        %plotSignalEnc = hammingSignal;
        plotSignalEnc = cbcSignal;
        plotOOKEnc = OOKSignalEnc;
        plotOOKTransmitEnc = OOKTransmitEnc;
        plotLPfilterOOKEnc = LPfilterOOKEnc;
        plotfinalResultOOKEnc = finalResultOOKDecoded;
    end
    errorRateOOKEnc(i) = (avgErrorOOKEnc/numRuns)/nBits;
end

% Plotting Error + OOK
%figure(1);
%semilogy (SNR_dB,errorRateOOK,'k-*');
%hold off

figure(1);
semilogy (SNR_dB,errorRateOOKEnc,'k-*');
hold off

%figure(3);
%subplot(511);title('Data Generated: ');plot(plotSignal);
%subplot(512);title('Step - Modulation (OOK): ');plot(plotOOK,'k');
%subplot(513);title('Signal Recieved: ');plot(plotOOKTransmit, 'k')
%subplot(514);title('Step - Demodulation (OOK): ');plot(plotLPfilterOOK, 'k');
%subplot(515);title('Decoded: ');plot(plotfinalResultOOK);

figure(2);
subplot(511);title('Data Generated: ');plot(plotSignalEnc);
subplot(512);title('Step - Modulation (OOK): ');plot(plotOOKEnc,'k');
subplot(513);title('Signal Recieved: ');plot(plotOOKTransmitEnc, 'k')
subplot(514);title('Step - Demodulation (OOK): ');plot(plotLPfilterOOKEnc, 'k');
subplot(515);title('Decoded: ');plot(plotfinalResultOOKEnc);

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