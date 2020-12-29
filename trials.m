FILE = 'eric.wav';
[y, Fs]= audioread(FILE);
Y = fftshift(fft(y));

%low pass filter
N = length(y);
numberOfOnes = floor(N * 8000 / Fs);
numberOfZeros = floor((N - numberOfOnes) / 2);
remainder = mod((N - numberOfOnes), 2);
rect = ones(numberOfOnes, 1);
filter = padarray(rect, numberOfZeros, 'pre');
filter = padarray(filter, numberOfZeros + remainder, 'post');
filteredSignalSpectrum = filter .* Y;
filteredSignal = real(ifft(ifftshift(filteredSignalSpectrum)));
sound(filteredSignal);

message = filteredSignal;
fc = 100000;
carrierFs = 5 * fc;
resampledMessage = resample(message, carrierFs, Fs);
vectorLength = length(resampledMessage);
time = linspace(0, vectorLength/carrierFs, vectorLength).';
carrier = cos(2 * pi * fc * time);
DSBSC = carrier .* resampledMessage;

figure(1);
plot(abs(fftshift(fft(DSBSC))));

A = max(abs(resampledMessage));
DSBTC = A.*(1+0.5.*resampledMessage).*carrier;
figure(2);
plot(abs(fftshift(fft(DSBTC))));

%%demodulization ?
reveivedSignal = DSBSC.*carrier; %same carrier assuming no error in phase or frequency (coherent)
N = length(reveivedSignal);
numberOfOnes = floor(N * 8000 / Fs);
numberOfZeros = floor((N - numberOfOnes) / 2);
remainder = mod((N - numberOfOnes), 2);
rect = ones(numberOfOnes, 1);
filter = padarray(rect, numberOfZeros, 'pre');
filter = padarray(filter, numberOfZeros + remainder, 'post');
filteredSignalSpectrum = filter .* reveivedSignal;
filteredSignalDemod = real(ifft(ifftshift(filteredSignalSpectrum)));
demodulatedSignal = resample(filteredSignalDemod, Fs, 5*fc);