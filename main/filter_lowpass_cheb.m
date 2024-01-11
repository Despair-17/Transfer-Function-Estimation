% function Hd = filter_lowpass_cheb
%FILTER_LOWPASS_CHEB Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.13 and DSP System Toolbox 9.15.
% Generated on: 05-May-2023 19:54:35

% Chebyshev Type I Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in Hz.
Fs = 31250;  % Sampling Frequency

Fpass = 500;         % Passband Frequency
Fstop = 700;         % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 80;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'cheby1', 'MatchExactly', match);

%%
% Calculate the frequency response of the filter
[H, W] = freqz(Hd);

% Convert frequency to Hz
W = W / (2*pi) * Fs;

% Plot magnitude response (in dB) and phase response
subplot(2,1,1)
plot(W, 20*log10(abs(H)))
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Filter Frequency Response')
grid on

subplot(2,1,2)
plot(W, angle(H))
xlabel('Frequency (Hz)')
ylabel('Phase (rad)')
grid on