%% EMG ANALYSIS - starter script (Y2 PW1)
% Loads an EMG recording (.mat with Time and Signal), filters it, works out
% the envelope, detects each contraction, and measures:
%   - rest vs contraction, number of contractions
%   - onset, offset and duration of each contraction
%   - relative force (envelope amplitude) of each contraction
%   - muscle fatigue (median frequency over time) - RAW signals only
%
% HOW TO USE
%   1. Put this file and your recording in the same folder.
%   2. Change "filename" and "signalType" below.
%   3. Press Run. Results print in the Command Window and figures open.
% Needs the Signal Processing Toolbox (butter, filtfilt).

clear; close all; clc;

%% 1. LOAD THE SIGNAL
filename   = 'op_test2.mat';  % <-- change to your file
signalType = 'ENV';     % 'RAW' if the MyoWare RAW pin is wired to A0,
                        % 'ENV' if the ENV (envelope) pin is wired to A0

d  = load(filename);
t  = d.Time(:);
x  = d.Signal(:);
fs = round(1 / mean(diff(t)));             % sampling frequency (Hz)
fprintf('Loaded %s: %d samples, %.1f s, fs = %d Hz\n', filename, numel(x), numel(x)/fs, fs);

x = x - mean(x);   % remove the DC offset (the Arduino reads 0-1023, centred ~512)

%% 2. FILTER AND ENVELOPE
% RAW: the useful EMG band is about 20-450 Hz.
%   band-pass 20-450 Hz -> removes movement artefacts and baseline drift
%                          (below 20 Hz) and high-frequency noise
%   band-stop 50 Hz     -> removes UK mains interference
%   then rectify (absolute value) and take a moving RMS over 100 ms to get
%   the envelope - a smooth line showing how hard the muscle is working.
% ENV: the MyoWare has already rectified and smoothed it, so only light
%   smoothing is needed.
win = round(0.1 * fs);                     % 100 ms window

if strcmpi(signalType, 'RAW')
    if fs < 900
        warning('fs = %d Hz is too low for RAW EMG (needs ~1000 Hz). Upper cut-off reduced.', fs);
    end
    fHigh = min(450, 0.45*fs);             % keep below the Nyquist limit (fs/2)
    [b,  a ] = butter(4, [20 fHigh] / (fs/2), 'bandpass');
    [bn, an] = butter(2, [48 52] / (fs/2), 'stop');
    y = filtfilt(b, a, x);
    y = filtfilt(bn, an, y);
    env = sqrt(movmean(y.^2, win));        % moving RMS envelope
else
    y = x;
    env = movmean(x, win);
    env = env - min(env);                  % lowest point = 0
end

figure('Name', 'Filtering');
subplot(2,1,1); plot(t, x); grid on;
title('Raw signal (offset removed)'); xlabel('Time (s)'); ylabel('Amplitude (ADC counts)');
subplot(2,1,2); plot(t, y); hold on; plot(t, env, 'r', 'LineWidth', 1.5); grid on;
title('Filtered signal and envelope'); xlabel('Time (s)'); ylabel('Amplitude (ADC counts)');
legend('Filtered', 'Envelope (100 ms RMS)');

%% 3. REST LEVEL AND THRESHOLD
% Rest level  = average of the quietest 20% of the envelope
% Active level = average of the strongest 5% of the envelope
% The muscle counts as contracting when the envelope is more than 20% of
% the way from rest to active.
s = sort(env);
restLevel   = mean(s(1 : round(0.20*end)));
activeLevel = mean(s(round(0.95*end) : end));
thr = restLevel + 0.2 * (activeLevel - restLevel);

%% 4. DETECT CONTRACTIONS
% 1. mark every sample above the threshold as "active"
% 2. fill gaps shorter than 0.2 s (one contraction that briefly dipped)
% 3. ignore bursts shorter than 0.2 s (twitches or noise spikes)
active = env > thr;
minLen = round(0.2 * fs);

[on, off] = edges(active);
for k = 1:numel(on) - 1                    % fill short gaps
    if on(k+1) - off(k) < minLen
        active(off(k) : on(k+1)) = true;
    end
end
[on, off] = edges(active);
keep = (off - on) >= minLen;               % drop short bursts
for k = find(~keep)'
    active(on(k) : off(k)) = false;
end
on = on(keep);  off = off(keep);

nC = numel(on);

%% 5. MEASURE EACH CONTRACTION
dur  = (off - on) / fs;                    % duration (s)
amp  = zeros(nC, 1);                       % mean envelope = relative force
for k = 1:nC
    amp(k) = mean(env(on(k) : off(k)));
end
relForce = 100 * amp / max([amp; eps]);    % % of the strongest contraction
restGap  = (on(2:end) - off(1:end-1)) / fs; % rest between contractions (s)

%% 6. FATIGUE - MEDIAN FREQUENCY OVER TIME (RAW only)
% As a muscle tires, its signal shifts to lower frequencies. In 0.5 s
% windows, find the frequency that splits the power spectrum in half
% (the median frequency). A falling trend during a held contraction
% indicates fatigue.
mdfT = [];  mdf = [];  slope = NaN;
if strcmpi(signalType, 'RAW') && nC > 0
    [~, longest] = max(dur);               % analyse the longest contraction
    seg = y(on(longest) : off(longest));
    L = round(0.5 * fs);                   % window length
    for st = 1 : L : numel(seg) - L + 1
        w  = seg(st : st + L - 1);
        P  = abs(fft(w)).^2;
        P  = P(1 : floor(L/2));            % one-sided spectrum
        fr = (0 : floor(L/2) - 1)' * fs / L;
        c  = cumsum(P);
        mdf(end+1, 1)  = fr(find(c >= c(end)/2, 1));
        mdfT(end+1, 1) = t(on(longest)) + (st + L/2) / fs;
    end
    if numel(mdf) >= 3
        p = polyfit(mdfT - mdfT(1), mdf, 1);
        slope = p(1);                      % Hz per second
    end
end

%% 7. RESULTS
fprintf('\n----- RESULTS -----\n');
fprintf('Rest level / threshold : %.1f / %.1f\n', restLevel, thr);
fprintf('Contractions detected  : %d\n', nC);
fprintf('Time contracting       : %.0f%% of the recording\n', 100*mean(active));
for k = 1:nC
    fprintf('  #%d  onset %5.2f s  offset %5.2f s  duration %4.2f s  force %3.0f%%\n', ...
            k, t(on(k)), t(off(k)), dur(k), relForce(k));
end
if nC > 1
    fprintf('Mean rest between      : %.2f s\n', mean(restGap));
end
if ~isnan(slope)
    fprintf('Median frequency       : %.0f Hz -> %.0f Hz (slope %.2f Hz/s)\n', mdf(1), mdf(end), slope);
    if slope < -0.5
        fprintf('Fatigue                : median frequency falling - signs of fatigue\n');
    else
        fprintf('Fatigue                : no clear fatigue trend\n');
    end
end

figure('Name', 'Contractions');
plot(t, env, 'k'); hold on; grid on;
yl = [0, 1.1*max(env)];
for k = 1:nC
    patch(t([on(k) off(k) off(k) on(k)]), yl([1 1 2 2]), [1 0.85 0.85], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t, env, 'k');
plot(t([1 end]), [thr thr], 'r--');
plot(t(on), env(on), 'g>', 'MarkerFaceColor', 'g');
plot(t(off), env(off), 'm<', 'MarkerFaceColor', 'm');
ylim(yl);
title(sprintf('%d contractions detected', nC));
xlabel('Time (s)'); ylabel('Envelope (ADC counts)');

if ~isempty(mdf)
    figure('Name', 'Fatigue');
    plot(mdfT, mdf, 'o-'); hold on; grid on;
    if ~isnan(slope), plot(mdfT, polyval(p, mdfT - mdfT(1)), 'r--'); end
    title('Median frequency during the longest contraction');
    xlabel('Time (s)'); ylabel('Median frequency (Hz)');
end

%% HELPER FUNCTION
function [on, off] = edges(active)
% start and end sample of each run of "true" values
    d   = diff([false; active(:); false]);
    on  = find(d == 1);
    off = find(d == -1) - 1;
end