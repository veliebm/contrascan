function [canonical4BOLD] = contrascan_canonical(onsets)

% onsets is a column vector of onset times

time = 0:0.01:750; % a new time vector that extends at least over the time of the session, for now we are not worried that it is sampled higher than the actual BOLD (TR = 2s)

contrast = zeros(1, length(time)); % this will be poplulated with the contrast cacnpincal HRF over time, to be correlated with BOLD

kernel = []; % This will be the kernel for the contrast (i.e. the shape of the contrast) 

for x = 1:410
    kernel(x) = exp(x./100);
end

for onset = 1: length(onsets)
    onsetindexvec(onset) = find(round(time.*100-onsets(onset)*100) == 0);
    contrast(onsetindexvec(onset):onsetindexvec(onset)+409) = kernel; 
end

onsets(1:4)                             % these are for testing that everything is as expected 
plot(time, contrast)                   % we should see the contrast function ramping up at the onset times
plot(time(1:10000), contrast(1:10000)) 

% now we could be done, but that stuff will never correlate with BOLD
% because it is at the wrong times (BOLD happens much later) and has this
% ridiculous sample rate, right, and we have BOLD every 2 secs. So now
% what? 

% we concolve the contrast time series with a canonical hemodynamic response function 
canonicalHRF = gampdf(0.0:0.01:18, 6); % 18 seconds long, gamma function, sampled at the same rate as time above
canonical = conv(contrast,canonicalHRF, 'same'); 
plot(0.0:0.01:18, canonicalHRF)
plot(canonical)

% now we resample to match actual TR
canonical4BOLD = resample(canonical,1,200); 
plot(canonical4BOLD)
