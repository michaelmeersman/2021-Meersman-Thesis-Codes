function [Tu, Vel_RMS, power_mag, power_freq, ...
    IntegralLengthScale, DisssipationLengthScale, IntegralTimeScale, epsilon, ReynoldsShearStressUU] = ...
    calcTurbValues(nu, VELOCITY, SAMPLEFREQ, WINDOW_TYPE,  WINDOW_OVERLAP, NFFT)
% The purpose of this function is to do standard turbulence calculations
% from hotwire velocity data. The data should be in column format where
% each column corresponds to a different flow speed or test parameter.
%
% Inputs: REQUIRED: VELOCITY (as column),  SAMPLEFREQ of the data
%
% Optional: WINDOW_TYPE (default hanning, 5 Hz resolution) for pwelch, 
% WINDOW_OVERLAP (default 50%) for pwelch, NFFT (number of
% FFT points (default is the default size of the hanning window for 5 Hz
% resolution).
%
% Outputs: Tu (turbulence intensity as percent), RMS of velocity, the
% magintude from the power spectrum, the frequency from the power spectrum,
% the integral length scale, the dissipation length scale and the integral
% time scale from cross-covariance.
%
% Sources:
% https://www.mathworks.com/help/signal/ug/power-spectral-density-estimates-using-fft.html
% Integral Length Scale from Procedure for Determining Turbulence Length Scales Using Hotwire Anemometry NASA Page 6
% Roach - The generation of nearly isotropic turbulence by means of grids
% Eqn 15
% Dissapation Length Scale
% Roach - Eqn 8
% http://www.mit.edu/course/1/1.061/www/dream/SEVEN/SEVENTHEORY.PDF page
% 2-3
% and Hotwire Practical Guide Page 33 1/(N-1)
% https://www.mathworks.com/matlabcentral/answers/298598-how-to-calculate-integral-time-scales-in-matlab
% http://web.aeromech.usyd.edu.au//15afmc/proceedings/papers/AFMC00064.pdf
%
% Written by Mark Agate 2/27/2018: magate@email.arizona.edu

switch nargin
    case 0
        error('Not enough inputs');
    case 1
        error('Not enough inputs');
    case 2
        error('Not enough inputs');
    case 3
        WINDOW_TYPE = hanning(2^(ceil(log2((SAMPLEFREQ/5)))));
        WINDOW_OVERLAP = 2^(ceil(log2((SAMPLEFREQ/5))))/2;
        NFFT = 2^(ceil(log2((SAMPLEFREQ/5))));
    case 4
        WINDOW_OVERLAP = size(WINDOW_TYPE, 1)/2;
        NFFT = size(WINDOW_TYPE, 1);
    case 5
        NFFT = size(WINDOW_TYPE, 1);
    otherwise
        % good to proceed
end

% Begin calculating values
Vel_mean = mean(VELOCITY, 1); % take the column mean

% RMS and fluctuations
% http://www.mit.edu/course/1/1.061/www/dream/SEVEN/SEVENTHEORY.PDF page
% 2-3
% and Hotwire Practival Guide Page 33 1/(N-1)
Vel_RMS = zeros(1,size(VELOCITY, 2));
Vel_prime = zeros(size(VELOCITY));

for j = 1:size(VELOCITY,2)
    Vel_prime(:,j) = (VELOCITY(:,j) - Vel_mean(1,j));
    Vel_RMS(1,j) = sqrt(1/(size(VELOCITY,1)-1)*sum(Vel_prime(:,j).^2));
end

ReynoldsShearStressUU = mean(Vel_prime.^2); % 1D Reynoldsstresses
Tu = Vel_RMS./Vel_mean*100; % Tu as percent for each speed

% Power spectra
% https://www.mathworks.com/help/signal/ug/power-spectral-density-estimates-using-fft.html
power_mag = zeros(size(WINDOW_TYPE, 1)/2+1, size(VELOCITY,2));
power_freq = zeros(size(WINDOW_TYPE, 1)/2+1, size(VELOCITY,2));

for j = 1:size(VELOCITY, 2)
    [power_mag(:,j), power_freq(:,j)] = ...
        pwelch(Vel_prime(:,j), WINDOW_TYPE, WINDOW_OVERLAP, NFFT, SAMPLEFREQ);
%    [power_mag(:,j), power_freq(:,j)] = pwelch(Vel_prime(:,j), WINDOW_TYPE, WINDOW_OVERLAP, NFFT, SAMPLEFREQ,'power'); %- Lutz
end

% Integral Length Scale from
% Procedure for Determining Turbulence Length Scales Using Hotwire Anemometry
% NASA Page 6
% Roach - The generation of nearly isotropic turbulence by means of grids
% Eqn 15
Ef0 = zeros(1, size(power_mag, 2));
IntegralLengthScale = zeros(1, size(power_mag, 2));
DisssipationLengthScale = zeros(1, size(power_freq, 2));
Y = zeros(size(power_freq));
Z = zeros(1, size(power_freq, 2)-1);

for j = 1:size(power_mag, 2)
    Ef0(1,j) = mean(power_mag(2:100,j));
    IntegralLengthScale(1,j) = Ef0(1,j)*Vel_mean(1,j)/(4*std(VELOCITY(:,j))^2);
end

% Dissapation Length Scale
% Roach - Eqn 8
for j = 1:size(power_freq, 2)
    for i = 1:size(power_freq,1)
        Y(i,j) = power_freq(i,j)*power_freq(i,j)*power_mag(i,j)/SAMPLEFREQ; % freq*freq*Energy where Energy is Power/Fs
    end
    Z(1,j) = trapz(power_freq(:,j), Y(:,j));
    Z(1,j) = Z(1,j)*2*pi^2/(Vel_mean(1,j)^2*mean(Vel_prime(:,j))^2);
    DisssipationLengthScale(1,j) = sqrt(1/Z(1,j));
end

% integral time scales from XCOV
% https://www.mathworks.com/matlabcentral/answers/298598-how-to-calculate-integral-time-scales-in-matlab
% http://web.aeromech.usyd.edu.au//15afmc/proceedings/papers/AFMC00064.pdf
% https://www.mathworks.com/help/signal/ref/xcov.html?s_tid=srchtitle
IntegralTimeScale = zeros(1, size(VELOCITY, 2));

for j = 1:size(Vel_prime, 2)
    [acor, ~] = xcorr(Vel_prime(:,j)-mean(Vel_prime(:,j)), 'coeff'); % auto correlation of time series scaled to 1
    n = length(Vel_prime(:,j));
    acor = acor(n:end); %only positive autocorrelation since it is two-sided
    id = find(acor<0, 1, 'first'); % first zero-crossing
    IntegralTimeScale(1,j) = trapz(linspace(0, (1/SAMPLEFREQ)*(id-1), id-1), acor(1:id-1)); % integral time scale (-1 accounts for last index where autocorrelation is less than zero)
end

epsilon = nu./IntegralTimeScale.^2;
end










