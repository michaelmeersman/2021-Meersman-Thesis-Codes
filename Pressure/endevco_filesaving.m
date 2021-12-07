clear; close

Re = 300;
for AOA = 0:10
    
    endevco = importdata(['C:\Users\micha\Documents\1.RESEARCH\'...
        '2.WIND TUNNEL\WIND_TUNNEL_DATA\2020_07_06_Pressure_Endevco\07_06_Endevco\'...
        '07_06_' num2str(Re) 'k_Endevco_' num2str(AOA) '.txt'],'\t',4);
    
    fftsize = 2^12;
    WINDOW_TYPE = hanning(fftsize);
    WINDOW_OVERLAP = fftsize/2;
    SAMPLEFREQ = 32760;
    
    lowpass = 2000;
    d1 = designfilt('lowpassiir','FilterOrder',8, ...
        'HalfPowerFrequency',lowpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
    
    E1 = endevco.data(:,2);
    E2 = endevco.data(:,3);
    E3 = endevco.data(:,4);
    
    E1 = filtfilt(d1,E1);
    E2 = filtfilt(d1,E2);
    E3 = filtfilt(d1,E3);
    
    
    test = [E1 E2 E3];
    for i =1:3
        
        sample = test(:,i);
        [Tu, Vel_RMS, power_mag, power_freq] = ...
            calcTurbValues(1.65E-5,sample,SAMPLEFREQ,WINDOW_TYPE,WINDOW_OVERLAP,fftsize);
        
        
        power_mags(i,:) = power_mag;
        power_freqs(i,:) = power_freq;
        
        
    end
    
    save(['Endevco_RE_' num2str(Re) 'k_AOA_' num2str(AOA) '_filt.mat'])
end