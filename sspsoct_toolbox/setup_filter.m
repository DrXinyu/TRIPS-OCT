function hf = setup_filter(afrq)

    Fstop1 = afrq-0.02;  % First Stopband Frequency
    Fpass1 = afrq-0.01;  % First Passband Frequency
    Fpass2 = afrq+0.01;  % Second Passband Frequency
    Fstop2 = afrq+0.02;  % Second Stopband Frequency
    Astop1 = 100;    % First Stopband Attenuation (dB)
    Apass  = 1;     % Passband Ripple (dB)
    Astop2 = 100;    % Second Stopband Attenuation (dB)


    h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
        Fpass2, Fstop2, Astop1, Apass, Astop2);

    Hd = design(h, 'equiripple', ...
        'MinOrder', 'any');

    [h,t] = impz(Hd);
    hf.hfft = fft(h,4096);
    hf.states = Hd.States;
end

