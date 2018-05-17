function [xs,f] = sab_fft(x, sample_freq, cut_freq)
    n = numel(x);
    sf = fft(x);
    sfabs = abs(sf / n );
    xs = 2*sfabs( 1:(1 + n/2) );
    xs(1) = 0.5*xs(1);
    f = 0 : sample_freq/n : cut_freq;
    xs = xs(1:numel(f));
end