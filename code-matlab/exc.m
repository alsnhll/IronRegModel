 function e=exc(bi,e0)
        %e=e0;
        B2=0.001; %"fermi temperature" of excretion versus body iron level. baseline 0.001
        e=e0./(1+exp(-bi/B2));
    end