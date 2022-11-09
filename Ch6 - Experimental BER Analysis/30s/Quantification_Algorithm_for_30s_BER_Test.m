clear all
clc

BER_5s = 'ber30sec (RAW).xlsx';

A = xlsread(BER_5s);

% Time vector conversion

T1 = datevec(A(1:44985,5));
T1_new = 60.*60.*T1(:,4) + 60.*T1(:,5) + T1(:,6);

% QUANTIFICATION ALGORITHM ------------------------------------------------

w = 1;

n = 0;
step = 1;

for i = 1:1:length(T1_new)

    if (T1_new(i) <= n + step) && (T1_new(i) >= n)

            Y(w,n + 1) = A(i,3);

    else
        
        n = n + step
        
        NNz = sum(squeeze(sum(Y~=0,3)));

        S = sum(Y);

    end
    
    w = w + 1;
    
end

% Data Plotting ----------------------------------------------------------

avg = S./NNz;

plot(avg)

% END OF ALGORITHM --------------------------------------------------------