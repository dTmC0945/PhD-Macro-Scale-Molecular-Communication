% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- Experimental Bit Error Rate by Daniel T. McGuiness -------------------
% -------------------------------------------------------------------------

clear variables
clc

filename = 'ber5sec (Quantised).xlsx';

A = xlsread(filename);



duration = 5;
bitLength = 200;
errorCounter = 0;

discreteTime = zeros(1,bitLength);
discreteData = zeros(1,bitLength);
decodedValues = zeros(1,20);

for i= 1:bitLength

    discreteData(i) = A(18 + duration*(i)); 
    discreteTime(i) =   18 + duration*(i);

end

% Relative Decoding: the decoding is based on the preceeding value, if the
% value is less than its previous symbol value designate it as 0, else 
% designate it as 1.

for n = 1:bitLength
    
    if n == 1
        
        if discreteData(n) < 1 * 10^-10
        
        decodedValues(n) = 0;
        
        end
        
    elseif discreteData(n) > discreteData(n-1)
            
        if abs(discreteData(n-1) - discreteData(n)) < 0.1 * discreteData(n)
            
            decodedValues(n) = decodedValues(n-1);
        
        else
            
            decodedValues(n) = 1;
        
        end
    
    elseif discreteData(n) < discreteData(n-1)
        
        if abs(discreteData(n-1) - discreteData(n)) < 0.1 * discreteData(n)
            
            decodedValues(n) = decodedValues(n-1);
        
        else
        
            decodedValues(n) =  0;
        
        end
    end
end
        
experimentalTransmittedSignal =    [0 1 0 1	1 1 0 1 0 1	0 1	0 0	0 0	0 0	1 1	0 0	1 0	1 1	0 0	1 1	0 0	...
                                    1 0	1 0	1 0	0 1	0 0	0 1	0 1	1 0	0 0	0 0	1 1	1 1	1 0	0 0	1 1	1 1	...
                                    0 0	0 1	1 0	1 1	0 1	0 0	0 1	1 0	1 0	0 0	0 1	0 0	1 0	1 0	0 0	1 1 ...
    	                            1 0	0 1	0 1	1 0	0 0	0 1	0 1	0 1	1 0	1 1	0 0	1 1	0 1	0 0	1 1	0 0	...
                                    0 1	1 1	1 1	0 0	1 1	1 1	0 1	0 0	1 1	1 1	1 1	0 0	1 0	0 1	0 0	1 0 ...
    	                            1 1	1 0	0 0	1 0	1 1	1 0	0 0	1 0	0 1	0 0	1 1	1 1	0 1	1 0	0 1	0 1	...
                                    1 1	1 1	0 1	1 0];

for u = 1:bitLength

    if decodedValues(u) == experimentalTransmittedSignal(u)
        
        fprintf('%d is correct \n'  , u);
        
    else
        
        fprintf('%d is wrong \n'    , u);
        errorCounter = errorCounter + 1;
        
    end    
end

% ======================== END OF CODE ====================================