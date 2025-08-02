% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- Experimental Bit Error Rate by Daniel T. McGuiness -------------------
% -------------------------------------------------------------------------

function errorCounter = experimentalBERDataAnalysis(filename, duration, bitLength, start)

    
    A = xlsread(filename);
        
    errorCounter = 0;
    
    discreteTime    = zeros(1, bitLength);
    discreteData    = zeros(1, bitLength);
    decodedValues   = zeros(1, bitLength);
    
    for i= 1:bitLength
    
        discreteData(i) = A(start + duration*i); 
        discreteTime(i) =   start + duration*i;
    
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
    
        if decodedValues(u) ~= experimentalTransmittedSignal(u)
                         
            errorCounter = errorCounter + 1;
            
        end    
    end

end
% ======================== END OF CODE ====================================