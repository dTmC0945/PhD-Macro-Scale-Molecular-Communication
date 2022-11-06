function [inputMatrix, probabilityMatrix, dataPoints] = decodingSimulation2D(levelCounter, bitLength, symbolDuration, noisedSignal, normal, randomBits1, randomBits2)
    
    dataPoints          = zeros(1, bitLength);
    Decoded             = zeros(1, bitLength);                               % Decoded
    Cecoded             = zeros(1, bitLength);        
    inputAlphabet       = zeros(levelCounter, levelCounter);       % Input Alph.
    Alphabet            = zeros(levelCounter, levelCounter);       % Alpabet    
    inputMatrix         = zeros(0 + 1, levelCounter^2 + 1);
    probabilityMatrix   = zeros(levelCounter^2 + 1, levelCounter^2 + 1);

    
    % Generating discrete received values from the transmission

    for message_operator = 2:1:bitLength
                dataPoints(message_operator) = ...
                    noisedSignal(symbolDuration*message_operator)*normal...
                  + complex(levelCounter + 1, levelCounter + 1); 
    end

    % Generation of the symbol values for use in communications

    for b = 1:1:levelCounter^2 % Number of symbol values       
        Alphabet(b) = b - 1;        
    end

    % Naming of the symbol for each constellation values

    alphabetCounter = 1;

    for b = 1:1:levelCounter % In-phase value        
        for a = 1:1:levelCounter  % Out-phase value          
            inputAlphabet(b,a) = Alphabet(alphabetCounter);            
            alphabetCounter = alphabetCounter + 1;
        end
    end

    % Decoded Values

    dataPointsReal = real(dataPoints); % real
    dataPointsImag = imag(dataPoints); % imaginary


    for c = 2:1:bitLength    
        for b = 1:1:levelCounter   
            for a = 1:1:levelCounter 
                if dataPointsReal(c) < 3 && ...
                   dataPointsImag(c) < 3
                    Decoded(c) = inputAlphabet(1,1);
                elseif dataPointsReal(c) > 2*levelCounter - 1 && ...
                       dataPointsImag(c) > 2*levelCounter - 1
                    Decoded(c) = ...
                        inputAlphabet(levelCounter, levelCounter);
                elseif dataPointsReal(c) > 2*levelCounter - 1 && ...
                       dataPointsImag(c) < 3
                    Decoded(c) = ...
                        inputAlphabet(levelCounter, 1);
                elseif dataPointsReal(c) < 3 && ...
                       dataPointsImag(c) > 2*levelCounter - 1
                    Decoded(c) = ...
                        inputAlphabet(1, levelCounter);
                elseif dataPointsReal(c) > 2*b - 1 && ...
                       dataPointsReal(c) < 2*b + 1 && ...
                       dataPointsImag(c) > 2*levelCounter + 1
                    Decoded(c) = ...
                        inputAlphabet(b ,levelCounter);
                elseif dataPointsReal(c) > 2*b - 1 && ...
                       dataPointsReal(c) < 2*b + 1 && ...
                       dataPointsImag(c) < 3
                    Decoded(c) = ...
                        inputAlphabet(b ,1);
                elseif dataPointsImag(c) > 2*a - 1 && ...
                       dataPointsImag(c) < 2*a + 1 && ...
                       dataPointsReal(c) > 2*levelCounter + 1
                    Decoded(c) = ...
                        inputAlphabet(levelCounter, a);
                elseif dataPointsImag(c) > 2*a - 1 && ...
                       dataPointsImag(c) < 2*a + 1 && ...
                       dataPointsReal(c) < 3
                    Decoded(c) = inputAlphabet(1, a); 
                elseif dataPointsReal(c) > 2*b-1 && ...
                       dataPointsImag(c) > 2*a-1 && ...
                       dataPointsReal(c) < 2*b+1 && ...
                       dataPointsImag(c) < 2*a+1 
                    Decoded(c) = inputAlphabet(b,a);                 
                end
            end
        end

        for b = 1:1:levelCounter
            for a = 1:1:levelCounter
                if  randomBits1(c) >= (2*b ...
                        - 1) ...
                 && randomBits1(c) <= (2*b ...
                        + 1) ...
                 && randomBits2(c) >= (2*a ...
                        - 1) ...
                 && randomBits2(c) <= (2*a ...
                        + 1)                    
                        Cecoded(c) = inputAlphabet(b,a);
                end               
            end
        end        
    end     % decision algorithm

    for c = 1:1:bitLength     
        for b = 1:1:levelCounter^2 + 1
            if Cecoded(c) == (b - 1)    
                inputMatrix(b) = inputMatrix(b) + 1; 
            end

            for a = 1:1:levelCounter^2 + 1          
                if Decoded(c) == (b - 1) && Cecoded(c) == (a - 1)                    
                    probabilityMatrix(b,a) = probabilityMatrix(b,a) + 1;
                end
            end
        end
    end

end