% ========================= BEGIN CODE ====================================

clear variables
clc

% Variables ---------------------------------------------------------------

x = 0; v = inf; D = 1e-50; symbolDuration= 1;      

min = -10;  mid = 2;    max = 40;  

bitMax = 5; bitLength = 10000; 

% Memory Allocation -------------------------------------------------------

Error   = zeros((max-min)/2 + 1, bitMax);
Channel = zeros((max-min)/2 + 1, bitMax);
I_XY    = zeros((max-min)/2 + 1, bitMax);

% -------------------------------------------------------------------------

for bitCounter = 1:1:bitMax
    
    levelCounter = 2^(bitCounter);    % Number of levels in the channel
        
    counter = 1;
    
    for dB = min:mid:max                 % Counter operator (dB)

        SNR = 10^(dB/10);
            
        % Memory Allocation -----------------------------------------------

        H_YX            = zeros(levelCounter,levelCounter); % Contional(YX)
        H_Y             = zeros(1,levelCounter);            % Entropy  (Y)
        Enput           = zeros(levelCounter,levelCounter); % Probability 
    
        % Generating random bits ------------------------------------------
        
        randomBits          = randi([0 levelCounter-1],1,bitLength);  
                
        % Delegating the random bits to two distinct matrices -------------

        [bitSequence, ...
         bitRepetition]     = splitter(randomBits);
    
        % Simulating physical transmission --------------------------------

        Cm                  = transmission(D, v, x, ...
                                           bitSequence, ...
                                           bitRepetition, symbolDuration);
    
        % Adding noise to the signal --------------------------------------
        
        mu_a = 0; sigma_a   = sqrt(mean(Cm.^2)/SNR);

        CmNoised            = addNoise(Cm, mu_a, sigma_a);
             
        % -----------------------------------------------------------------

        [Decoded, Cecoded] = decoder(bitLength, ...
                                     levelCounter, ...
                                     randomBits, ...
                                     CmNoised, ...
                                     symbolDuration);     

        % Calculation of the error values per transmission ----------------

        [EValue, ...
         inputMatrix, ...
         Adj_Probability_Matrix, ...
         probabilityCheck, ...
         inputCheck]    = errorCalculator(Decoded, ...
                                          Cecoded, ...
                                          levelCounter, ...
                                          bitLength, ...
                                          counter, ...
                                          bitCounter);
        
        % Channel calculation function ------------------------------------
       
        [chValue, IValue] = channelCalculations(inputMatrix, ...
                                              Decoded, ...
                                              levelCounter, ...
                                              bitLength, ...
                                              Adj_Probability_Matrix, ...
                                              counter, ...
                                              bitCounter, ...
                                              dB);

        % Assigning the variables to matrices------------------------------
 
        Error(counter,bitCounter)   = EValue;   % Error Value
        I_XY(counter,bitCounter)    = IValue;   % Mutual Information
        Channel(counter,bitCounter) = chValue;  % Shannon Channel

        % -----------------------------------------------------------------
        
        counter = counter + 1;                  % counter iterator
        
    end
end

% -------------------------------------------------------------------------

plotDataChannelSER(Error, I_XY, Channel, bitCounter, min, mid, max);

% ======================== END OF CODE ====================================