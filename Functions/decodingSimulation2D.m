function [inputMatrix, probabilityMatrix] = decodingSimulation2D(levelCounter, bitLength, symbolDuration, noisedSignal, normal)
    
           % Generating discrete received values from the transmission

        for message_operator = 2:1:bitLength*symbolDuration
                    Rn(message_operator) = noisedSignal(symbolDuration*message_operator)*normal + (levelCounter + 1)*(1 + 1i); 
        end

        % Generation of the symbol values for use in communications

        for a = 1:1:levelCounter^2 % Number of symbol values       
            Alphabet(a) = a - 1;        
        end

        % Naming of the symbol for each constellation values

        alpha_counter = 1;

        for a = 1:1:levelCounter % In-phase value        
            for b = 1:1:levelCounter  % Out-phase value          
                Input_Alphabet(a,b) = Alphabet(alpha_counter);            
                alpha_counter = alpha_counter + 1;
            end
        end

        Rp = Prime_seed_1 + 1j*Prime_seed_2;

        % Decoded Values

        Rn_re = real(Rn); % real
        Rn_im = imag(Rn); % imaginary


        for dec = 2:1:bitLength    
            for a = 1:1:levelCounter   
                for b = 1:1:levelCounter 
                    if Rn_re(dec) < 3 && ...
                       Rn_im(dec) < 3
                        Decoded(dec) = Input_Alphabet(1,1);
                    elseif Rn_re(dec) > 2*levelCounter - 1 && ...
                           Rn_im(dec) > 2*levelCounter - 1
                        Decoded(dec) = ...
                            Input_Alphabet(levelCounter  ,levelCounter  );
                    elseif Rn_re(dec) > 2*levelCounter - 1 && ...
                           Rn_im(dec) < 3
                        Decoded(dec) = ...
                            Input_Alphabet(levelCounter  ,1);
                    elseif Rn_re(dec) < 3 && ...
                           Rn_im(dec) > 2*levelCounter - 1
                        Decoded(dec) = ...
                            Input_Alphabet(1  ,levelCounter);
                    elseif Rn_re(dec) > 2*a - 1 && ...
                           Rn_re(dec) < 2*a + 1 && ...
                           Rn_im(dec) > 2*levelCounter + 1
                        Decoded(dec) = ...
                            Input_Alphabet(a  ,levelCounter);
                    elseif Rn_re(dec) > 2*a - 1 && ...
                           Rn_re(dec) < 2*a + 1 && ...
                           Rn_im(dec) < 3
                        Decoded(dec) = ...
                            Input_Alphabet(a  ,1);
                    elseif Rn_im(dec) > 2*b - 1 && ...
                           Rn_im(dec) < 2*b + 1 && ...
                           Rn_re(dec) > 2*levelCounter + 1
                        Decoded(dec) = ...
                            Input_Alphabet(levelCounter  ,b);
                    elseif Rn_im(dec) > 2*b - 1 && ...
                           Rn_im(dec) < 2*b + 1 && ...
                           Rn_re(dec) < 3
                        Decoded(dec) = Input_Alphabet(1  ,b); 
                    elseif Rn_re(dec) > 2*a-1 && ...
                           Rn_im(dec) > 2*b-1 && ...
                           Rn_re(dec) < 2*a+1 && ...
                           Rn_im(dec) < 2*b+1 
                        Decoded(dec) = Input_Alphabet(a,b);                 
                    end
                end
            end

             for a = 1:1:levelCounter
                for b = 1:1:levelCounter
                    if  Prime_seed_1(dec) >= (2*a ...
                            - 1) ...
                     && Prime_seed_1(dec) <= (2*a ...
                            + 1) ...
                     && Prime_seed_2(dec) >= (2*b ...
                            - 1) ...
                     && Prime_seed_2(dec) <= (2*b ...
                            + 1)                    
                            Cecoded(dec) = Input_Alphabet(a,b);
                    end               
                end
             end        
        end     % decision algorithm

        for c = 1:1:length(Decoded)     
            for a = 1:1:levelCounter^2 + 1

                if Cecoded(c) == (a - 1)    
                    inputMatrix(a) = inputMatrix(a) + 1; 
                end

                for b = 1:1:levelCounter^2 + 1          
                    if Decoded(c) == (a - 1) && Cecoded(c) == (b - 1)                    
                        probabilityMatrix(a,b) = probabilityMatrix(a,b) + 1;
                    end
                end
            end
        end

end