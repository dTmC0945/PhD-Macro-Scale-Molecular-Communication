function [Output1, Output2] = decoder(bit_length, level_counter,randomBits, Cm1,Tt)


    Decoded         = zeros(1,bit_length);             % Decoded
    Cecoded         = zeros(1,bit_length);
    Input_Alphabet  = zeros(1,level_counter);
    Alphabet        = zeros(1,level_counter);

    % Generating discrete received values from the transmission
    
    for i = 2:1:bit_length

               R1(i) = Cm1(Tt*i);   % message1 received values

    end
    
    % Generation of the symbol values for use in communications
    
    for a = 1:1:level_counter % Number of symbol values
        
        Alphabet(a) = a - 1; 
        
    end
    
    % Naming of the symbol for each constellation values
    
    alpha_counter = 1;
    
    for a = 1:1:level_counter % In-phase value
            
            Input_Alphabet(a) = Alphabet(alpha_counter);
            
            alpha_counter = alpha_counter + 1;
    end

    
    % Decoded Values
    
    for i = 2:1:bit_length
        
        for a = 1:1:level_counter
                
                if R1(i) < Input_Alphabet(1,1) - 0.5 
                    
                    Decoded(i) = 0;
                    
                elseif R1(i) > Input_Alphabet(1,level_counter) + 0.5
                    
                    Decoded(i) = level_counter - 1;
                    
                elseif  R1(i) >= (a - 1 - 1/2) ...
                     && R1(i) <= (a - 1 + 1/2)
                    Decoded(i) = Input_Alphabet(1,a);
                    
                end
                
        end

        
        for a = 1:1:level_counter

                if  randomBits(i) >= (a  - 1 ...
                        - 1/2) ...
                 && randomBits(i) <= (a  - 1 ...
                        + 1/2) 
                    
                        Cecoded(i) = Input_Alphabet(1,a);
                end
                
       
        end
        
    end

    Output1 = Decoded;
    Output2 = Cecoded;

end