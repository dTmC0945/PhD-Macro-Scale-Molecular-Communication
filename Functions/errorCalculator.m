function [Error, ...
          Input_Matrix, ...
          Adj_Probability_Matrix, ...
          Probability_check, ...
          Input_check] = errorCalculator(Decoded, ...
                                                    Cecoded, ...
                                                    level_counter, ...
                                                    bit_length, ...
                                                    counter, ...
                                                    bit_counter)
    
    Input_Matrix            = zeros(1,level_counter);
    Probability_Matrix      = zeros(level_counter,level_counter);
    Adj_Probability_Matrix  = zeros(level_counter,level_counter);
    Error_Matrix            = zeros(level_counter,level_counter); % Error

    for c = 1:1:bit_length
            
        for a = 1:1:level_counter + 1
                
            if Cecoded(c) == (a - 1)
                    
                Input_Matrix(a) = Input_Matrix(a) + 1;
            end
                
            for b = 1:1:level_counter + 1
                 
                if Decoded(c) == (a - 1) && Cecoded(c) == (b - 1)
                        
                    Probability_Matrix(a,b) = Probability_Matrix(a,b) + 1;
                 end
             end
        end
    end
        
    % Check values for the probability matrix and the input matrix
    % Both of them must equal to 1
    
    Probability_check = sum(sum(Probability_Matrix))./bit_length;

    Input_check       = sum(sum(Input_Matrix))      ./bit_length;
    
    % Generation the amount of input a single output received.
    
    Probability_length = sum(Probability_Matrix,1);

    % Adjusting the probability matrix by dividing the probability matrix
    % with the probability length
    
    for a = 1:1:level_counter
            
        for b = 1:1:level_counter
            
            Adj_Probability_Matrix(a,b) = ...
                Probability_Matrix(a,b)./Probability_length(b);
            
        end
    end
        
    % Error Count ---------------------------------------------------------
    
    Error_Input_Matrix = (Adj_Probability_Matrix - eye(level_counter)...
                        .*Adj_Probability_Matrix);
    
    for a = 1:1:level_counter
            
        for b = 1:1:level_counter
            
            Error_Matrix(a,b) = ...
                Error_Input_Matrix(b,a)*Input_Matrix(a)/bit_length ;
            
        end
    end
               
    Error = sum(sum(Error_Matrix));

end