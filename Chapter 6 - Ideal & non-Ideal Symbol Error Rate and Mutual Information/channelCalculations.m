function [chValue, IValue] = channelCalculations(inputMatrix, ...
                                               Decoded, ...
                                               levelCounter, ...
                                               bitLength, ...
                                               Adj_Probability_Matrix, ...
                                               counter, ...
                                               bitCounter, ...
                                               dB)

        for c = 1:1:length(Decoded)
            
            for a = 1:1:levelCounter
                
                for b = 1:1:levelCounter
                    
                    Enput(a,b) =  inputMatrix(a)...
                        /bitLength*Adj_Probability_Matrix(b,a);
                             
                end
            end
        end
        
        Enput_check = sum(sum(Enput));      % P_Y Check the sum 
                                            % all probabilities equal to 1
                                            
        P_Y_adj = sum(Enput);
        
        for a = 1:1:levelCounter
            
                H_Y(a) = - P_Y_adj(a)*mylog2(P_Y_adj(a));
                
        end
        
        Adj_Input = inputMatrix./bitLength;
        
        Adj_H_Y = sum(H_Y);
        
        for a = 1:1:levelCounter
            for b = 1:1:levelCounter
                
                H_YX(a,b) = Enput(a,b)*mylog2((Adj_Input(1,a))/Enput(a,b));
    
            end
        end
        
        H_YX(isnan(H_YX)) = 0;
        
        Adj_H_YX = sum(sum(H_YX));

        IValue = Adj_H_Y - Adj_H_YX;

        chValue = 1/2*log2(1 + 10^(dB/10));

end