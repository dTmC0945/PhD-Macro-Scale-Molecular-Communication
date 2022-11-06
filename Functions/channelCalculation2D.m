function [outputArg1,outputArg2] = channelCalculation2D(inputArg1,inputArg2)

     

        % Generation the amount of input a single output received.

        Probability_length = sum(probabilityMatrix,1);
        T_Probability_length = ones(1,modLevelCounter^2)*bitLength/modLevelCounter^2;

        S_H_YX          = zeros(modLevelCounter^2 ,modLevelCounter^2);   % Contional
        T_H_YX          = zeros(modLevelCounter^2 ,modLevelCounter^2);   % Contional
        S_H_Y           = zeros(1, modLevelCounter^2);                      % Ent. (Y)
        T_H_Y           = zeros(1, modLevelCounter^2);                 % Ent. (Y)
        S_Enput         = zeros(modLevelCounter^2 ,modLevelCounter^2);   % Probability
        T_Enput         = zeros(modLevelCounter^2 ,modLevelCounter^2);   % Probability 
        S_Error_Matrix  = zeros(modLevelCounter^2 ,modLevelCounter^2);   % Error
        T_Error_Matrix  = zeros(modLevelCounter^2 ,modLevelCounter^2);   % Error
                              % Cecoded

        S_Adj_Probability_Matrix = zeros(modLevelCounter^2,modLevelCounter^2);
        T_Adj_Probability_Matrix = zeros(modLevelCounter^2,modLevelCounter^2);

    for a = 1:1:modLevelCounter^2          
        for b = 1:1:modLevelCounter^2           
            S_Adj_Probability_Matrix(a,b) = ...
                probabilityMatrix(a,b)./Probability_length(b);           
           T_Adj_Probability_Matrix(a,b) = ...
           Theo_Probability_Matrix(a,b)./T_Probability_length(b);
        end
    end

    % Error Count ---------------------------------------------------------

    S_Error_Input_Matrix = (S_Adj_Probability_Matrix - eye(modLevelCounter^2)...
                            .*S_Adj_Probability_Matrix);
    T_Error_Input_Matrix = (T_Adj_Probability_Matrix - eye(modLevelCounter^2)...
                            .* T_Adj_Probability_Matrix);

    for a = 1:1:modLevelCounter^2           
        for b = 1:1:modLevelCounter^2            
            S_Error_Matrix(a,b) = ...
                S_Error_Input_Matrix(b,a)*inputMatrix(a)/bitLength;
            T_Error_Matrix(a,b) = ...
              T_Error_Input_Matrix(b,a)*T_Probability_length(a)/bitLength;            
        end
    end

    S_Error(counter,varIterator) = sum(sum(S_Error_Matrix));
    T_Error(counter,varIterator) = sum(sum(T_Error_Matrix));  

    S_Correct(counter,varIterator) = sum(sum(eye(modLevelCounter^2)...
                            .* S_Adj_Probability_Matrix));
    T_Correct(counter,varIterator) = sum(sum(eye(modLevelCounter^2)...
                            .* T_Adj_Probability_Matrix));
    % ---------------------------------------------------------------------

    for c = 1:1:bitLength       
        for a = 1:1:modLevelCounter^2            
            for b = 1:1:modLevelCounter^2
                T_Enput(a,b) = T_Probability_length(a)/bitLength*T_Adj_Probability_Matrix(b,a);
                S_Enput(a,b) =  inputMatrix(a)/bitLength*S_Adj_Probability_Matrix(b,a);                         
            end
        end
    end

    S_Enput_check = sum(sum(S_Enput));      % P_Y Check the sum 
    T_Enput_check = sum(sum(T_Enput));      % P_Y Check the sum 

    S_P_Y_adj = sum(S_Enput);
    T_P_Y_adj = sum(T_Enput);

    for a = 1:1:modLevelCounter^2       
            S_H_Y(a) = - S_P_Y_adj(a)*mylog2(S_P_Y_adj(a));
            T_H_Y(a) = - T_P_Y_adj(a)*mylog2(T_P_Y_adj(a));         
    end

    S_Adj_Input = inputMatrix./bitLength;
    T_Adj_Input = T_Probability_length./bitLength;

    S_Adj_H_Y = sum(S_H_Y);
    T_Adj_H_Y = sum(T_H_Y);

    for a = 1:1:modLevelCounter^2
        for b = 1:1:modLevelCounter^2         
           S_H_YX(a,b) = S_Enput(a,b)*mylog2((S_Adj_Input(1,a))/S_Enput(a,b));
           T_H_YX(a,b) = T_Enput(a,b)*mylog2((T_Adj_Input(1,a))/T_Enput(a,b));
        end
    end

    S_H_YX(isnan(S_H_YX)) = 0;
    T_H_YX(isnan(T_H_YX)) = 0;

    S_Adj_H_YX = sum(sum(S_H_YX));
    T_Adj_H_YX = sum(sum(T_H_YX));

    S_I_XY(counter,varIterator) = S_Adj_H_Y - S_Adj_H_YX;
    T_I_XY(counter,varIterator) = T_Adj_H_Y - T_Adj_H_YX;

end