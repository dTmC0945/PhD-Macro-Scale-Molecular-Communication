function [message1,k1] = splitter(bit_transmission)

    Re1 = [];
    Pe1 = [];
    
    i = 0;
    
    w = 1;
    
    for n = 1:1:length(bit_transmission)

        if n == 1

            i = i + 1;

            Re1(w) = bit_transmission(n);

            Pe1(w) = i;

        elseif bit_transmission(n-1) == bit_transmission(n)

            i = i + 1;

            Pe1(w) = i;

            Re1(w) = bit_transmission(n);

        elseif bit_transmission(n-1) ~= bit_transmission(n)

            Pe1(w) = i;

            w = w + 1;

            Re1(w) = bit_transmission(n);

            i = 1;

            if n == length(bit_transmission)

                Pe1(w) = 1;

                w = w + 1;

             end
        end
    end

    sum_P1 = sum(Pe1);      % Sum check (must be equal to bit_length)

    message1 = [0 Re1];     % reassign message1 (must add 0)

    k1  = [0 Pe1];          % Reassign bit frequency for message1

end