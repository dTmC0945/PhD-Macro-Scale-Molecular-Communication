function Cm = transmission( D, v, x, message, k, Tt)

% This function calculates the transmitted signal using theory.
% Inputs are;
% D: Coefficient of diffusivity of the chemical used in transmission (m2/s)
% v: The advective flow present in the medium (m/s)
% x: Distance of transmission (m)
% k: the symbol which repeats until there is a change (it is an array)
% message: the transmitted bits in alternating form (it is an array)
% Tt: the symbol length.

    Ts = 0;     Cm = 0; 

    for j = 2:1:length(message)
            
            if message(j) > message(j-1)
                
                M_old = Cm(end);
                
                for t = 1:1:k(j)*Tt
                    
                    Q(t) = M_old + (message(j) - message(j-1)) ...
                            - ((message(j) - message(j-1))*...
                            (erf((x-(t)*v)/(2*D*sqrt(1/(D*(t)))*(t))) + ...
                            erf((x+(t)*v)/(2*D*sqrt(1/(D*(t)))*(t)))))/2 ;
                        
                    Cm(t + Ts) = Q(t); 
    
                    
                end
                
            else
                    M_r = Cm(end); 
                
                for t = 1:1:k(j)*Tt
                    
                    Q(t) = (abs(M_r - message(j)))...
                        *(erf((x-(t)*v)/(2*D*sqrt(1/(D*(t)))*(t))) ...
                        + erf((x+(t)*v)/(2*D*sqrt(1/(D*(t)))*(t))))/2 ...
                        + message(j);
                        
                    Cm(t + Ts) = Q(t); 
                    
                end                
            end
            
            Ts = k(j)*Tt + Ts;
        
    end
end