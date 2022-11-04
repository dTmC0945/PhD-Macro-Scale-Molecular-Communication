function CmNoised = addNoise(Cm, mu_a, sigma_a)

WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive 
    
CmNoised = Cm + WGN; 

end