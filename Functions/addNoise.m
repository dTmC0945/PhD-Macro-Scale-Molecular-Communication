function singalNoised = addNoise(Cm, mu_a, sigma_a)

WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive 
    
singalNoised = Cm + WGN; 

end