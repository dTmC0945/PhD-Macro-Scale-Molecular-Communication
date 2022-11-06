% ========================= BEGIN CODE ====================================

clear variables
clc

% Variables ---------------------------------------------------------------

x           = 2.5;                      % Transmission Distance (cm)
v           = 0.1;                      % Diffusivity coefficient (cm^2/s)
M           = 1;                        % Injected mass (kg)
D_R         = 1e-5;
bitLength   = 100;                      % Transmitted Random bit length
Tt          = 20;                       % Bit duration (counter) (sec)
r           = 0.3;

% loop parameters ---------------------------------------------------------

l_min  = 25;                     	% maximum value
l_max  = 25;                        % mininum value
l_mid  = 1;                        % incrimental value

% -------------------------------------------------------------------------

level_iterator = 1;
levelCounter = 2;
    
        S_I_XY          = zeros(l_max./l_mid+1 ,3);   % Error
        T_I_XY          = zeros(l_max./l_mid+1 ,3);   % Error
        S_Error         = zeros(l_max./l_mid+1 ,3);   % Error    
        T_Error         = zeros(l_max./l_mid+1 ,3);   % Error
    	T_Correct       = zeros(l_max./l_mid+1 ,3);   % Error
       	S_Correct       = zeros(l_max./l_mid+1 ,3);   % Error

for  D_L = 0.1:0.1:0.4
                                  % Advective Flow (cm)
    counter = 1;

    for Es_N0 = l_min:l_mid:l_max                                       % Counter operator (dB)

        Alphabet        = zeros(levelCounter   ,levelCounter);       % Alpabet    
        S_H_YX          = zeros(levelCounter^2 ,levelCounter^2);   % Contional
        T_H_YX          = zeros(levelCounter^2 ,levelCounter^2);   % Contional
        S_H_Y           = zeros(1              ,levelCounter^2);                 % Ent. (Y)
        T_H_Y           = zeros(1              ,levelCounter^2);                 % Ent. (Y)
        S_Enput         = zeros(levelCounter^2 ,levelCounter^2);   % Probability
        T_Enput         = zeros(levelCounter^2 ,levelCounter^2);   % Probability 
        S_Error_Matrix  = zeros(levelCounter^2 ,levelCounter^2);   % Error
        T_Error_Matrix  = zeros(levelCounter^2 ,levelCounter^2);   % Error
        Input_Alphabet  = zeros(levelCounter   ,levelCounter);       % Input Alph.
        Decoded         = zeros(1               ,bitLength);                               % Decoded
        Cecoded         = zeros(1               ,bitLength);                               % Cecoded
        mu_m            = zeros(levelCounter^2 ,2);
        Rn              = zeros(1, bitLength);
        T_Prob          = zeros(levelCounter^2, levelCounter^2);

        S_Adj_Probability_Matrix = zeros(levelCounter^2,levelCounter^2);
        T_Adj_Probability_Matrix = zeros(levelCounter^2,levelCounter^2);
        % Cleainng the input matrix and the probability matrix before each
        % iteration

        Probability_Matrix  = zeros(levelCounter^2 + 1,levelCounter^2 + 1);
        Input_Matrix        = zeros(000000000000000 + 1,levelCounter^2 + 1);


        Cm1 = zeros(1, 1); Cm2 = zeros(1, 1); % First value of transmission

        Ts      = 0;                        % Transmission time                 
                                            % (keep it 0)
        Q_neg   = 0;                        % Flushed chemicals                 
                                            % (keep it 0)
        Q1 = zeros(1, 1); Q2 = zeros(1, 1);
        
        randomBits1 = 2*(randperm(length(...
                         randi(levelCounter,1,bitLength))));
        randomBits2 = 2*(randperm(length(...
                         randi(levelCounter,1,bitLength))));


        % Delegating the random bits to two distinct matrices -------------

        [bitSequence1, ...
         bitRepetition1]     = splitter(randomBits1);

        [bitSequence2, ...
         bitRepetition2]     = splitter(randomBits2);

        for con = 2:1:length(message1)

                if message1(con) > message1(con-1)

                    M_old = Cm1(end);

                    for t = 1:1:k1(con)*Tt

                        Q1(t) = M_old + (message1(con) - message1(con-1)) ...
                                 - (message1(con) - message1(con-1))/2*...
                                 (erf((x-v*t)/sqrt(4*D_L*t)) + ...
                                  erf((x+v*t)/sqrt(4*D_L*t)));
                   
                        Cm1(t + Ts) = Q1(t); 

                    end

                else
                        M_r = Cm1(end); 

                    for t = 1:1:k1(con)*Tt

                        Q1(t) = (abs(M_r - message1(con)))/2 * ...
                           (erf((x-v*t)/sqrt(4*D_L*t)) + ...
                            erf((x+v*t)/sqrt(4*D_L*t))) + ...
                            message1(con);

                        Cm1(t + Ts) = Q1(t); 

                    end

                end

                Ts = k1(con)*Tt + Ts;

        end

        Ts = 0; Mr = 0; M_old = 0; con = 0;

        for con = 2:1:length(message2)

                if message2(con) > message2(con-1)

                    M_old = Cm2(end);

                    for t = 1:1:k2(con)*Tt

                        Q2(t) = M_old + (message2(con) - message2(con-1)) ...
                              - (message2(con) - message2(con-1))/2 * ...
                                (erf((x-v*t)/sqrt(4*D_L*t)) + ...
                                 erf((x+v*t)/sqrt(4*D_L*t)));

                        Cm2(t + Ts) = Q2(t); 

                    end

                else
                        M_r = Cm2(end); 

                    for t = 1:1:k2(con)*Tt

                        Q2(t) = (abs(M_r - message2(con)))/2 * ...
                        (erf((x-v*t)/sqrt(4*D_L*t)) + ...
                         erf((x+v*t)/sqrt(4*D_L*t))) + ...
                         message2(con);

                        Cm2(t + Ts) = Q2(t); 

                    end

                end

                Ts = k2(con)*Tt + Ts;

        end    

        S = Cm1 + Cm2*1i - (levelCounter + 1)*(1 + 1i);
        normal = sqrt(sum(abs(S .^2)) / (Tt*bitLength));

        noise = 10^(-Es_N0/20)/sqrt(2)*(randn(1,length(Cm1)) + 1i*randn(1,length(Cm1))); % white guassian noise, 0dB variance
        R = S/normal + noise; % additive white gaussian noise

        K = snr(R,noise);
        int_count = 1;

        theo_normalization = 1.219877889386052;
        
        for a = 1:1:levelCounter
            for b = 1:1:levelCounter
                mu_m(int_count,:) = ([2*b 2*a]- (levelCounter + 1))/theo_normalization;
                int_count = int_count + 1;
            end
        end
   
        min = (-levelCounter + 1)/theo_normalization;  mu_x = mu_m(:,1);
        max = (+levelCounter - 1)/theo_normalization;  mu_y = mu_m(:,2);

        cov_matrix = [1 0;0 1]*10^(-Es_N0/10)/2;

        for a = 1:1:levelCounter^2

            for b = 1:1:levelCounter^2  
                
                if mu_x(b) > min ...                    % Inner Cases
                && mu_x(b) < max ...
                && mu_y(b) < max ...
                && mu_y(b) > min  
                    
                    T_Prob(a,b) =  mvncdf([mu_x(b) - 1/theo_normalization , ...
                                           mu_y(b) - 1/theo_normalization], ...
                                          [mu_x(b) + 1/theo_normalization , ...
                                           mu_y(b) + 1/theo_normalization], ...
                                          [mu_m(a,1) mu_m(a,2)], ...
                                          cov_matrix);

                elseif  mu_x(b) == min ...
                     && mu_y(b) == min         % Corner Cases
                     T_Prob(a,b) = mvncdf([-10000, -10000], ...
                                          [min + 1/theo_normalization,min + 1/theo_normalization],...
                                          [mu_x(a) mu_y(a)],...
                                           cov_matrix);
                elseif  mu_x(b) == max ...
                     && mu_y(b) == min
                     T_Prob(a,b) = mvncdf([max-1/theo_normalization,...
                         -10000],...
                         [10000, ...
                         min+1/theo_normalization]...
                         ,[mu_x(a) mu_y(a)],...
                         cov_matrix);
                elseif  mu_x(b) == min ...
                     && mu_y(b) == max
                     T_Prob(a,b) = mvncdf([-10000, max-1/theo_normalization],[min+1/theo_normalization,  10000],[mu_x(a) mu_y(a)],cov_matrix); 
                elseif  mu_x(b) == max ...
                     && mu_y(b) == max
                     T_Prob(a,b) = mvncdf([max-1/theo_normalization, max-1/theo_normalization],[10000,   10000],[mu_x(a) mu_y(a)],cov_matrix);  
                elseif  mu_x(b) == min ...
                     && mu_y(b) <  max ...
                     && mu_y(b) >  min  % Edge Cases
                    T_Prob(a,b) = mvncdf([-10000, mu_y(b)-1/theo_normalization],[mu_x(b)+1/theo_normalization, mu_y(b)+1/theo_normalization],[mu_x(a) mu_y(a)],cov_matrix);
                elseif  mu_x(b) == max ...
                     && mu_y(b) <  max ...
                     && mu_y(b) >  min
                    T_Prob(a,b) = mvncdf([mu_x(b)-1/theo_normalization,  mu_y(b)-1/theo_normalization],[10000, mu_y(b)+1/theo_normalization],[mu_x(a) mu_y(a)],cov_matrix);         
                elseif  mu_y(b) == min ...
                     && mu_x(b) <  max ...
                     && mu_x(b) >  min
                    T_Prob(a,b) = mvncdf([mu_x(b)-1/theo_normalization,  -10000],[mu_x(b)+1/theo_normalization, mu_y(b)+1/theo_normalization],[mu_x(a) mu_y(a)],cov_matrix);
                elseif  mu_y(b) == max ...
                     && mu_x(b) <  max ...
                     && mu_x(b) >  min
                    T_Prob(a,b) = mvncdf([mu_x(b)-1/theo_normalization, mu_y(b)-1/theo_normalization],[mu_x(b)+1/theo_normalization, 10000],[mu_x(a) mu_y(a)],cov_matrix);
                end
            end
        end  

        % Generating discrete received values from the transmission

        for message_operator = 2:1:length(length_op)
                    Rn(message_operator) = R(Tt*message_operator)*normal + (levelCounter + 1)*(1 + 1i); 
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
                    Input_Matrix(a) = Input_Matrix(a) + 1; 
                end

                for b = 1:1:levelCounter^2 + 1          
                    if Decoded(c) == (a - 1) && Cecoded(c) == (b - 1)                    
                        Probability_Matrix(a,b) = Probability_Matrix(a,b) + 1;
                    end
                end
            end
        end

        % Check values for the probability matrix and the input matrix
        % Both of them must equal to 1

        Theo_Probability_Matrix = T_Prob.*bitLength./levelCounter^2;
        Probability_check = sum(sum(Probability_Matrix))./bitLength;
        Theo_check        = sum(sum(Theo_Probability_Matrix))./bitLength;
        Input_check       = sum(sum(Input_Matrix))      ./bitLength;

        % Generation the amount of input a single output received.

        Probability_length = sum(Probability_Matrix,1);
        T_Probability_length = ones(1,levelCounter^2)*bitLength/levelCounter^2;

        % Adjusting the probability matrix by dividing the probability matrix
        % with the probability length

        for a = 1:1:levelCounter^2          
            for b = 1:1:levelCounter^2           
                S_Adj_Probability_Matrix(a,b) = ...
                    Probability_Matrix(a,b)./Probability_length(b);           
               T_Adj_Probability_Matrix(a,b) = ...
               Theo_Probability_Matrix(a,b)./T_Probability_length(b);
            end
        end

        % Error Count ---------------------------------------------------------

        S_Error_Input_Matrix = (S_Adj_Probability_Matrix - eye(levelCounter^2)...
                                .*S_Adj_Probability_Matrix);
        T_Error_Input_Matrix = (T_Adj_Probability_Matrix - eye(levelCounter^2)...
                                .* T_Adj_Probability_Matrix);

        for a = 1:1:levelCounter^2           
            for b = 1:1:levelCounter^2            
                S_Error_Matrix(a,b) = ...
                    S_Error_Input_Matrix(b,a)*Input_Matrix(a)/bitLength;
                T_Error_Matrix(a,b) = ...
                  T_Error_Input_Matrix(b,a)*T_Probability_length(a)/bitLength;            
            end
        end

        S_Error(counter,level_iterator) = sum(sum(S_Error_Matrix));
        T_Error(counter,level_iterator) = sum(sum(T_Error_Matrix));  

        S_Correct(counter,level_iterator) = sum(sum(eye(levelCounter^2)...
                                .* S_Adj_Probability_Matrix));
        T_Correct(counter,level_iterator) = sum(sum(eye(levelCounter^2)...
                                .* T_Adj_Probability_Matrix));
        % ---------------------------------------------------------------------

        for c = 1:1:length(Decoded)       
            for a = 1:1:levelCounter^2            
                for b = 1:1:levelCounter^2
                    T_Enput(a,b) = T_Probability_length(a)/bitLength*T_Adj_Probability_Matrix(b,a);
                    S_Enput(a,b) =  Input_Matrix(a)/bitLength*S_Adj_Probability_Matrix(b,a);                         
                end
            end
        end

        S_Enput_check = sum(sum(S_Enput));      % P_Y Check the sum 
        T_Enput_check = sum(sum(T_Enput));      % P_Y Check the sum 

        S_P_Y_adj = sum(S_Enput);
        T_P_Y_adj = sum(T_Enput);

        for a = 1:1:levelCounter^2       
                S_H_Y(a) = - S_P_Y_adj(a)*mylog2(S_P_Y_adj(a));
                T_H_Y(a) = - T_P_Y_adj(a)*mylog2(T_P_Y_adj(a));         
        end

        S_Adj_Input = Input_Matrix./bitLength;
        T_Adj_Input = T_Probability_length./bitLength;

        S_Adj_H_Y = sum(S_H_Y);
        T_Adj_H_Y = sum(T_H_Y);

        for a = 1:1:levelCounter^2
            for b = 1:1:levelCounter^2         
               S_H_YX(a,b) = S_Enput(a,b)*mylog2((S_Adj_Input(1,a))/S_Enput(a,b));
               T_H_YX(a,b) = T_Enput(a,b)*mylog2((T_Adj_Input(1,a))/T_Enput(a,b));
            end
        end

        S_H_YX(isnan(S_H_YX)) = 0;
        T_H_YX(isnan(T_H_YX)) = 0;

        S_Adj_H_YX = sum(sum(S_H_YX));
        T_Adj_H_YX = sum(sum(T_H_YX));

        S_I_XY(counter,level_iterator) = S_Adj_H_Y - S_Adj_H_YX;
        T_I_XY(counter,level_iterator) = T_Adj_H_Y - T_Adj_H_YX;

       counter = counter + 1;

    end
    
     Kr(level_iterator,:) = real(Rn);
     Ki(level_iterator,:) = imag(Rn);   
     
    level_iterator = level_iterator + 1;
    
   
end

figure
subtightplot(3,3,[4 5 7 8],[0.01,0.01],0.15,0.15)


scatter(Kr(1,:),Ki(1,:),'filled','MarkerEdgeColor', [51 77 92]/255,'MarkerFaceColor', [51 77 92]/255)
hold on
scatter(Kr(2,:),Ki(2,:),'filled','MarkerEdgeColor', [69 178 157]/255,'MarkerFaceColor', [69 178 157]/255)
hold on
scatter(Kr(3,:),Ki(3,:),'filled','MarkerEdgeColor', [239 201 76]/255,'MarkerFaceColor', [239 201 76]/255)
hold on
scatter(Kr(4,:),Ki(4,:),'filled','MarkerEdgeColor', [226 122 63]/255,'MarkerFaceColor', [226 122 63]/255)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 14        , ...
  'TickDir'     , 'out'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XMinorGrid'  , 'on'      , ...
  'YMinorGrid'  , 'on'      , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.1         );

% -------------------------------------------------------------------------
xlabel('In-phase Chemical ($m_I$)','Interpreter','Latex')
ylabel('Out-phase Chemical ($m_Q$)','Interpreter','Latex')
xticks([1,2,3])
yticks([1,2,3])


MI_legend = legend('$D = 0.1\, \mathrm{cm^2/s}$',...
                   '$D = 0.2\, \mathrm{cm^2/s}$',...
                   '$D = 0.3\, \mathrm{cm^2/s}$',...
                   '$D = 0.4\, \mathrm{cm^2/s}$');
set(MI_legend,'Interpreter','Latex','box','off');
set(MI_legend,'position',[.71 .71 .1 .1]);
ylim([0.5 3.5])
xlim([0.5 3.5])




subtightplot(3,3,[1 2],[0.01,0.01],0.10,0.15)
nbins = 150;
h1 = histogram(Kr(1,:),nbins,'EdgeColor', 'none','FaceColor', [51 77 92]/255,'FaceAlpha',1)
hold on
h2 = histogram(Kr(2,:),nbins,'EdgeColor','none','FaceColor', [69 178 157]/255,'FaceAlpha',1)
hold on
h3 = histogram(Kr(3,:),nbins,'EdgeColor', 'none','FaceColor', [239 201 76]/255,'FaceAlpha',1)
hold on
h4 = histogram(Kr(4,:),nbins,'EdgeColor', 'none','FaceColor', [226 122 63]/255,'FaceAlpha',1)
%set(gca,'visible','off')
xlim([0.5 3.5])
% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 14        , ...
  'TickDir'     , 'out'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XMinorGrid'  , 'on'      , ...
  'YMinorGrid'  , 'on'      , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );
set(gca,'xticklabel',{[]})
%set(gca,'yticklabel',{[]})

%set(gca,'visible','off')
% -------------------------------------------------------------------------

subtightplot(3,3,[6 9],[0.01,0.01],0.15,0.10)
nbins = 150;
h5 = histogram(Ki(1,:),nbins,'EdgeColor', 'none','FaceColor', [51 77 92]/255,'FaceAlpha',1)
hold on
h6 = histogram(Ki(2,:),nbins,'EdgeColor','none','FaceColor', [69 178 157]/255,'FaceAlpha',1)
hold on
h7 = histogram(Ki(3,:),nbins,'EdgeColor', 'none','FaceColor', [239 201 76]/255,'FaceAlpha',1)
hold on
h8 = histogram(Ki(4,:),nbins,'EdgeColor', 'none','FaceColor', [226 122 63]/255,'FaceAlpha',1)
xlim([0.5 3.5])

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 14        , ...
  'TickDir'     , 'out'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YAxisLocation','right',...
  'XMinorGrid'  , 'on'      , ...
  'YMinorGrid'  , 'on'      , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );
set(gca,'CameraUpVector',[1,0,0]);
set(gca, 'XDir','reverse')
set(gca,'xticklabel',{[]})

%set(gca,'yticklabel',{[]})
%set(gca,'visible','off')
% ======================== END OF CODE ====================================