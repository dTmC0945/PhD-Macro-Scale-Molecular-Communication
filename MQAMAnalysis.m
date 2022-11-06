% ========================= BEGIN CODE ====================================

clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")

% Variables ---------------------------------------------------------------

x           = 2.5;                      % Transmission Distance (cm)
v           = 0.1;                      % Diffusivity coefficient (cm^2/s)
M           = 1;                        % Injected mass (kg)
D_R         = 1e-5;
bitLength   = 100;                      % Transmitted Random bit length
symbolDuration          = 20;                       % Bit duration (counter) (sec)
r           = 0.3;
levelIterator   = 1;
levelCounter    = 2;

% loop parameters ---------------------------------------------------------

l_min  = 25;                     	% maximum value
l_max  = 25;                        % mininum value
l_mid  = 1;                        % incrimental value

% -------------------------------------------------------------------------
    
S_I_XY          = zeros(l_max./l_mid+1 ,3);   % Error
T_I_XY          = zeros(l_max./l_mid+1 ,3);   % Error
S_Error         = zeros(l_max./l_mid+1 ,3);   % Error    
T_Error         = zeros(l_max./l_mid+1 ,3);   % Error
T_Correct       = zeros(l_max./l_mid+1 ,3);   % Error
S_Correct       = zeros(l_max./l_mid+1 ,3);   % Error

% -------------------------------------------------------------------------

for  D_L = 0.1:0.1:0.4

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
        Rn              = zeros(1, bitLength);

        S_Adj_Probability_Matrix = zeros(levelCounter^2,levelCounter^2);
        T_Adj_Probability_Matrix = zeros(levelCounter^2,levelCounter^2);
        % Cleainng the input matrix and the probability matrix before each
        % iteration

        Probability_Matrix  = zeros(levelCounter^2 + 1,levelCounter^2 + 1);
        Input_Matrix        = zeros(000000000000000 + 1,levelCounter^2 + 1);

        randomBits1 = 2*(randperm(length(...
                         randi(levelCounter,1,bitLength))));
        randomBits2 = 2*(randperm(length(...
                         randi(levelCounter,1,bitLength))));


        % Delegating the random bits to two distinct matrices -------------

        [bitSequence1, ...
         bitRepetition1]     = splitter(randomBits1);

        [bitSequence2, ...
         bitRepetition2]     = splitter(randomBits2);

        CmP = transmission2D(D_R, D_L, v, x, r,...
                             bitSequence1, bitRepetition1,...
                             symbolDuration);
        
        CmQ = transmission2D(D_R, D_L, v, x, r,...
                             bitSequence2, bitRepetition2,...
                             symbolDuration);

        signal2D = complex(CmP - (levelCounter + 1),...
                           CmQ - (levelCounter + 1));
        
        normal = sqrt(sum(abs(signal2D .^2)) ...
                   / (symbolDuration*bitLength));

        noise = 10^(-Es_N0/20)/sqrt(2)*(randn(1,bitLength*symbolDuration) + 1i*randn(1,bitLength*symbolDuration)); % white guassian noise, 0dB variance
        
        %signalNoised = addNoise(signal2D, 0, 1./noise);

        noisedSignal = signal2D/normal + noise; % additive white gaussian noise

        K = snr(noisedSignal,noise);
        
        T_Prob = decodingTheory2D(levelCounter, Es_N0);

        [inputMatrix, probabilityMatrix] = decodingSimulation2D(levelCounter, bitLength, symbolDuration, noisedSignal, normal);

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

        S_Error(counter,levelIterator) = sum(sum(S_Error_Matrix));
        T_Error(counter,levelIterator) = sum(sum(T_Error_Matrix));  

        S_Correct(counter,levelIterator) = sum(sum(eye(levelCounter^2)...
                                .* S_Adj_Probability_Matrix));
        T_Correct(counter,levelIterator) = sum(sum(eye(levelCounter^2)...
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

        S_I_XY(counter,levelIterator) = S_Adj_H_Y - S_Adj_H_YX;
        T_I_XY(counter,levelIterator) = T_Adj_H_Y - T_Adj_H_YX;

       counter = counter + 1;

    end
    
     Kr(levelIterator,:) = real(Rn);
     Ki(levelIterator,:) = imag(Rn);   
     
    levelIterator = levelIterator + 1;
    
   
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