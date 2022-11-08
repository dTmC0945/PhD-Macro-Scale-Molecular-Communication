% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- Symbol Error Rate and Mutual Information by Daniel T. McGuiness ------
% -------------------------------------------------------------------------

clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")

% Variables ---------------------------------------------------------------

x = 0; v = inf; D = 1e-50; symbolDuration= 1;      

min = -10;  mid = 2;    max = 40;  

bitMax = 5; bitLength = 10000; 

% Memory Allocation -------------------------------------------------------

Error   = zeros((max-min)/2 + 1, bitMax);
Channel = zeros((max-min)/2 + 1, bitMax);
I_XY    = zeros((max-min)/2 + 1, bitMax);

% -------------------------------------------------------------------------

for bitCounter = 1:1:bitMax
    
    levelCounter = 2^(bitCounter);    % Number of levels in the channel
        
    counter = 1;
    
    for dB = min:mid:max                 % Counter operator (dB)

        SNR = 10^(dB/10);
            
        % Memory Allocation -----------------------------------------------

        H_YX            = zeros(levelCounter,levelCounter); % Contional(YX)
        H_Y             = zeros(1,levelCounter);            % Entropy  (Y)
        Enput           = zeros(levelCounter,levelCounter); % Probability 
    
        % Generating random bits ------------------------------------------
        
        randomBits          = randi([0 levelCounter-1],1,bitLength);  
                
        % Delegating the random bits to two distinct matrices -------------

        [bitSequence, ...
         bitRepetition]     = splitter(randomBits);
    
        % Simulating physical transmission --------------------------------

        Cm                  = transmission(D, v, x, ...
                                           bitSequence, ...
                                           bitRepetition, symbolDuration);
    
        % Adding noise to the signal --------------------------------------
        
        mu_a = 0; sigma_a   = sqrt(mean(Cm.^2)/SNR);

        CmNoised            = addNoise(Cm, mu_a, sigma_a);
             
        % -----------------------------------------------------------------

        [Decoded, Cecoded] = decoder(bitLength, ...
                                     levelCounter, ...
                                     randomBits, ...
                                     CmNoised, ...
                                     symbolDuration);     

        % Calculation of the error values per transmission ----------------

        [EValue, ...
         inputMatrix, ...
         Adj_Probability_Matrix, ...
         probabilityCheck, ...
         inputCheck]    = errorCalculator(Decoded, ...
                                          Cecoded, ...
                                          levelCounter, ...
                                          bitLength, ...
                                          counter, ...
                                          bitCounter);
        
        % Channel calculation function ------------------------------------
       
        [chValue, IValue] = channelCalculations(inputMatrix, ...
                                              Decoded, ...
                                              levelCounter, ...
                                              bitLength, ...
                                              Adj_Probability_Matrix, ...
                                              counter, ...
                                              bitCounter, ...
                                              dB);

        % Assigning the variables to matrices------------------------------
 
        Error(counter,bitCounter)   = EValue;   % Error Value
        I_XY(counter,bitCounter)    = IValue;   % Mutual Information
        Channel(counter,bitCounter) = chValue;  % Shannon Channel

        % -----------------------------------------------------------------
        
        counter = counter + 1;                  % counter iterator
        
    end
end

% -------------------------------------------------------------------------

plotDataChannelSER(Error, I_XY, Channel, bitCounter, min, mid, max);

function plotDataChannelSER(Error,I_XY, Channel, bitCounter, min, mid, max)

% -------------------------------------------------------------------------
% First Plot (Symbol Error Rate) ------------------------------------------
% -------------------------------------------------------------------------

dB_axis = (min:mid:max);                        % dB values for the x-axis 

semilogy(dB_axis,Error,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'      , ...
  'LineWidth'   , 1.5         );

% Label information -------------------------------------------------------

ylabel('Symbol Error Rate [SER]','Interpreter','Latex')
xlabel('SNR [dB]','Interpreter','Latex')
ylim([1e-4,1e-0]);

% Legend information ------------------------------------------------------

modLegend = cell(1,bitCounter + 1);

for i = 1:1:bitCounter
    modLegend{i} = join(compose("M = %d", 2^ i));

    if i == bitCounter
        modLegend{i+1} = "AWGN Channel";
    end
end

set(legend(modLegend{1:bitCounter}),...
    'Interpreter','Latex','Location','Southwest')

% -------------------------------------------------------------------------
% Second Plot (Mutual Information) ----------------------------------------
% -------------------------------------------------------------------------

figure

for i = 1:1:bitCounter
    plot(dB_axis,I_XY(:,i),'Linewidth',2)
    hold on
end

plot(dB_axis,Channel,'--k','Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'     , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'     , ...
  'LineWidth'   , 1.5         );

% Label information -------------------------------------------------------

ylabel('Mutual Information [bit/sym]','Interpreter','Latex')
xlabel('SNR [dB]','Interpreter','Latex')

set(legend(modLegend),...
    'Interpreter','Latex','Location','Northwest')

end

% ======================== END OF CODE ====================================