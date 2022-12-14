% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- Experimental and Theoretical BER Analysis by Daniel T. McGuiness -----
% -------------------------------------------------------------------------

clear variables
clc

% Variables ---------------------------------------------------------------

distance        = 2.5;                      % Transmission Distance (cm)
advectiveFlow   = 0.1;                      % Advective Flow        (cm/s)
diffusion       = 0.124;                    % Diffusivity           (cm^2/s)
mass            = 1;                        % Injected Mass         (ng)

experimentalBitLength   = 200;              % Experimentally sent bit count
bitLength               = 10000;            % Transmitted Random bit length
bitCounter              = 1;                % Bits per symbol
levelCounter            = 2^(bitCounter);   % Number of levels in the channel

startTime   = [18 18 16 20 20 20];          % start time for experimental measurements
duration    = [5  10 15 20 25 30];          % symbol duration of experiments

% Generating the filenames ------------------------------------------------

filename = {'Ber5secExperimentalData.xlsx', 'Ber10secExperimentalData.xlsx', ...
            'Ber15secExperimentalData.xlsx','Ber20secExperimentalData.xlsx', ...
            'Ber25secExperimentalData.xlsx','Ber30secExperimentalData.xlsx'};

% Memory Allocation -------------------------------------------------------

experimentalData    = zeros(1, 6);  Error       = zeros(1, 6);
I_XY                = zeros(1, 6);  Channel     = zeros(1, 6);
corrTheoretical     = zeros(1, 6);

% Acquiring the BER experimental data -------------------------------------

for i = 1:1:6
experimentalData(i) = experimentalBERDataAnalysis(filename{i}, ...
                                                  duration(i), ...
                                                  experimentalBitLength, ...
                                                  startTime(i))/experimentalBitLength;
end

% Generating the theoretical Values ---------------------------------------

counter = 1;     theoreticalIterator = 1;

theoryDuration = 1:1:30;

for i = 1:1:30             % Counter operator (dB)

    symbolDuration = theoryDuration(i);
    
    dB  = 25;
    
    SNR = 10^(dB/10);
   
    % Memory Allocation ---------------------------------------------------

        H_YX            = zeros(levelCounter,levelCounter); % Contional(YX)
        H_Y             = zeros(1,levelCounter);            % Entropy  (Y)
        Enput           = zeros(levelCounter,levelCounter); % Probability 

    % Cleainng the input matrix and the probability matrix before each
    % iteration
    
    Probability_Matrix  = zeros(levelCounter,levelCounter);
    Input_Matrix        = zeros(1,levelCounter);

    % Generating random bits ----------------------------------------------
        
    randomBits          = randi([0 levelCounter-1],1,bitLength);  
            
    % Delegating the random bits to two distinct matrices -----------------

    [bitSequence, ...
     bitRepetition]     = splitter(randomBits);

    % Simulating physical transmission ------------------------------------

    Cm                  = transmission(diffusion, advectiveFlow, distance, ...
                                       bitSequence, ...
                                       bitRepetition, symbolDuration);

    % Adding noise to the signal ------------------------------------------
    
    mu_a = 0; sigma_a   = sqrt(mean(Cm.^2)/SNR);

    CmNoised            = addNoise(Cm, mu_a, sigma_a);
                                             
    % Generating discrete received values from the transmission -----------
    
        [Decoded, Cecoded] = decoder(bitLength, ...
                                     levelCounter, ...
                                     randomBits, ...
                                     CmNoised, ...
                                     symbolDuration);     
    
    % Calculation of the error values per transmission --------------------

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
    
   % Channel calculation function -----------------------------------------
       
        [chValue, IValue] = channelCalculations(inputMatrix, ...
                                              Decoded, ...
                                              levelCounter, ...
                                              bitLength, ...
                                              Adj_Probability_Matrix, ...
                                              counter, ...
                                              bitCounter, ...
                                              dB);

    % Assigning the variables to matrices----------------------------------

    Error(counter)   = EValue;      % Error Value
   
    if mod(i,5) == 0
        corrTheoretical(theoreticalIterator) = Error(i);
                        theoreticalIterator  = theoreticalIterator + 1;
    end

    I_XY(counter)    = IValue;      % Mutual Information
    Channel(counter) = chValue;     % Shannon Channel

    % ---------------------------------------------------------------------
        
    counter = counter + 1;

end

figure % plotting SER for theory vs. experiment ---------------------------

semilogy(theoryDuration, Error, 'Linewidth', 2)
hold on
semilogy(duration, experimentalData,'s','Linewidth',2,...
                'MarkerFaceColor',[0.8500    0.3250    0.0980]   , ...
                'MarkerEdgeColor',[0.8500    0.3250    0.0980]   , ...
                'MarkerSize', 8)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YScale'      , 'log'     , ...
  'Box'         , 'on'      , ...
  'LineWidth'   , 1.5         );

% Label information -------------------------------------------------------

ylim([1e-4,1e-0]);

xlabel('Symbol Duration [s]','Interpreter','Latex')
ylabel('Symbol Error Rate [SER]','Interpreter','Latex')

set(legend('Theoretical Model', 'Experimental Results'), ...
           'Interpreter', 'Latex', ...
           'Orientation', 'Vertical', ...
           'Location'   , 'Southwest')


figure    % ---------------------------------------------------------------

plot(theoryDuration, I_XY, 'Linewidth', 2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'      , ...
  'LineWidth'   , 1.5         );

% Label information -------------------------------------------------------

xlabel('Symbol Duration [s]'            , 'Interpreter', 'Latex')
ylabel('Mutual Information [bit/sym]'   , 'Interpreter', 'Latex')

% Correlation -------------------------------------------------------------

correlation = findCorrelation(corrTheoretical, experimentalData);

% ======================== END OF CODE ====================================