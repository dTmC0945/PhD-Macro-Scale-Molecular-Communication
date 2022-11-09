% MOLECULAR SINGLE VELOCITY ANALYSIS -----------------------------------------
% Daniel McGuiness 06.07.2018 ---------------------------------------------
% -------------------------------------------------------------------------

clear all
clc

% Variables ---------------------------------------------------------------

x = 2.5;                                    % Transmission Distance(cm)
v = 0.1;                                    % Advective Flow       (cm/s)
r = 0.3;                                    % Detector Radius      (cm)

D_L = 0.124;                                  % Longitudional coefficient of 
                                            % Diffusivity          (cm^2/s)

M = 1;                                      % Injected Mass        (ng)
D_R = 0.001;                                % Radial coefficient of 
                                            % Diffusivity          (cm^2/s)

%x = 2.5;                                % Transmission Distance (cm)
D = 0.124;                              % Advective Flow (cm)
%M = 0.9;                                % Injected mass (kg)
v = 0.135;
bitLength       = 1000;                    % Transmitted Random bit length
bitCounter      = 1;
levelCounter    = 2^(bitCounter);    % Number of levels in the channel

A1 = experimentalBERDataAnalysis('Ber5secExperimentalData.xlsx',  5, 200, 18);
A2 = experimentalBERDataAnalysis('Ber10secExperimentalData.xlsx', 10, 200, 18);
A3 = experimentalBERDataAnalysis('Ber15secExperimentalData.xlsx', 15, 200, 16);
A4 = experimentalBERDataAnalysis('Ber20secExperimentalData.xlsx', 20, 200, 20);
A5 = experimentalBERDataAnalysis('Ber25secExperimentalData.xlsx', 25, 200, 20);
A6 = experimentalBERDataAnalysis('Ber30secExperimentalData.xlsx', 30, 200, 20);

CBER_Error  = [26 22 28 29; 20 17 11 10 ; 12 1 4 1; 11 5 4 0; 2 3 2 0]./100;


    Error               = zeros(1, 5);
    I_XY                = zeros(1, 5);
    Channel             = zeros(1, 5);


counter = 1;

for symbolDuration = [5 10 15 20 25 30]             % Counter operator (dB)
    
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

    Cm                  = transmission(D, v, x, ...
                                       bitSequence, ...
                                       bitRepetition, symbolDuration);

    % Adding noise to the signal ------------------------------------------
    
    mu_a = 0; sigma_a   = sqrt(mean(Cm.^2)/SNR);

    CmNoised            = addNoise(Cm, mu_a, sigma_a);
         
    % ---------------------------------------------------------------------
                                    
    t = 1:1:length(CmNoised);    % Transmission time x-axis value

    % Generating discrete received values from the transmission
    
    % ---------------------------------------------------------------------

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
    I_XY(counter)    = IValue;      % Mutual Information
    Channel(counter) = chValue;     % Shannon Channel

    % ---------------------------------------------------------------------
        
    counter = counter + 1;

end


y_pos = [3.33333333333334,3.66666666666667,1,0.666666666666667,0]./100;
y_neg = [4.66666666666666,2.33333333333333,1,0.333333333333333,0]./100;

%BERE =[38.6666666666667,19.3333333333333,6,4.33333333333333,0.333333333333333]./100;
BERE =[A1,A2,A3,A4,A5, A6]./200;

dB_axis = [5 10 15 20 25 30];
exp_axis = [5 10 15 20 25 30];


hold on
plot(exp_axis, BERE,'o','Linewidth',1.2)
% hold on
semilogy(dB_axis,Error,'s','Linewidth',2,...
                'MarkerFaceColor',[0.8500    0.3250    0.0980]   , ...
                'MarkerEdgeColor',[0.8500    0.3250    0.0980]   , ...
                'MarkerSize', 8)
            
% set(he1                        , ...
%   'LineWidth'       , 1.5      , ...
%   'Color'           , [0    0.4470    0.7410]   , ...
%   'MarkerSize'      , 8        , ...
%   'MarkerEdgeColor' , [0    0.4470    0.7410]  , ...
%   'MarkerFaceColor' , [0    0.4470    0.7410]   );

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

ylabel('Symbol Error Rate [SER]','Interpreter','Latex')
xlabel('Symbol Duration [s]','Interpreter','Latex')
%title({'$D$ = 0.124 $\mathrm{cm^2/s}$, $u_x$ = 0.136 $\mathrm{cm/s}$'; '$x$ = 2.5 $\mathrm{cm}$, $\rho$ = 0.94'},'Interpreter','Latex','Fontsize',16)

label = legend('Experimental Results','Theoretical Model');
set(label,'Interpreter','Latex','Orientation','Vertical','Location','Southwest')


figure

plot(dB_axis,I_XY,'-s','Linewidth',2,...
                'MarkerFaceColor',[0    0.4470    0.7410]   , ...
                'MarkerEdgeColor',[0    0.4470    0.7410])
hold on


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

ylabel('Mutual Information [bit/sym]','Interpreter','Latex')
xlabel('Symbol Duration [s]','Interpreter','Latex')

BER_Error = [Error(1) Error(2) Error(3) Error(4) Error(5) Error(6)];

corrcoef(BER_Error,BERE)