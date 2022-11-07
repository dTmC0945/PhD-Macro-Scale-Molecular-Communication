% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- Bulk Flow Analysis by Daniel T. McGuiness ----------------------------
% -------------------------------------------------------------------------

clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")
addpath("./Experimental Data")

% Parameters --------------------------------------------------------------
 
distance        = 2.5;
diffusion       = 0.14;
mass            = 2.4;
symbolDuration  = 10; 

mu_a            = 1.21*10^-3;              % Mean
sigma_a         = sqrt(0.0960*10^-6);   % Standard deviation

% Data Reading ------------------------------------------------------------

A = xlsread('BulkFlowExperimentalData.xlsx');

% Memory Allocation -------------------------------------------------------

correlation                         = zeros(1, 19);

experimentalSignalEnergy            = zeros(1, 19);
experimentalSignal                  = zeros(19, 450);
experimentalVariance                = zeros(1, 19);
experimentalMaximumAmplitude        = zeros(1, 19);

theoreticalSignal                   = zeros(1, 300);
theoreticalSignalEnergy             = zeros(1, 19);

discreteTheoreticalMaximumAmplitude = zeros(1,19);
discreteTheoreticalSignalEnergy     = zeros(1, 19);
discreteTheoreticalSNR              = zeros(1, 19);

CmNoised                            = zeros(300, 500);

signalLegendText1                   = cell(1, 9);
signalLegendText2                   = cell(1, 9);

% Adjustments for NaN data points ------------------------------------------

A(4, 297) = (A(4, 301) +  A(4, 293))/2;
A(4, 298) = (A(4, 302) +  A(4, 294))/2;
A(4, 299) = (A(4, 303) +  A(4, 295))/2;
A(4, 300) = (A(4, 304) +  A(4, 296))/2;
A(7, 174) = (A(7, 173) +  A(7, 175))/2;

% -------------------------------------------------------------------------

for i = 1:1:19
    experimentalSignal(i,:)         = A(i,1:450)* 10^9;
    experimentalMaximumAmplitude(i) = max(experimentalSignal(i,:));
    experimentalVariance(i)         = var(experimentalSignal(i,90:390));
    experimentalSignalEnergy(i)     = sum(abs(experimentalSignal(i,1:450)).^2);
end

% -------------------------------------------------------------------------

figure

for i = 1:2:19
    plot(experimentalSignal(i,:),'Linewidth',2);
    hold on;
end

xlim([0 450]);  ylim([0 1.2]);

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   , ...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'XTick'               , 0:60:480  , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlabel('Time [s]'           ,'Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

for i = 1:1:9
    signalLegendText1{i} = join(compose("%d ml/min", i*500));
end
       
set(legend(signalLegendText1),'Interpreter' ,'Latex', ...
                              'Location'    ,'Northwest')

figure

for i = 2:2:18
    plot(experimentalSignal(i,:), 'Linewidth', 2);
    hold on;
end

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   , ...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'XTick'               , 0:60:480  , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlabel('Time [s]'           ,'Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

for i = 1:1:9
    signalLegendText2{i} = join(compose("%d ml/min", i*500 + 250));
end

set(legend(signalLegendText2), ...
          'Interpreter' , 'Latex', ...
          'Location'    , 'Northwest')

xlim([0 450]); ylim([0 0.8]);

% -------------------------------------------------------------------------

figure

for i = 1:1:300
    
    advectiveFlow   = i/1000; 
    
    adjustedMass    = mass*expint(advectiveFlow*distance/diffusion);
    
    bitSequence     = [0 0  1  0]*adjustedMass;
    bitRepetition   = [0 11 28 11];


    Cm              = transmission(diffusion, advectiveFlow, distance, ...
                                   bitSequence, bitRepetition, ...
                                    symbolDuration);

    CmNoised(i,:)   = addNoise(Cm, mu_a, sigma_a);          

end

theoreticalXAxis    = 0.1:0.1:30;
experimentalXAxis   = 3:1:21;

for i = 1:1:19
    theoreticalSignal(i)        = max(CmNoised(20 + 10*i, :));
end

for i = 1:1:300
    theoreticalSignalEnergy(i)  = sum(abs(transpose(CmNoised(i, :)).^2));
end

theoreticalMaximumAmplitude     = max(transpose(CmNoised));

plot(experimentalXAxis, experimentalMaximumAmplitude,   'o', ...
                          'LineWidth'         , 1.2, ...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410], ...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(theoreticalXAxis, theoreticalMaximumAmplitude,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   , ...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'YScale'              , 'log'     , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlabel('Carrier Flow [ml/min]'  ,'Interpreter','Latex')
ylabel('Signal Current [nA]'    ,'Interpreter','Latex')

set(legend('Experimental Results', 'Theoretical Model'),...
           'Interpreter','Latex' , ...
           'Location'   ,'Northeast')

xlim([1 29]); xticks([1 5 9 13 17 21 25 29]);
xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000', '7000'});

% Signal Energy Analysis ==================================================

figure

plot(experimentalXAxis, experimentalSignalEnergy,  'o',  ...
                                            'LineWidth'         , 1.2, ...
                                            'MarkerFaceColor'   ,[0 0.4470 0.7410], ...
                                            'MarkerEdgeColor'   ,[0 0.4470 0.7410])
hold on
plot(theoreticalXAxis, theoreticalSignalEnergy,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'TickDir'     , 'out'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YScale'      , 'log'     , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Carrier Flow [ml/min]'  ,'Interpreter','Latex')
ylabel('Signal Energy [nW]'     ,'Interpreter','Latex')

xlim([1 29]); xticks([1 5 9 13 17 21 25 29]);
xticklabels({'0', '1000', '2000', '3000', '4000', '5000', '6000', '7000'});

set(legend('Experimental Results', 'Theoretical Model'), ...
           'Interpreter','Latex' , ...
           'Location'   ,'Northeast')

% SNR Calculation and Plot ------------------------------------------------

figure

experimentalSNR = 10*log10(experimentalSignalEnergy/sigma_a^2);
theoreticalSNR  = 10*log10(theoreticalSignalEnergy /sigma_a^2);

for i = 1:1:19
    discreteTheoreticalMaximumAmplitude(i)  = theoreticalMaximumAmplitude(20 + 10*i);
    discreteTheoreticalSignalEnergy(i)      = theoreticalSignalEnergy(20 + 10*i);
    discreteTheoreticalSNR(i)               = theoreticalSNR(20 + 10*i);
end

plot(experimentalXAxis, experimentalSNR, 'o','Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(theoreticalXAxis, theoreticalSNR,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   ,...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlabel('Carrier Flow [ml/min]'      ,'Interpreter','Latex')
ylabel('Signal-to-Noise ratio [dB]' ,'Interpreter','Latex')

xlim([1 29])
xticks([1 5 9 13 17 21 25 29])
xticklabels({'0','1000','2000', '3000','4000','5000','6000','7000'})

set(legend('Experimental Results', 'Theoretical Model'),...
           'Interpreter','Latex' , ...
           'Location'   ,'Northeast')

% Correlation of Experimental Signal Plot ---------------------------------

for i = 1:1:19
    correlation(i) = findCorrelation(experimentalSignal(1,1:450), ...
                                     experimentalSignal(i,1:450));
end

figure  

plot(correlation,'o', ...
                 'Linewidth', 1.2, ...
                 'MarkerFaceColor', [0 0.4470 0.7410], ...
                 'MarkerEdgeColor', [0 0.4470 0.7410])
hold on

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   ,...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

ylabel('$\rho(A, B)$','Interpreter','Latex')
xlabel('Correlated pairs ($Q_{500},Q_n$)','Interpreter','Latex')

% Signal shape comparison -------------------------------------------------

figure

for i = 1:1:5

    plot(CmNoised(20 + 10*i, :),    'LineWidth' ,    2  , ...
                                    'LineStyle' ,   '--')
    hold on
    plot(experimentalSignal(i, :),  'LineWidth' ,    2);
    hold on

end

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   ,...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

xlim([0 450])

P = legend('500 ml/min (T)'   ,'500 ml/min (E)' , ...
           '750 ml/min (T)'   ,'750 ml/min (E)' , ...
           '1000 ml/min (T)'  ,'1000 ml/min (E)' , ...
           '1250 ml/min (T)'  ,'1250 ml/min (E)' , ...
           '1500 ml/min (T)'  ,'1500 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')

figure

for i = 6:1:10

    plot(CmNoised(20 + 10*i, :),   'LineWidth' ,    2  , ...
                                   'LineStyle' ,   '--')
    hold on
    plot(experimentalSignal(i, :), 'LineWidth' ,    2);
    hold on

end

hold off
                       
% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   ,...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

P = legend('1750 ml/min (T)'  ,'1750 ml/min (E)' , ...
           '2000 ml/min (T)'  ,'2000 ml/min (E)' , ...
           '2250 ml/min (T)'  ,'2250 ml/min (E)' , ...
           '2500 ml/min (T)'  ,'2500 ml/min (E)' , ...
           '2750 ml/min (T)'  ,'2750 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')

xlim([0 450]); ylim([0 0.3]);

figure

for i = 11:1:15

    plot(CmNoised(20 + 10*i, :),        'LineWidth' ,    2  , ...
                                   'LineStyle' ,   '--')
    hold on
    plot(experimentalSignal(i, :), 'LineWidth' ,    2);
    hold on

end

hold off          
            
% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   ,...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')
xlim([0 450]);  ylim([0 0.15]);

P = legend('3000 ml/min (T)'  ,'3000 ml/min (E)' , ...
           '3250 ml/min (T)'  ,'3250 ml/min (E)' , ...
           '3500 ml/min (T)'  ,'3500 ml/min (E)' , ...
           '3750 ml/min (T)'  ,'3750 ml/min (E)' , ...
           '4000 ml/min (T)'  ,'4000 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')

figure

for i = 16:1:19

    plot(CmNoised(20 + 10*i, :),   'LineWidth' ,    2  , ...
                                   'LineStyle' ,   '--')
    hold on
    plot(experimentalSignal(i, :), 'LineWidth' ,    2);
    hold on

end

hold off
            
% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   ,...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

xlim([0 450]); ylim([0 0.05]);

P = legend('4250 ml/min (T)'  ,'4250 ml/min (E)' , ...
           '4500 ml/min (T)'  ,'4500 ml/min (E)' , ...
           '4750 ml/min (T)'  ,'4750 ml/min (E)' , ...
           '5000 ml/min (T)'  ,'5000 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')

% Correlation between theory and experiment -------------------------------

correlationAmplitude= findCorrelation(experimentalMaximumAmplitude, ...
                                      discreteTheoreticalMaximumAmplitude);
correlationEnergy   = findCorrelation(experimentalSignalEnergy, ...
                                      discreteTheoreticalSignalEnergy);
correlationSNR      = findCorrelation(experimentalSNR, ...
                                      discreteTheoreticalSNR);

% ======================== END OF CODE ====================================