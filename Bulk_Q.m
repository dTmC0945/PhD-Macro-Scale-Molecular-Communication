% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- Bulk Flow Analysis by Daniel T. McGuiness ----------------------------
% -------------------------------------------------------------------------

clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")

% Noise Parameters --------------------------------------------------------

mu_a    = 1.21*10^-3;              % Mean
sigma_a = sqrt(0.0960*10^-6);   % Standard deviation

% Data Reading ------------------------------------------------------------

A = xlsread('Bulk_Data_Q.xlsx');

% Memory Allocation -------------------------------------------------------

correlation                     = zeros(1, 19);
experimentalSignalEnergy        = zeros(1, 19);
experimentalSignal              = zeros(19, 450);
experimentalVariance            = zeros(1, 19);
experimentalMaximumAmplitude    = zeros(1, 19);
discreteTheoreticalSNR          = zeros(1, 19);
theoreticalAmplitude            = zeros(1, 19);
theoreticalSignalEnergy         = zeros(1, 19);
discreteTheoreticalSignalEnergy = zeros(1, 19);

signalLegendText1               = cell(1, 9);
signalLegendText2               = cell(1, 9);

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

xlabel('Time [s]','Interpreter','Latex')
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

 
x = 2.5;
diffusion = 0.14;
d_l = 0.001;
u = 0;
w = 0;
r = 0.3;
x_e = x;
x_d = x;
m = 2.4;

for i = 1:1:300
    
    v = i/1000; 
    
    Em = m*expint(v*x/diffusion);
    
    message = [0 0 1  0]*Em;
    k  =      [0 11 28 11];

    Tt      = 10; 

    Cm = transmission(diffusion, v, x, message, k, Tt);
       
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
CmA(i,:) = Cm + WGN;     

end

x_axis_theo = 0.1:0.1:30;
x_axis_exp  = 3:1:21;

for i = 1:1:19
    theoreticalAmplitude(i) = max(CmA(20 + 10*i, :));
end

for i = 1:1:300
    theoreticalSignalEnergy(i) = sum(abs(transpose(CmA(i, :)).^2));
end

%corr_amp = corrcoef( theo_data_amp,   max_amp_exp);             
max_amp_theo = max(transpose(CmA));

%plot(x_axis_exp, max_amp_exp,   'o',  'Linewidth'         ,    1.2,...
%                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
%                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(x_axis_theo, max_amp_theo,'Linewidth',2)

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

xlim([1 29])
xticks([1 5 9 13 17 21 25 29])
xticklabels({'0','1000','2000', '3000','4000','5000','6000','7000'})

% Signal Energy Analysis ==================================================

figure

plot(x_axis_exp, experimentalSignalEnergy,        'o',  'Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(x_axis_theo, theoreticalSignalEnergy,'Linewidth',2)

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

xlabel('Carrier Flow [ml/min]','Interpreter','Latex')
ylabel('Signal Energy [nW]','Interpreter','Latex')

xlim([1 29])
xticks([1 5 9 13 17 21 25 29])
xticklabels({'0','1000','2000', '3000','4000','5000','6000','7000'})

set(legend('Experimental Results', 'Theoretical Model'),...
           'Interpreter','Latex' , ...
           'Location'   ,'Northeast')

figure

experimentalSNR = 10*log10(experimentalSignalEnergy/sigma_a^2);
theoreticalSNR  = 10*log10(theoreticalSignalEnergy /sigma_a^2);

for i = 1:1:19
    discreteTheoreticalSignalEnergy(i)  = theoreticalSignalEnergy(20 + 10*i);
    discreteTheoreticalSNR(i)           = theoreticalSNR(20 + 10*i);
end

correlationEnergy   = findCorrelation(experimentalSignalEnergy, ...
                                      discreteTheoreticalSignalEnergy);
correlationSNR      = findCorrelation(experimentalSNR, ...
                                      discreteTheoreticalSNR);

plot(x_axis_exp, experimentalSNR, 'o','Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(x_axis_theo, theoreticalSNR,'Linewidth',2)

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

figure

for i = 1:1:19
    correlation(i) = findCorrelation(experimentalSignal(1,1:450), ...
                                     experimentalSignal(i,1:450));
end

figure  

plot(correlation,'o','Linewidth',1.2,'MarkerFaceColor',[0    0.4470    0.7410],'MarkerEdgeColor',[0    0.4470    0.7410])
hold on  

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
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

ylabel('$\sigma_Q^2$','Interpreter','Latex')
xlabel('Carrier Flow [ml/min]','Interpreter','Latex')

xlim([500 5000])

figure

plot(correlation,'o','Linewidth',1.2,'MarkerFaceColor',[0    0.4470    0.7410],'MarkerEdgeColor',[0    0.4470    0.7410])
hold on

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
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

ylabel('$\rho(A, B)$','Interpreter','Latex')
xlabel('Correlated pairs ($Q_{500},Q_n$)','Interpreter','Latex')

% Theory Analysis =========================================================

figure

for i = 1:1:5

    plot(CmA(20 + 10*i, :),        'LineWidth' ,    2  , ...
                                   'LineStyle' ,   '--')
    hold on
    plot(experimentalSignal(i, :), 'LineWidth' ,    2);
    hold on

end
         


%corrcoef(CmA(30,1:480),P0500(1:480))
%corrcoef(CmA(40,1:480),P0750(1:480))
%corrcoef(CmA(50,1:480),P1250(1:480))
%corrcoef(CmA(60,1:480),P1500(1:480))
%corrcoef(CmA(70,1:480),P1750(1:480))
%corrcoef(CmA(80,1:480),P2000(1:480))
%corrcoef(CmA(90,1:480),P2250(1:480))
%corrcoef(CmA(100,1:480),P2500(1:480))
%corrcoef(CmA(110,1:480),P2750(1:480))
%corrcoef(CmA(120,1:480),P3000(1:480))
%corrcoef(CmA(130,1:480),P3250(1:480))
%corrcoef(CmA(140,1:460),P3500(1:460))
%corrcoef(CmA(150,1:460),P3750(1:460))
%corrcoef(CmA(160,1:460),P4000(1:460))
%corrcoef(CmA(170,1:440),P4250(1:440))
%corrcoef(CmA(180,1:460),P4500(1:460))
%corrcoef(CmA(190,1:460),P4750(1:460))
%corrcoef(CmA(200,1:460),P5000(1:460))

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
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

P = legend('500 ml/min (T)'   ,'500 ml/min (E)' , ...
           '750 ml/min (T)'   ,'750 ml/min (E)' , ...
           '1000 ml/min (T)'  ,'1000 ml/min (E)' , ...
           '1250 ml/min (T)'  ,'1250 ml/min (E)' , ...
           '1500 ml/min (T)'  ,'1500 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')

figure

for i = 6:1:11

    plot(CmA(20 + 10*i, :),        'LineWidth' ,    2  , ...
                                   'LineStyle' ,   '--')
    hold on
    plot(experimentalSignal(i, :), 'LineWidth' ,    2);
    hold on

end

hold off
                       
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
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

P = legend('1750 ml/min (T)'  ,'1750 ml/min (E)' , ...
           '2000 ml/min (T)'  ,'2000 ml/min (E)' , ...
           '2250 ml/min (T)'  ,'2250 ml/min (E)' , ...
           '2500 ml/min (T)'  ,'2500 ml/min (E)' , ...
           '2750 ml/min (T)'  ,'2750 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')

xlim([0 450])

figure

for i = 12:1:16

    plot(CmA(20 + 10*i, :),        'LineWidth' ,    2  , ...
                                   'LineStyle' ,   '--')
    hold on
    plot(experimentalSignal(i, :), 'LineWidth' ,    2);
    hold on

end

hold off
            
            
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
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')
ylim([0 0.15])

P = legend('3000 ml/min (T)'  ,'3000 ml/min (E)' , ...
           '3250 ml/min (T)'  ,'3250 ml/min (E)' , ...
           '3500 ml/min (T)'  ,'3500 ml/min (E)' , ...
           '3750 ml/min (T)'  ,'3750 ml/min (E)' , ...
           '4000 ml/min (T)'  ,'4000 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')

figure

% Color definitions -------------------------------------------------------

CT_M = cbrewer('qual', 'Paired', 10);

% -------------------------------------------------------------------------

for i = 17:1:19

    plot(CmA(20 + 10*i, :),        'LineWidth' ,    2  , ...
                                   'LineStyle' ,   '--')
    hold on
    plot(experimentalSignal(i, :), 'LineWidth' ,    2);
    hold on

end

hold off
            
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
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

P = legend('4250 ml/min (T)'  ,'4250 ml/min (E)' , ...
           '4500 ml/min (T)'  ,'4500 ml/min (E)' , ...
           '4750 ml/min (T)'  ,'4750 ml/min (E)' , ...
           '5000 ml/min (T)'  ,'5000 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')


