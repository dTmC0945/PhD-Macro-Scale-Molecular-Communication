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

experimentalSignalEnergy        = zeros(1, 19);
experimentalSignal              = zeros(19, 481);
experimentalVariance            = zeros(1, 19);
experimentalMaximumAmplitude    = zeros(1, 19);
x_label                         = zeros(1, 19);

% -------------------------------------------------------------------------

for i = 1:1:19
    experimentalSignal(i,:)         = A(i,:)* 10^9;
    experimentalMaximumAmplitude(i) = max(experimentalSignal(i,:));
    experimentalVariance(i)         = var(experimentalSignal(i,90:390));
    experimentalSignalEnergy(i)     = sum(abs(experimentalSignal(i,1:450)).^2);
    x_label(i)                      = 500 + 250*(i-1);
end

% Adjustments for no data points ------------------------------------------

experimentalSignal(4, 297) = (experimentalSignal(4, 301) ...
                           +  experimentalSignal(4, 293))/2;
experimentalSignal(4, 298) = (experimentalSignal(4, 302) ...
                           +  experimentalSignal(4, 294))/2;
experimentalSignal(4, 299) = (experimentalSignal(4, 303) ...
                           +  experimentalSignal(4, 295))/2;
experimentalSignal(4, 300) = (experimentalSignal(4, 304) ...
                           +  experimentalSignal(4, 296))/2;
experimentalSignal(7, 173) = (experimentalSignal(7, 173) ...
                           +  experimentalSignal(7, 175))/2;

% -------------------------------------------------------------------------
                   
figure

for i = 1:2:19
    plot(experimentalSignal(i,:),'Linewidth',2);
    hold on;
end

xlim([0 480]);  ylim([0 1.2]);

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
    plot(experimentalSignal(i,:),'Linewidth',2);
    hold on;
end

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   , ...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XTick'       ,0:60:480   , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

for i = 1:1:9
    signalLegendText2{i} = join(compose("%d ml/min", i*500 + 250));
end

set(legend(signalLegendText2), ...
          'Interpreter' , 'Latex', ...
          'Location'    , 'Northwest')

xlim([0 480]); ylim([0 0.8]);

figure

 
for i = 1:1:19
    x(i) = 500 + 250*(i-1);

end
 
x = 2.5;
d = 0.14;
d_l = 0.001;
u = 0;
w = 0;
r = 0.3;
x_e = x;
x_d = x;
m = 2.4;



for i = 1:1:300
    
    v = i/1000; 
    
    Em = m*expint(v*x/d);
    
    message = [0 0 1  0]*Em;
    k  =      [0 11 28 11];

    Tt      = 10;  Ts      = 0;  Cm      = 0;

for j = 2:1:length(message)
        
        if message(j) > message(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + Em - Em/2 * (erf(r/sqrt(4*d_l*t)))^2*(erf((x-v*t)/sqrt(4*d*t)) + erf((x+v*t)/sqrt(4*d*t)));
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (M_r)/2 * (erf(r/sqrt(4*d_l*t)))^2*(erf((x-v*t)/sqrt(4*d*t)) + erf((x+v*t)/sqrt(4*d*t)));
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
CmA(i,:) = Cm + WGN;     

end

x_axis_theo = 0.1:0.1:30;
x_axis_exp  = 3:1:21;

for i = 1:1:19
    theoreticalAmplitude(i) = max(CmA(20 + 10*i, :));
    theoreticalSignalEnergy(i) = sum(abs(transpose(CmA(20 + 10*i, :)).^2));
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
xlim([1 29])
xticks([1 5 9 13 17 21 25 29])
xticklabels({'0','1000','2000', '3000','4000','5000','6000','7000'})

xlabel('Carrier Flow [ml/min]'  ,'Interpreter','Latex')
ylabel('Signal Current [nA]'    ,'Interpreter','Latex')

set(legend('Experimental Results', 'Theoretical Model'),...
           'Interpreter','Latex' , ...
           'Location'   ,'Northeast')

% Signal Energy Analysis ==================================================

figure

corr_eng = corrcoef(theo_data_energy, Energy);

plot(x_axis_exp, Energy,        'o',  'Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(x_axis_theo, Energy_Theory,'Linewidth',2)

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

SNR = 10*log10(Energy/sigma_a^2);

Smode = 10*log10(Energy_Theory/sigma_a^2);

SNR_data_theo = 10*log10(theo_data_energy/sigma_a^2);

corr_snr = corrcoef(SNR_data_theo, SNR);


plot(x_axis_exp, SNR, 'o','Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(x_axis_theo, Smode,'Linewidth',2)


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

xlabel('Carrier Flow [ml/min]'      ,'Interpreter','Latex')
ylabel('Signal-to-Noise ratio [dB]' ,'Interpreter','Latex')

xlim([1 29])
xticks([1 5 9 13 17 21 25 29])
xticklabels({'0','1000','2000', '3000','4000','5000','6000','7000'})

set(legend('Experimental Results', 'Theoretical Model'),...
           'Interpreter','Latex' , ...
           'Location'   ,'Northeast')

figure

for i = 2:1:19
    correlation = findCorrelation(experimentalSignal(1,1:458), ...
                                  experimentalSignal(i,1:458));
end

figure  

plot(x_label,Vr,'o','Linewidth',1.2,'MarkerFaceColor',[0    0.4470    0.7410],'MarkerEdgeColor',[0    0.4470    0.7410])
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

% Color definitions -------------------------------------------------------

% -------------------------------------------------------------------------

plot(CmA(30,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(1,:))
hold on
plot(P0500,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(2,:))
hold on
plot(CmA(40,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(3,:))
hold on
plot(P0750,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(4,:))
hold on
plot(CmA(50,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(5,:))
hold on
plot(P1000,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(6,:))
hold on
plot(CmA(60,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(7,:))
hold on
plot(P1250,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(8,:))
hold on
plot(CmA(70,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(9,:))
hold on
plot(P1500,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(10,:))
            

corrcoef(CmA(30,1:480),P0500(1:480))
corrcoef(CmA(40,1:480),P0750(1:480))
corrcoef(CmA(50,1:480),P1250(1:480))
corrcoef(CmA(60,1:480),P1500(1:480))
%corrcoef(CmA(70,1:480),P1750(1:480))
%corrcoef(CmA(80,1:480),P2000(1:480))
%corrcoef(CmA(90,1:480),P2250(1:480))
%corrcoef(CmA(100,1:480),P2500(1:480))
%corrcoef(CmA(110,1:480),P2750(1:480))
%corrcoef(CmA(120,1:480),P3000(1:480))
%corrcoef(CmA(130,1:480),P3250(1:480))
%corrcoef(CmA(140,1:460),P3500(1:460))
corrcoef(CmA(150,1:460),P3750(1:460))
corrcoef(CmA(160,1:460),P4000(1:460))
corrcoef(CmA(170,1:440),P4250(1:440))
corrcoef(CmA(180,1:460),P4500(1:460))
corrcoef(CmA(190,1:460),P4750(1:460))
corrcoef(CmA(200,1:460),P5000(1:460))
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

% Color definitions -------------------------------------------------------

CT_M = cbrewer('qual', 'Paired', 10);

% -------------------------------------------------------------------------

plot(CmA(80,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(1,:))
hold on
plot(P1750,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(2,:))
hold on
plot(CmA(90,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(3,:))
hold on
plot(P2000,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(4,:))
hold on
plot(CmA(100,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(5,:))
hold on
plot(P2250,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(6,:))
hold on
plot(CmA(110,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(7,:))
hold on
plot(P2500,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(8,:))
hold on
plot(CmA(120,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(9,:))
hold on
plot(P2750,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(10,:))
            
            
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

figure

% Color definitions -------------------------------------------------------

CT_M = cbrewer('qual', 'Paired', 10);

% -------------------------------------------------------------------------

plot(CmA(130,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(1,:))
hold on
plot(P3000,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(2,:))
hold on
plot(CmA(140,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(3,:))
hold on
plot(P3250,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(4,:))
hold on
plot(CmA(150,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(5,:))
hold on
plot(P3500,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(6,:))
hold on
plot(CmA(160,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(7,:))
hold on
plot(P3750,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(8,:))
hold on
plot(CmA(170,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(9,:))
hold on
plot(P4000,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(10,:))
            
            
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

plot(CmA(180,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(1,:))
hold on
plot(P4250,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(2,:))
hold on
plot(CmA(190,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(3,:))
hold on
plot(P4500,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(4,:))
hold on
plot(CmA(200,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(5,:))
hold on
plot(P4750,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(6,:))
hold on
plot(CmA(210,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(7,:))
hold on
plot(P5000,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(8,:))
            
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


