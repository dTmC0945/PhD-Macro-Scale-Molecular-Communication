% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- Open Distance Transmission by Daniel T. McGuiness --------------------
% -------------------------------------------------------------------------

clear variables
clc

experimentalData = xlsread('Master_Open_Distance.xlsx');   % Excel Data

environmentalNoise = xlsread('Noise.xlsx');              % Noise Data


% Maximum Amplitude Values from experimental measurements -----------------

maxamp25    = range([9.9243e-10,...         % 2.5 cm 
                     1.0950e-09,...
                     1.1000e-09])*10^9;
maxamp50    = range([8.2725e-10,...         % 5.0 cm 
                     8.8774e-10,...
                     8.0549e-10])*10^9;
maxamp75    = range([5.3547e-10,...         % 7.5 cm 
                     4.2969e-10,...
                     5.1658e-10])*10^9;
maxamp100   = range([5.1797e-11,...         % 10.0 cm 
                     7.8666e-11])*10^9;
maxamp125   = range([8.1053e-12,...         % 12.5 cm 
                     2.9265e-11])*10^9;
maxamp150   = range([6.6418e-12,...         % 15.0 cm 
                     7.3241e-12,...
                     2.7424e-12])*10^9;

% Errorbars for maximum Amplitude -----------------------------------------                 
  

errBar = [maxamp25 maxamp50 maxamp75 maxamp100 maxamp125 maxamp150]./2;

% -------------------------------------------------------------------------

% Amplitude vs. Distance Figure

for i = 1:1:6 % Amplitude values (max)
    xAxis(i) = 2.5*i;
    Amp(i) = max(experimentalData(:,i));
    E(i) = sum(abs(experimentalData(1:180,i)).^2);
end
     


% Signal Energy -----------------------------------------------------------

N1 = sum(abs(environmentalNoise(1:180)./10^3).^2);

errBarEnergy = [0.7952, 0.9876, 1.4246, 0.1599, 1.982e-05, 1.8698e-04]./2; 

% Noise Probability Analysis ----------------------------------------------

mu          = 1.21815;
sigma       = sqrt(0.0960);

[fCDF,xCDF] = ecdf(environmentalNoise);
    
gaussianCDF = 1/2.*(1 + erf((xCDF - mu)./(sigma.*sqrt(2))));

ks          = max(abs(gaussianCDF - fCDF));

[h,p]       = kstest((fCDF-mu)/sigma);

figure

check = chi2gof(environmentalNoise);     % needs to be 0

plot(xCDF, fCDF,        'LineWidth', 2);    hold on;
plot(xCDF, gaussianCDF, 'LineWidth', 2);    hold on;

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   ,...
  'Fontsize'            , 15        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'XGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       );

% -------------------------------------------------------------------------

xlim([0 2])

xlabel('$x$','Interpreter','Latex','Fontsize',15)
ylabel('CDF','Interpreter','Latex','Fontsize',15)

set(legend('Experimental CDF','Gaussian Distribution CDF'), ...
           'Interpreter', 'Latex', ...
           'Location'   , 'Northwest');


% SIGNAL MODELLING --------------------------------------------------------
% 11.07.2018 --------------------------------------------------------------


v = 0.18;                                       % Transmission Distance (cm)
d = 0.15;                                       % Advective Flow (cm)
M = 1.2;                                        % Injected mass (gr)

x_length = 0:0.01:20;

cc = 1;

for i = 0:0.01:20                          % Decay generator
    
    K(cc) = 7.5e-05*i^(2.6);
    %b = 2.616
    cc = cc + 1;
end

counter = 1;
    
for x = 0:0.01:20                           % Open Distance (cm)
        
    Tt      = 60;                           % Bit duration (counter) (sec)
    
    Cm      = 0;                            % First value of transmission
                                            % (keep it 0)

    M_old   = 0;                            % Leftover chemicals in the system  
                                            % (keep it 0)
    Ts      = 0;                            % Transmission time                 
                                            % (keep it 0)
    Q_neg   = 0;                            % Flushed chemicals                 
                                            % (keep it 0)
    
    message = [0 0 1 0]*M;                  % reassign message1 (must add 0)

    kf  = [0 1 1 1];                        % Reassign bit frequency for message1

    C_time = zeros(1,15);
    
    for j = 2:1:length(message)

            if message(j) > message(j-1)

                M_old = Cm(end);

                for t = 1:1:kf(j)*Tt
                    
                    Q(t) = M_old +  exp(-K(counter)*t)*(message(j) -  ...
                                 ((message(j))*...
                        (erf((x - (t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                        1))/2) ;

                    Cm(t + Ts) = Q(t); 

                end

            else
                    M_r = Cm(end); 

                for t = 1:1:kf(j)*Tt

                    Q(t) = (abs(M_r - message(j)))*...
                        exp(-K(counter)*t)*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                         1)/2;

                    Cm(t + Ts) = Q(t); 

                end

            end

            Ts = kf(j)*Tt + Ts;

    end
    
    C_delay = [C_time Cm];
    
    %   Noise parameters
    
    mu_a = 1.21*10^-3;              % mean
    
    sigma_a = sqrt(0.0960*10^-6);   % standard deviation
    
    WGN = transpose(sigma_a.*randn(length(C_delay),1) + mu_a);    % Additive
    
    C_final = C_delay + WGN;     

    Signal(counter,:) = C_final;
    
    Signal_Amplitude(counter) = max(C_final);
    
    counter = counter + 1;
end

figure
subplot(3,2,1)
plot(experimentalData(:,1),'Linewidth',1.2)
hold on
plot(Signal(250,:),'--','Linewidth',1.2)
hold on


% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YTick'       , 0:0.3:1.2   , ...
  'XTick'       , 0:60:180   , ...
  'LineWidth'   , 1.5         );

% -------------------------------------------------------------------------
xlim([0 180])
ylim([0 1.2])
title('(a) 2.5 cm','Interpreter','Latex')

S_label = legend('Experimental Results','Theoretical Model');
set(S_label,'Interpreter','Latex','Orientation','Vertical')


subplot(3,2,2)
plot(experimentalData(:,2),'Linewidth',1.2)
hold on
plot(Signal(500,:),'--','Linewidth',1.2)
hold on


% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
    'YAxisLocation'   ,'right',...
  'YTick'       , 0:0.25:1, ...
  'XTick'       , 0:60:180  , ...
  'LineWidth'   , 1.5         );

% -------------------------------------------------------------------------
xlim([0 180])
ylim([0 1])
title('(b) 5 cm','Interpreter','Latex')

subplot(3,2,3)
plot(experimentalData(:,3),'Linewidth',1.2)
hold on
plot(Signal(750,:),'--','Linewidth',1.2)
hold on


% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YTick'       , 0:0.15:0.6, ...
  'XTick'       , 0:60:180  , ...
  'LineWidth'   , 1.5         );

% -------------------------------------------------------------------------
xlim([0 180])
ylim([0 0.6])
title('(c) 7.5 cm','Interpreter','Latex')


subplot(3,2,4)
plot(experimentalData(:,4),'Linewidth',1.2)
hold on
plot(Signal(1000,:),'--','Linewidth',1.2)
hold on


% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YTick'       , 0:0.04:0.16, ...
    'YAxisLocation'   ,'right',...
  'XTick'       , 0:60:180  , ...
  'LineWidth'   , 1.5         );

% -------------------------------------------------------------------------
xlim([0 180])
ylim([0 0.16])
title('(d) 10 cm','Interpreter','Latex')

subplot(3,2,5)
plot(experimentalData(:,5),'Linewidth',1.2)
hold on
plot(Signal(1250,:),'--','Linewidth',1.2)
hold on

title('(e) 12.5 cm','Interpreter','Latex')

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YTick'       , 0:0.01:0.03, ...
  'XTick'       , 0:60:180  , ...
  'LineWidth'   , 1.5         );

% -------------------------------------------------------------------------
xlim([0 180])
ylim([0 0.03])

subplot(3,2,6)
plot(experimentalData(:,6),'Linewidth',1.2)
hold on
plot(Signal(1500,:),'--','Linewidth',1.2)
hold on

title('(f) 15 cm','Interpreter','Latex')

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'        , 15        , ...
  'XGrid'           , 'on'      , ...
  'YGrid'           , 'on'      , ...
  'GridLineStyle'   ,'--'       , ...
  'TickDir'         , 'in'      , ...
  'TickLength'      , [.02 .02] , ...
  'XMinorTick'      , 'on'      , ...
  'YMinorTick'      , 'on'      , ...
  'YAxisLocation'   ,'right',...
  'YTick'           , 0:0.0015:0.006, ...
  'XTick'           , 0:60:180  , ...
  'LineWidth'       , 1.5         );

% -------------------------------------------------------------------------
xlim([0 180])
ylim([0 0.006])


% Amplitude Modelling -----------------------------------------------------

Noise_gen = transpose(sigma_a.*randn(length(Signal),1) + mu_a);

figure

errorbar(xAxis,Amp,errBar,'o','Linewidth',1.2)
hold on
plot(x_length,Signal_Amplitude,'-','Linewidth',2)
hold on
plot(x_length,Noise_gen,'-','Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'Yscale'      ,'log'      , ...
  'GridLineStyle','--'      , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1.5         );
set(gca, 'YTick', [10.^-4 10.^-3 10.^-2 10.^-1 10.^0 10.^1])
% -------------------------------------------------------------------------

xlabel('Transmission Distance [cm]','Interpreter','Latex')
ylabel('Signal Amplitude [nA]','Interpreter','Latex')

E_label = legend('Experimental Results','Theoretical Model','Environment Noise');
set(E_label,'Interpreter','Latex')
xlim([0 20])
ylim([1e-4,1e1])

Amp_model_data = [max(Signal(250,:))  max(Signal(500,:))  max(Signal(750,:)) ...
                  max(Signal(1000,:)) max(Signal(1250,:)) max(Signal(1500,:))];

rho_amp = corrcoef(Amp_model_data,Amp);
              



Energy_model   = sum(abs(Signal).^2,2);
Spectral_noise = sum(abs(Noise_gen(1:195)).^2)*ones(1,length(x_length));


figure

errorbar(xAxis,E,errBarEnergy,'o','Linewidth',1.2)
hold on
plot(x_length,Energy_model,'-','Linewidth',2)
hold on
plot(x_length,Spectral_noise,'--k','Linewidth',2)


Energy_model_data = [Energy_model(250) Energy_model(500) Energy_model(750) ...
                     Energy_model(1000) Energy_model(1250) Energy_model(1500)];

rho_energy = corrcoef(Energy_model_data,E);                 
                 
% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'Yscale'      ,'log'      , ...
  'GridLineStyle','--'      , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1.5         );
set(gca, 'YTick', [10.^-4 10.^-3 10.^-2 10.^-1 10.^0 10.^1 10.^2])

% -------------------------------------------------------------------------
xlabel('Transmission Distance [cm]','Interpreter','Latex')
ylabel('Signal Energy [nW]','Interpreter','Latex')

E_label = legend('Experimental Results','Theoretical Model','$N_{0}$');
set(E_label,'Interpreter','Latex')
xlim([0 20])

s1 = corrcoef(Signal(250,1:180),experimentalData(1:180,1));
s2 = corrcoef(Signal(500,1:180),experimentalData(1:180,2));
s3 = corrcoef(Signal(750,1:180),experimentalData(1:180,3));
s4 = corrcoef(Signal(1000,1:180),experimentalData(1:180,4));
s5 = corrcoef(Signal(1250,1:180),experimentalData(1:180,5));
s6 = corrcoef(Signal(1500,1:180),experimentalData(1:180,6));

% SNR Test ------------------------------------------------------------

snr_n = 10*log10(E/N1);

snr_m = 10*log10(Energy_model/Spectral_noise(1));

figure
plot(xAxis,snr_n,'o','Linewidth',1.2)
hold on
plot(x_length,snr_m,'-','Linewidth',2)

% Nice plot code ----------------------------------------------------------
set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------

xlabel('Transmission Distance [cm]','Interpreter','Latex')
ylabel('Signal Energy [nW]','Interpreter','Latex')

K_values = [K(250) K(500) K(750) K(1000) K(1250) K(1500)];

SNR_label = legend('Experimental Results','Theoretical Model');
set(SNR_label,'Interpreter','Latex')