clear all
clc

A = xlsread('Impulse_Analysis_Q.xlsx');

% Array Arrangement -------------------------------------------------------

p90     =   A(1,:)*10^9;
p120    =   A(3,:)*10^9;
p150    =   A(5,:)*10^9;
p180    =   A(7,:)*10^9;
p210    =   A(9,:)*10^9;
p240    =   A(11,:)*10^9;
p270    =   A(13,:)*10^9;
p300    =   A(15,:)*10^9;
p330    =   A(17,:)*10^9;
p360    =   A(19,:)*10^9;

% NaN Filtering -----------------------------------------------------------

p90(75)     = (p90(76)  + p90(74))/2;
p90(96)     = (p90(95)  + p90(97))/2;
p90(122)    = (p90(121) + p90(123))/2;
p90(125)    = (p90(124) + p90(126))/2;
p90(130)    = (p90(128)  + p90(132))/2;
p90(131)    = (p90(129)  + p90(133))/2;
p90(148)    = (p90(147)  + p90(149))/2;
p90(156)    = (p90(155)  + p90(157))/2;
p90(200)    = (p90(199)  + p90(201))/2;

p150(195)   = (p150(194) + p150(196))/2;

p180(360)   =  p180(359);

p210(47)    = (p210(46) + p210(48))/2;

p300(480)   =  p300(479);

p330(64)    = (p330(68) + p330(60))/2;
p330(65)    = (p330(69) + p330(61))/2;
p330(66)    = (p330(70) + p330(62))/2;
p330(67)    = (p330(71) + p330(63))/2;

%find(all(isnan(p360),1))

x = 90:30:360;      % Common x-axis for use in pulse width graphs

% Color definitions -------------------------------------------------------

CT = cbrewer('qual', 'Dark2', 10);

% -------------------------------------------------------------------------

plot(p90,   'Linewidth',    2,  'Color',    CT(1,:))
hold on
plot(p120,  'Linewidth',    2,  'Color',    CT(2,:))
hold on
plot(p150,  'Linewidth',    2,  'Color',    CT(3,:))
hold on
plot(p180,  'Linewidth',    2,  'Color',    CT(4,:))
hold on
plot(p210,  'Linewidth',    2,  'Color',    CT(5,:))
hold on
plot(p240,  'Linewidth',    2,  'Color',    CT(6,:))
hold on
plot(p270,  'Linewidth',    2,  'Color',    CT(7,:))
hold on
plot(p300,  'Linewidth',    2,  'Color',    CT(8,:))
hold on
plot(p330,  'Linewidth',    2,  'Color',    CT(9,:))
hold on
plot(p360,  'Linewidth',    2,  'Color',    CT(10,:))


% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'TickDir'     , 'out'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'XTick'       , 0:60:540  , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------
 
ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Time [s]','Interpreter','Latex')
L = legend( '90s','120s','150s','180s','210s', ...
            '240s','270s','300s','330s','360s');

set(L, 'Interpreter','Latex','Location','Northwest')

xlim([0 540])

figure   
   
x = 2.5;
v = 0.03;
d = 0.12;
d_l = 0.001;
u = 0;
w = 0;
r = 0.3;
x_e = x;
x_d = x;

% Noise Parameters ========================================================

mu_a = 1.21*10^-3;              % Mean
    
sigma_a = sqrt(0.0960*10^-6);   % Standard deviation

% =========================================================================

for i = 1:1:10
    
    Em = 0.36 +0.09*i;
    
    message = [0 0  1  0]*Em;
    k  =      [0 11 6 + 3*i - 2 40 - 3*i];

    Tt      = 10;  Ts      = 0;  Cm      = 0;

for j = 2:1:length(message)
        
        if message(j) > message(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + Em - Em/2 * ...
                    (erf(r/sqrt(4*d_l*t)))^2 * ...
                    (erf((x-v*t)/sqrt(4*d*t)) + erf((x+v*t)/sqrt(4*d*t)));
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (M_r)/2 * ...
                    (erf(r/sqrt(4*d_l*t)))^2 * ...
                    (erf((x-v*t)/sqrt(4*d*t)) + erf((x+v*t)/sqrt(4*d*t)));
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end

WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
CmA(i,:) = Cm + WGN;

end


% Color definitions -------------------------------------------------------

CT_M = cbrewer('qual', 'Paired', 10);

% -------------------------------------------------------------------------


plot(CmA(1,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(1,:))
hold on
plot(p90,   'Linewidth',    2,  'Color',    CT_M(2,:))
hold on
plot(CmA(3,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(3,:))
hold on
plot(p150,  'Linewidth',    2,  'Color',    CT_M(4,:))
hold on
plot(CmA(5,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(5,:))
hold on
plot(p210,  'Linewidth',    2,  'Color',    CT_M(6,:))
hold on
plot(CmA(7,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(7,:))
hold on
plot(p270,  'Linewidth',    2,  'Color',    CT_M(8,:))
hold on
plot(CmA(9,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(9,:))
hold on
plot(p330,  'Linewidth',    2,  'Color',    CT_M(10,:))
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

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Time [s]','Interpreter','Latex')

P = legend('90 s (T)'   ,'90 s (E)' , ...
           '150 s (T)'  ,'150 s (E)' , ...
           '210 s (T)'  ,'210 s (E)' , ...
           '270 s (T)'  ,'270 s (E)' , ...
           '330 s (T)'  ,'330 s (E)');

set(P,'Interpreter','Latex','Location','Northwest')


figure

plot(CmA(2,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(1,:))
hold on
plot(p120,  'Linewidth',    2,  'Color',    CT_M(2,:))
hold on
plot(CmA(4,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(3,:))
hold on
plot(p180,  'Linewidth',    2,  'Color',    CT_M(4,:))
hold on
plot(CmA(6,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(5,:))
hold on
plot(p240,  'Linewidth',    2,  'Color',    CT_M(6,:))
hold on
plot(CmA(8,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(7,:))
hold on
plot(p300,  'Linewidth',    2,  'Color',    CT_M(8,:))
hold on
plot(CmA(10,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(9,:))
hold on
plot(p360,  'Linewidth',    2,  'Color',    CT_M(10,:))

%corrcoef(p90(1:270),CmA(1,1:270))
%corrcoef(p120(1:300),CmA(2,1:300))
%corrcoef(p150(1:330),CmA(3,1:330))
%corrcoef(p180(1:360),CmA(4,1:360))
%corrcoef(p210(1:390),CmA(5,1:390))
%corrcoef(p240(1:420),CmA(6,1:420))
%corrcoef(p270(1:450),CmA(7,1:450))
%corrcoef(p300(1:480),CmA(8,1:480))
%corrcoef(p330(1:510),CmA(9,1:510))
%corrcoef(p360(1:540),CmA(10,1:540))

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

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Time [s]','Interpreter','Latex')

P = legend('120 s (T)'   ,'120 s (E)' , ...
           '180 s (T)'  ,'180 s (E)' , ...
           '240 s (T)'  ,'240 s (E)' , ...
           '300 s (T)'  ,'300 s (E)' , ...
           '360 s (T)'  ,'360 s (E)');

set(P,'Interpreter','Latex','Location','Northwest')

%Signal Amplitude ---------------------------------------------------------

figure

max_amp_exp = [max(p90)  max(p120) max(p150) max(p180) max(p210) ...
       max(p240) max(p270) max(p300) max(p330) max(p360)];

max_amp_theo = max(transpose(CmA));

plot(max_amp_exp,   'o',  'Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(max_amp_theo,'Linewidth',2)

%corrcoef(max_amp_exp,max_amp_theo)


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

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Pulse Width [s]','Interpreter','Latex')

xticks([1 2 3 4 5 6 7 8 9 10])
xticklabels({'90','120','150', '180','210','240','270','300','330','360'})
xlim([1 10])

figure

Background = [p90(270)  p120(300) p150(330) p180(360) p210(390) ...
              p240(420) p270(450) p300(480) p330(510) p360(540)];

CmA_Background = [CmA(1,270)  CmA(2,300) CmA(3,330) CmA(4,360) CmA(5,390) ...
              CmA(6,420) CmA(7,450) CmA(8,480) CmA(9,510) CmA(10,540)] ;         
          
plot(Background,   'o',  'Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(CmA_Background,'Linewidth',2)

%corrcoef(Background,CmA_Background)

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

ylabel('Leftover Signal [nA]','Interpreter','Latex')
xlabel('Pulse Width [s]','Interpreter','Latex')

xticks([1 2 3 4 5 6 7 8 9 10])
xticklabels({'90','120','150', '180','210','240','270','300','330','360'})
xlim([1 10])


% Signal Energy Analysis ==================================================

figure

Energy = [ sum(abs(transpose(p90(1:270)).^2))  ...
           sum(abs(transpose(p120(1:300)).^2)) ...
           sum(abs(transpose(p150(1:330)).^2)) ...
           sum(abs(transpose(p180(1:360)).^2)) ...
           sum(abs(transpose(p210(1:390)).^2)) ...
           sum(abs(transpose(p240(1:420)).^2)) ...
           sum(abs(transpose(p270(1:450)).^2)) ...
           sum(abs(transpose(p300(1:480)).^2)) ...
           sum(abs(transpose(p330(1:510)).^2)) ...
           sum(abs(transpose(p360(1:540)).^2)) ];

Energy_Theory = sum(abs(transpose(CmA).^2));

plot(Energy,        'o',  'Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(Energy_Theory,'Linewidth',2)

%corrcoef(Energy,Energy_Theory)

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

xticks([1 2 3 4 5 6 7 8 9 10])
xticklabels({'90','120','150', '180','210','240','270','300','330','360'})
xlim([1 10])

xlabel('Pulse Width [s]','Interpreter','Latex')
ylabel('Signal Energy [nW]','Interpreter','Latex')

K = legend('Experimental Results','Theoretical Model');

set(K,'Interpreter','Latex','Location','Northwest')

figure

SNR = 10*log10(Energy/sigma_a^2);

Smode = 10*log10(Energy_Theory/sigma_a^2);

plot(SNR,           'o',  'Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(Smode,'Linewidth',2)

%corrcoef(SNR,Smode)

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

xlabel('Pulse Width [s]','Interpreter','Latex')
ylabel('Signal-to-Noise ratio [dB]','Interpreter','Latex')

xticks([1 2 3 4 5 6 7 8 9 10])
xticklabels({'90','120','150', '180','210','240','270','300','330','360'})
xlim([1 10])

K = legend('Experimental Results','Theoretical Model');

set(K,'Interpreter','Latex','Location','Northwest')