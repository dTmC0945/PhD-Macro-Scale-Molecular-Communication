% QUANTITY ANALYSIS -------------------------------------------------------

clear all 
clc

A = xlsread('Quantity_Q.xlsx');         % Experimental Data

q1 = A(1,:)*10^9;                       % 1 ml/min
q2 = A(3,:)*10^9;                       % 2 ml/min
q3 = A(5,:)*10^9;                       % 3 ml/min
q4 = A(7,:)*10^9;                       % 4 ml/min
q5 = A(9,:)*10^9;                       % 5 ml/min
q6 = A(11,:)*10^9;                      % 6 ml/min
q7 = A(13,:)*10^9;                      % 7 ml/min
q8 = A(15,:)*10^9;                      % 8 ml/min
q9 = A(17,:)*10^9;                      % 9 ml/min
q10 = A(19,:)*10^9;                     % 10 ml/min

% Color definitions -------------------------------------------------------

CT = cbrewer('qual', 'Dark2', 10);

% -------------------------------------------------------------------------

plot(q1,	'Linewidth',    2,  'Color',    CT(1,:))
hold on
plot(q2,	'Linewidth',    2,  'Color',    CT(2,:))
hold on
plot(q3,	'Linewidth',    2,  'Color',    CT(3,:))
hold on
plot(q4,	'Linewidth',    2,  'Color',    CT(4,:))
hold on
plot(q5,	'Linewidth',    2,  'Color',    CT(5,:))
hold on
plot(q6,	'Linewidth',    2,  'Color',    CT(6,:))
hold on
plot(q7,	'Linewidth',    2,  'Color',    CT(7,:))
hold on
plot(q8,	'Linewidth',    2,  'Color',    CT(8,:))
hold on
plot(q9,	'Linewidth',    2,  'Color',    CT(9,:))
hold on
plot(q10,	'Linewidth',    2,  'Color',    CT(10,:))
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
  'XTick'       ,0:60:480   , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% ------------------------------------------------------------------------- 

xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

La = legend('1 ml/min', '2 ml/min', ...
            '3 ml/min', '4 ml/min', ...
            '5 ml/min', '6 ml/min', ...
            '7 ml/min', '8 ml/min', ...
            '9 ml/min', '10 ml/min');
       
set(La,'Interpreter','Latex','Location','Northwest')

xlim([0 480])


V = [ var(q1(90:390)) var(q2(90:390)) var(q3(90:390)) ...
      var(q4(90:390)) var(q5(90:390)) var(q6(90:390)) ...
      var(q7(90:390)) var(q8(90:390)) var(q9(90:390)) ...
      var(q10(90:390))];

figure % Signal variance plot

plot(V,'o',     'Linewidth'      ,1.2                       , ...
                'MarkerFaceColor',[0    0.4470    0.7410]   , ...
                'MarkerEdgeColor',[0    0.4470    0.7410])

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

% Signal variance plot fitting function

for i= 1:1:100
    
    Vr(i) = 0.001*(i/10)^2;
    
end

vr = 0.1:0.1:10;    % Signal variance plot function x-axis

hold on
plot(vr, Vr, 'Linewidth',2)

xlabel('Signal Flow [ml/min]','Interpreter','Latex')
ylabel('$\sigma_q^2$','Interpreter','Latex')

figure % Signal correlation plot

COR_12 = corrcoef(q1,q2);
COR_23 = corrcoef(q1,q3);
COR_34 = corrcoef(q1,q4);
COR_45 = corrcoef(q1,q5);
COR_56 = corrcoef(q1,q6);
COR_67 = corrcoef(q1,q7);
COR_78 = corrcoef(q1,q8);
COR_89 = corrcoef(q1,q9);
COR_910 = corrcoef(q1,q10);

COR = [ COR_12(1,2) COR_23(1,2) ...
        COR_34(1,2) COR_45(1,2) ...
        COR_56(1,2) COR_67(1,2) ...
        COR_78(1,2) COR_89(1,2) ...
        COR_910(1,2) ];

plot(COR,'o',   'Linewidth'         ,1.2                        , ...
                'MarkerFaceColor'   ,[0    0.4470    0.7410]    , ...
                'MarkerEdgeColor'   ,[0    0.4470    0.7410])

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

xlabel('Correlated Pairs ($q_1$, $q_{n+1}$)','Interpreter','Latex')
ylabel('$\rho (A, B)$','Interpreter','Latex')

figure

x = 2.5;
v = 0.02;
d = 0.15;
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
    
    Em = 0.12*i;
    
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

% Color definitions -------------------------------------------------------

CT_M = cbrewer('qual', 'Paired', 10);

% -------------------------------------------------------------------------

plot(CmA(2,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(1,:))
hold on
plot(q2,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(2,:))
hold on
plot(CmA(4,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(3,:))
hold on
plot(q4,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(4,:))
hold on
plot(CmA(6,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(5,:))
hold on
plot(q6,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(6,:))
hold on
plot(CmA(8,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(7,:))
hold on
plot(q8,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(8,:))
hold on
plot(CmA(10,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(9,:))
hold on
plot(q10,        'Linewidth' ,    2  , ...
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

P = legend('2 ml/min (T)'  ,'2 ml/min (E)' , ...
           '4 ml/min (T)'  ,'4 ml/min (E)' , ...
           '6 ml/min (T)'  ,'6 ml/min (E)' , ...
           '8 ml/min (T)'  ,'8 ml/min (E)' , ...
           '10 ml/min (T)' ,'10 ml/min (E)');

set(P,'Interpreter','Latex','Location','Northwest')

figure

plot(CmA(1,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(1,:))
hold on
plot(q1,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(2,:))
hold on
plot(CmA(3,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(3,:))
hold on
plot(q3,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(4,:))
hold on
plot(CmA(5,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(5,:))
hold on
plot(q5,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(6,:))
hold on
plot(CmA(7,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...
                'Color'     , CT_M(7,:))
hold on
plot(q7,        'Linewidth' ,    2  , ...
                'Color'     , CT_M(8,:))
hold on
plot(CmA(9,:),  'Linewidth' ,    2  , ...
                'LineStyle' ,   '--', ...   
                'Color'     , CT_M(9,:))
hold on
plot(q9,        'Linewidth' ,    2  , ...
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

%corrcoef(CmA(1,1:480),q1(1:480))
%corrcoef(CmA(2,1:480),q2(1:480))
%corrcoef(CmA(3,1:480),q3(1:480))
%corrcoef(CmA(4,1:480),q4(1:480))
%corrcoef(CmA(5,1:480),q5(1:480))
%corrcoef(CmA(6,1:480),q6(1:480))
%corrcoef(CmA(7,1:480),q7(1:480))
%corrcoef(CmA(8,1:480),q8(1:480))
%corrcoef(CmA(9,1:480),q9(1:480))
%corrcoef(CmA(10,1:480),q10(1:480))


xlabel('Time [s]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

P = legend('1 ml/min (T)'  ,'1 ml/min (E)' , ...
           '3 ml/min (T)'  ,'3 ml/min (E)' , ...
           '5 ml/min (T)'  ,'5 ml/min (E)' , ...
           '7 ml/min (T)'  ,'7 ml/min (E)' , ...
           '9 ml/min (T)'  ,'9 ml/min (E)');


set(P,'Interpreter','Latex','Location','Northwest')

figure

max_amp_exp = [max(q1) max(q2) ...
               max(q3) max(q4) ...
               max(q5) max(q6) ...
               max(q7) max(q8) ...
               max(q9) max(q10)];

i = 1;

for mass = 1:1:10
    
    Em = 0.12*mass;
    
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

i = i + 1;

end

x_axis = 1:1:10;

max_amp_theo = max(transpose(CmA));

plot(max_amp_exp,   'o',  'Linewidth'         ,    1.2                ,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(x_axis, max_amp_theo,'Linewidth',2)

corrcoef(max_amp_exp,max_amp_theo)


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

xlabel('Signal Flow [ml/min]','Interpreter','Latex')
ylabel('Signal Current [nA]','Interpreter','Latex')

K = legend('Experimental Results','Theoretical Model');

set(K,'Interpreter','Latex','Location','Northwest')

% Energy Analysis ---------------------------------------------------------

Energy_Theory = sum(abs(transpose(CmA).^2));

Energy = [sum(abs(q1.^2)) sum(abs(q2.^2)) ...
          sum(abs(q3.^2)) sum(abs(q4.^2)) ...
          sum(abs(q5.^2)) sum(abs(q6.^2)) ...
          sum(abs(q7.^2)) sum(abs(q8.^2)) ...
          sum(abs(q9.^2)) sum(abs(q10.^2))];


figure

plot(Energy,        'o',  'Linewidth'         ,    1.2,...
                          'MarkerFaceColor'   ,[0    0.4470    0.7410],...
                          'MarkerEdgeColor'   ,[0    0.4470    0.7410])
hold on
plot(x_axis, Energy_Theory,'Linewidth',2)

%corrcoef(Energy,Energy_Theory)


% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 15        , ...
  'TickDir'     , 'out'     , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'off'     , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Signal Flow [ml/min]','Interpreter','Latex')
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
plot(x_axis, Smode,'Linewidth',2)

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

xlabel('Signal Flow [ml/min]','Interpreter','Latex')
ylabel('Signal-to-Noise ratio [dB]','Interpreter','Latex')

K = legend('Experimental Results','Theoretical Model');

set(K,'Interpreter','Latex','Location','Northwest')
