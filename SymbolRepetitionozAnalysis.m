% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- oz- Repetititon by Daniel T. McGuiness -------------------------------
% -------------------------------------------------------------------------

clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")              % path where all used function are
addpath("./Experimental Data")      % path where experimental data are

% Memory Allocation -------------------------------------------------------

leftoverTheoretical         = zeros(4, 11); 
leftoverExperimentalAverage = zeros(4, 11);
leftoverExperimental        = zeros(12, 11);

bitSequence                 = zeros(4, 14); bitRepetition   = zeros(4, 14);

CmNoisedoz                   = cell(4,1);    ozData          = cell(4,1);

% Model Parameters --------------------------------------------------------

diffusion       = 0.20;                     % diffusion coefficient (m2/s)
advectiveFlow   = 0.155;                    % advective flow (m/s)
distance        = 2.5;                      % transmission distance (m)
mass            = [0.18, 0.32, 0.40, 0.45]; % injected mass (ng)
symbolDuration  = 20;                       % symbol duration (s)

% Noise Parameters ========================================================

mu_a = 1.21*10^-3;              % Mean
    
sigma_a = sqrt(0.0960*10^-6);   % Standard deviation

for i = 1:1:4
    bitSequence(i,:)    = mass(i).*[0 0 1 0 1 0 1 0 1 0 1 0 1 0];
    bitRepetition(i,:)  =       1.*[0 5 5 i 5 i 5 i 5 i 5 i 5 5];

end

mu_a    = 1.21*10^-3;                   % Mean
sigma_a = sqrt(0.0960*10^-6);           % Standard deviation

% Colour Data =============================================================

CT = cbrewer('qual', 'Paired', 10);

% Experimental Values -----------------------------------------------------

A = table2array(dataRead("ozExperimentalData.xlsx"))*10^9; % reading data using dataread

% Allocating the rows and removing rows and columns with text 

ozExperimentalData = [A(2,15:end);  A(3,15:end) ; ...
                      A(5,15:end);  A(6,15:end) ; ...
                      A(8,15:end);  A(9,15:end) ; ...
                      A(11,15:end); A(12,15:end); ...
                      A(14,15:end); A(15,15:end)];  

% -------------------------------------------------------------------------

for ozValue = 1:1:4      % Comparing theoretical and experimental values ---

    % Generating the theoretical transmission using transmission function

    CmNoisedoz(ozValue,:) = {addNoise(transmission(diffusion, advectiveFlow, ...
                                                   distance, ....
                                                   bitSequence(ozValue,:), ...
                                                   bitRepetition(ozValue,:), ...
                                                   symbolDuration), ...
                                                   mu_a, sigma_a)};            
    
    % Getting the experimental data matching the kValue

    kData(ozValue,:) = {mean((ozExperimentalData(1 + 3*(ozValue-1): ...
                                                         3 + 3*(ozValue-1),:)))};

    % Comparing them by plotting

    kPlot(kData, CmNoisedoz, ...
          ozValue, ...
          diffusion, advectiveFlow, distance)

end

% oz = 1 Analysis =========================================================

message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0]*0.18;
k  =      [0 5 5 1 5 1 5 1 5 1 5 1 5 5];

Tt      = 20; Ts      = 0; Cm      = 0;

for j = 2:1:length(message)
        
        if message(j) > message(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + (message(j) - message(j-1)) - ((message(j) - message(j-1))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t)))))/2 ;
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (abs(M_r - message(j)))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))))/2 + message(j);
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
Cm = Cm + WGN;     

L_T_1 = Cm;

OZ_1 = (oz_1_1 + oz_1_2)/2;

figure

plot(OZ_1,'Linewidth',2)
hold on
plot(Cm,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'     , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'TickLength'  , [.02 .02] , ...
  'GridLineStyle','--'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'     , ...
  'XTick'       ,80:80:820 ,...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.2 $\mathrm{cm^2/s}$, $u_x$ = 0.155 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 0.18 $\mathrm{ng}$, $\rho$ = 0.95'},'Interpreter','Latex','Fontsize',16)

oz1_legend = legend('Experimental Results','Theoretical Model');
set(oz1_legend,'Interpreter','Latex','Location','Southeast')

xlim([80 820])

xticks([80:37:820])

xticklabels({'0','','','','148','','','','296',...
                 '','','','444','','','','592','','','','740'})

% oz = 2 Analysis =========================================================

message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0]*0.32;
k  =      [0 5 5 2 5 2 5 2 5 2 5 2 5 5];

Tt      = 20; Ts      = 0; Cm      = 0;

for j = 2:1:length(message)
        
        if message(j) > message(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + (message(j) - message(j-1)) - ((message(j) - message(j-1))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t)))))/2 ;
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (abs(M_r - message(j)))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))))/2 + message(j);
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
Cm = Cm + WGN;     

L_T_2 = Cm;

OZ_2 = (oz_2_1 + oz_2_2)/2;


figure

plot(OZ_2,'Linewidth',2)
hold on
plot(Cm,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'     , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'TickLength'  , [.02 .02] , ...
  'GridLineStyle','--'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'     , ...
  'XTick'       ,80:43:940 ,...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------

xlim([80 940])

xticks([80:43:940])

xticklabels({'0','','','','172','','','','344',...
                 '','','','516','','','','688','','','','860'})

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.2 $\mathrm{cm^2/s}$, $u_x$ = 0.155 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 0.32 $\mathrm{ng}$, $\rho$ = 0.97'},'Interpreter','Latex','Fontsize',16)

oz2_legend = legend('Experimental Results','Theoretical Model');
set(oz2_legend,'Interpreter','Latex','Location','Southeast')

% oz = 3 Analysis =========================================================

message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0]*0.4;
k  =      [0 5 5 3 5 3 5 3 5 3 5 3 5 5];

Tt      = 20; M_old   = 0; Ts      = 0; Q_neg   = 0; Cm      = 0;

for j = 2:1:length(message)
        
        if message(j) > message(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + (message(j) - message(j-1)) - ((message(j) - message(j-1))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t)))))/2 ;
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (abs(M_r - message(j)))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))))/2 + message(j);
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
Cm = Cm + WGN;     

L_T_3 = Cm;

OZ_3 = (oz_3_1 + oz_3_2)/2;


figure


plot(OZ_3,'Linewidth',2)
hold on
plot(Cm,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'     , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'TickLength'  , [.02 .02] , ...
  'GridLineStyle','--'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'     , ...
  'XTick'       ,80:49:1060 ,...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------

xlim([80 1060])

xticks([80:49:1060])

xticklabels({'0','','','','196','','','','392',...
                 '','','','588','','','','784','','','','980'})
             
             
ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.2 $\mathrm{cm^2/s}$, $u_x$ = 0.155 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 0.4 $\mathrm{ng}$, $\rho$ = 0.97'},'Interpreter','Latex','Fontsize',16)

oz3_legend = legend('Experimental Results','Theoretical Model');
set(oz3_legend,'Interpreter','Latex','Location','Southeast')
% oz = 4 Analysis =========================================================

message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0]*0.45;
k  =      [0 5 5 4 5 4 5 4 5 4 5 4 5 5];

Tt      = 20; M_old   = 0; Ts      = 0; Q_neg   = 0; Cm      = 0;

for j = 2:1:length(message)
        
        if message(j) > message(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + (message(j) - message(j-1)) - ((message(j) - message(j-1))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t)))))/2 ;
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (abs(M_r - message(j)))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))))/2 + message(j);
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
Cm = Cm + WGN;     

L_T_4 = Cm;

OZ_4 = (oz_4_1 + oz_4_2)/2;

figure


plot(OZ_4,'Linewidth',2)
hold on
plot(Cm,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'     , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'TickLength'  , [.02 .02] , ...
  'GridLineStyle','--'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'     , ...
  'XTick'       ,80:55:1180 ,...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------

xlim([80 1180])

xticks([80:55:1180])

xticklabels({'0','','','','220','','','','440',...
                 '','','','660','','','','880','','','','1100'})


oz1_legend = legend('Experimental Results','Theoretical Model');
set(oz1_legend,'Interpreter','Latex','Location','Southeast')

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.2 $\mathrm{cm^2/s}$, $u_x$ = 0.155 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 0.45 $\mathrm{ng}$, $\rho$ = 0.96'},'Interpreter','Latex','Fontsize',16)

% oz = 5 Analysis =========================================================

message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0]*0.32;
k  =      [0 5 5 5 5 5 5 5 5 5 5 5 5 5];

Tt      = 20; M_old   = 0; Ts      = 0; Q_neg   = 0; Cm      = 0;

for j = 2:1:length(message)
        
        if message(j) > message(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + (message(j) - message(j-1)) - ((message(j) - message(j-1))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t)))))/2 ;
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (abs(M_r - message(j)))*(erf((x-(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))) + ...
                                                                 erf((x+(t)*v)/(2*d*sqrt(1/(d*(t)))*(t))))/2 + message(j);
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
Cm = Cm + WGN;     

L_T_5 = Cm;

OZ_5 = (oz_5_1 + oz_5_2)/2;

figure

plot(OZ_5,'Linewidth',2)
hold on
plot(Cm,'Linewidth',2)
hold on
plot(WGN,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'     , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'TickLength'  , [.02 .02] , ...
  'GridLineStyle','--'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'     , ...
  'XTick'       ,80:61:1300 ,...
  'YTick'       ,0:0.05:035,...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------
xlim([80 1300])

xticks([80:61:1300])

xticklabels({'0','','','','244','','','','488',...
                 '','','','732','','','','976','','','','1220'})


ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.2 $\mathrm{cm^2/s}$, $u_x$ = 0.155 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 0.32 $\mathrm{ng}$, $\rho$ = 0.96'},'Interpreter','Latex','Fontsize',16)

oz5_legend = legend('Experimental Results','Theoretical Model','Environment Noise');
set(oz5_legend,'Interpreter','Latex','Location','Southeast')

for t = 1:1:7
    
    Leftover_E_1(t) = OZ_1(100 + 120*(t-1));
    Leftover_E_1_1(t) = oz_1_1(100 + 120*(t-1));
    Leftover_E_1_2(t) = oz_1_2(100 + 120*(t-1));
    
    Leftover_E_2(t) = OZ_2(100 + 140*(t-1));
    Leftover_E_2_1(t) = oz_2_1(100 + 140*(t-1));
    Leftover_E_2_2(t) = oz_2_2(100 + 140*(t-1));
    
    Leftover_E_3(t) = OZ_3(100 + 160*(t-1));
    Leftover_E_3_1(t) = oz_3_1(100 + 160*(t-1));
    Leftover_E_3_2(t) = oz_3_2(100 + 160*(t-1));
    
    Leftover_E_4(t) = OZ_4(100 + 180*(t-1));
    Leftover_E_4_1(t) = oz_4_1(100 + 180*(t-1));
    Leftover_E_4_2(t) = oz_4_2(100 + 180*(t-1));
    
    Leftover_E_5(t) = OZ_5(100 + 200*(t-1));
    Leftover_E_5_1(t) = oz_5_1(100 + 200*(t-1));
    Leftover_E_5_2(t) = oz_5_2(100 + 200*(t-1));
    
    Leftover_T_1(t) = L_T_1(100 + 120*(t-1));
    Leftover_T_2(t) = L_T_2(100 + 140*(t-1));
    Leftover_T_3(t) = L_T_3(100 + 160*(t-1));
    Leftover_T_4(t) = L_T_4(100 + 180*(t-1));
    Leftover_T_5(t) = L_T_5(100 + 200*(t-1));
end
x = 1:1:7;
err_oz_1 = range([Leftover_E_1_1 ; Leftover_E_1_2])./2;
err_oz_2 = range([Leftover_E_2_1 ; Leftover_E_2_2])./2;
err_oz_3 = range([Leftover_E_3_1 ; Leftover_E_3_2])./2;
err_oz_4 = range([Leftover_E_4_1 ; Leftover_E_4_2])./2;
err_oz_5 = range([Leftover_E_5_1 ; Leftover_E_5_2])./2;

figure
plot(Leftover_T_1,'Linewidth',2,'Color',CT(2,:))
hold on
plot(Leftover_T_2,'Linewidth',2,'Color',CT(4,:))
hold on
plot(Leftover_T_3,'Linewidth',2,'Color',CT(6,:))
hold on
plot(Leftover_T_4,'Linewidth',2,'Color',CT(8,:))
hold on
plot(Leftover_T_5,'Linewidth',2,'Color',CT(10,:))
hold on

he1 = errorbar(x,Leftover_E_1,err_oz_1,'o');
hold on
he2 = errorbar(x,Leftover_E_2,err_oz_2,'o');
hold on
he3 = errorbar(x,Leftover_E_3,err_oz_3,'o');
hold on
he4 = errorbar(x,Leftover_E_4,err_oz_4,'o');
hold on
he5 = errorbar(x,Leftover_E_5,err_oz_5,'o');
hold on


set(he1                        , ...
  'LineWidth'       , 1.5      , ...
  'Marker'          , 'o'      , ...
  'Color'           , CT(1,:)  , ...
  'MarkerSize'      , 8        , ...
  'MarkerEdgeColor' , CT(1,:)  , ...
  'MarkerFaceColor' , CT(1,:)   );
set(he2                        , ...
  'LineWidth'       , 1.5      , ...
  'Marker'          , 'o'      , ...
  'Color'           , CT(3,:)  , ...
  'MarkerSize'      , 8        , ...
  'MarkerEdgeColor' , CT(3,:)  , ...
  'MarkerFaceColor' , CT(3,:)  );
set(he3                        , ...
  'LineWidth'       , 1.5      , ...
  'Marker'          , 'o'      , ...
  'Color'           , CT(5,:)  , ...
  'MarkerSize'      , 8        , ...
  'MarkerEdgeColor' , CT(5,:)  , ...
  'MarkerFaceColor' , CT(5,:)  );
set(he4                        , ...
  'LineWidth'       , 1.5      , ...
  'Marker'          , 'o'      , ...
  'Color'           , CT(7,:)  , ...
  'MarkerSize'      , 8        , ...
  'MarkerEdgeColor' , CT(7,:)  , ...
  'MarkerFaceColor' , CT(7,:)  );
set(he5                        , ...
  'LineWidth'       , 1.5      , ...
  'Marker'          , 'o'      , ...
  'Color'           , CT(9,:)  , ...
  'MarkerSize'      , 8        , ...
  'MarkerEdgeColor' , CT(9,:)  , ...
  'MarkerFaceColor' , CT(9,:)  );

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
xlabel('Repetition','Interpreter','Latex')
xlim([1 7])
ylim([0 0.15])

Background_legend = legend('$o/z$ = 5/1','$o/z$ = 5/2','$o/z$ = 5/3','$o/z$ = 5/4','$o/z$ = 5/5');
set(Background_legend,'Interpreter','Latex','Location','Northwest')
%title({'$D$ = 0.2 $\mathrm{cm^2/s}$, $u_x$ = 0.155 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$'; '$M_{5/1}$ = 0.18 $\mathrm{ng}$, $M_{5/2}$ = 0.32 $\mathrm{ng}$';'$M_{5/3}$ = 0.4 $\mathrm{ng}$, $M_{5/4}$ = 0.45 $\mathrm{ng}$, $M_{5/5}$ = 0.32 $\mathrm{ng}$'},'Interpreter','Latex','Fontsize',16)


% Correlation Values ======================================================

rho_oz_1 = corrcoef(L_T_1,OZ_1(1:length(L_T_1)));
rho_oz_2 = corrcoef(L_T_2,OZ_2(1:length(L_T_2)));
rho_oz_3 = corrcoef(L_T_3,OZ_3(1:length(L_T_3)));
rho_oz_4 = corrcoef(L_T_4,OZ_4(1:length(L_T_4)));
rho_oz_5 = corrcoef(L_T_5,OZ_5(1:length(L_T_5)));