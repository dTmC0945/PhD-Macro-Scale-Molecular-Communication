clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")
addpath("./Experimental Data")

A = xlsread('kExperimentalData.xlsx');

% Model Parameters ========================================================

diffusion = 0.15; 
advectiveFlow = 0.17; 
distance = 2.5;      % Parameters
mass = 0.6;
symbolDuration      = 20; 

% Noise Parameters ========================================================

mu_a = 1.21*10^-3;              % Mean
sigma_a = sqrt(0.0960*10^-6);   % Standard deviation

% Experimental Values -----------------------------------------------------

kExperimentalData = [A(1,15:end) ;  A(2,15:end) ;  A(3,15:end) ; ...
                     A(5,15:end) ;  A(6,15:end) ;  A(7,15:end) ; ...
                     A(9,15:end) ;  A(10,15:end);  A(11,15:end); ...
                     A(13,15:end);  A(14,15:end);  A(15,15:end)];

% k = 1 Analysis ----------------------------------------------------------

bitSequence     = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]*mass;
bitRepetition   = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

Cmk1 = transmission(diffusion, advectiveFlow, ...
                    distance, ....
                    bitSequence, bitRepetition, ...
                    symbolDuration);

CmNoisedk1 = addNoise(Cmk1, mu_a, sigma_a);            

k1 = (kExperimentalData(1,:) ...
    + kExperimentalData(2,:) ...
    + kExperimentalData(3,:))/3;

plot(k1,'Linewidth',2)
hold on
plot(CmNoisedk1,'Linewidth',2)

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
  'XTick'       ,0:20:420 ,...
  'YTick'       ,0:0.1:0.7, ...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------

xlim([0 420])
ylim([0 0.7])


xticklabels({'0','','','60','','','120','','','180','','','240','','','300','','','360','','','420'})

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.15 $\mathrm{cm^2/s}$, $u_x$ = 0.17 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 0.6 $\mathrm{ng}$, $\rho$ = 0.96'},'Interpreter','Latex','Fontsize',16)

set(legend('Experimental Results','Theoretical Model'),'Interpreter','Latex','Location','Northeast')

% K = 2 Analysis ==========================================================

bitSequence     = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]*0.9;
bitRepetition   = [0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

Cmk2 = transmission(diffusion, advectiveFlow, ...
                    distance, ....
                    bitSequence, bitRepetition, ...
                    symbolDuration);

CmNoisedk2 = addNoise(Cmk2, mu_a, sigma_a);            

k2 = (kExperimentalData(4,:) ...
    + kExperimentalData(5,:) ...
    + kExperimentalData(6,:))/3;

figure
plot(k2,'Linewidth',2)
hold on
plot(CmNoisedk2,'Linewidth',2)

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
  'XTick'       ,40:40:840 ,...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------
xticks([0:40:420*2])

xticklabels({'0','','','120','','','240','','','360','','','480','','','600','','','720','','','840'})

xlim([0 840])
title({'$D$ = 0.15 $\mathrm{cm^2/s}$, $u_x$ = 0.17 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 0.9 $\mathrm{ng}$, $\rho$ = 0.97'},'Interpreter','Latex','Fontsize',16)
ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
k2_legend = legend('Experimental Results','Theoretical Model');
set(k2_legend,'Interpreter','Latex','Location','Northeast')

% K = 3 Analysis ==========================================================

bitSequence = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]*1.2;
k  =      [0 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];

Tt      = 20; M_old   = 0; Ts      = 0; Q_neg   = 0; Cm      = 0;

for j = 2:1:length(bitSequence)
        
        if bitSequence(j) > bitSequence(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + (bitSequence(j) - bitSequence(j-1)) - ((bitSequence(j) - bitSequence(j-1))*(erf((distance-(t)*advectiveFlow)/(2*diffusion*sqrt(1/(diffusion*(t)))*(t))) + ...
                                                                 erf((distance+(t)*advectiveFlow)/(2*diffusion*sqrt(1/(diffusion*(t)))*(t)))))/2 ;
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (abs(M_r - bitSequence(j)))*(erf((distance-(t)*advectiveFlow)/(2*diffusion*sqrt(1/(diffusion*(t)))*(t))) + ...
                                                                 erf((distance+(t)*advectiveFlow)/(2*diffusion*sqrt(1/(diffusion*(t)))*(t))))/2 + bitSequence(j);
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
Cm = Cm + WGN;     

L_T_3 = Cm;


K_3 = (k_3_1 + k_3_2 + k_3_3)/3;

figure

plot(K_3,'Linewidth',2)
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
  'XTick'       ,60:60:420*3 ,...
  'YTick'       ,0:0.2:1.4,...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------

xlim([0 420*3])
xticks([0:60:420*3])

xticklabels({'0','','','180','','','360','','','540','','','720','','','900','','','1080','','','1260'})

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.15 $\mathrm{cm^2/s}$, $u_x$ = 0.17 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 1.2 $\mathrm{ng}$, $\rho$ = 0.94'},'Interpreter','Latex','Fontsize',16)

k3_legend = legend('Experimental Results','Theoretical Model');
set(k3_legend,'Interpreter','Latex','Location','Northeast')

% K = 4 Analysis ==========================================================

bitSequence = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]*1.9;
k  =      [0 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4];

Tt      = 20; M_old   = 0; Ts      = 0; Q_neg   = 0; Cm      = 0;

for j = 2:1:length(bitSequence)
        
        if bitSequence(j) > bitSequence(j-1)
            
            M_old = Cm(end);
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = M_old + (bitSequence(j) - bitSequence(j-1)) - ((bitSequence(j) - bitSequence(j-1))*(erf((distance-(t)*advectiveFlow)/(2*diffusion*sqrt(1/(diffusion*(t)))*(t))) + ...
                                                                 erf((distance+(t)*advectiveFlow)/(2*diffusion*sqrt(1/(diffusion*(t)))*(t)))))/2 ;
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        else
                M_r = Cm(end); 
            
            for t = 1:1:k(j)*Tt
                
                Q(t) = (abs(M_r - bitSequence(j)))*(erf((distance-(t)*advectiveFlow)/(2*diffusion*sqrt(1/(diffusion*(t)))*(t))) + ...
                                                                 erf((distance+(t)*advectiveFlow)/(2*diffusion*sqrt(1/(diffusion*(t)))*(t))))/2 + bitSequence(j);
                    
                Cm(t + Ts) = Q(t); 
                
            end
            
        end
        
        Ts = k(j)*Tt + Ts;
       
end
        
WGN = transpose(sigma_a.*randn(length(Cm),1) + mu_a);    % Additive
    
Cm = Cm + WGN;     

L_T_4 = Cm;


K_4 = (k_4_1 + k_4_2 + k_4_3)/3;

figure

plot(K_4,'Linewidth',2)
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
  'XTick'       ,80:80:420*4 ,...
  'LineWidth'   , 1.5         );
% -------------------------------------------------------------------------
xticks([0:80:420*4])
xticklabels({'0','','','240','','','480','','','720','','','960','','','1200','','','1440','','','1680'})

xlim([0 420*4])



k4_legend = legend('Experimental Results','Theoretical Model');
set(k4_legend,'Interpreter','Latex','Location','Northeast')
ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.15 $\mathrm{cm^2/s}$, $u_x$ = 0.17 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 1.9 $\mathrm{ng}$, $\rho$ = 0.94'},'Interpreter','Latex','Fontsize',16)


for t = 1:1:11
    
    Leftover_E_1(t) = k1(20 + 40*(t-1));
    Leftover_E_1_1(t) = k_1_1(20 + 40*(t-1));
    Leftover_E_1_2(t) = k_1_2(20 + 40*(t-1));
    Leftover_E_1_3(t) = k_1_3(20 + 40*(t-1));
    
    Leftover_E_2(t) = K_2(40 + 80*(t-1));
    Leftover_E_2_1(t) = k_2_1(40 + 80*(t-1));
    Leftover_E_2_2(t) = k_2_2(40 + 80*(t-1));
    Leftover_E_2_3(t) = k_2_3(40 + 80*(t-1));
    
    Leftover_E_3(t) = K_3(60 + 120*(t-1));
    Leftover_E_3_1(t) = k_3_1(60 + 120*(t-1));
    Leftover_E_3_2(t) = k_3_2(60 + 120*(t-1));
    Leftover_E_3_3(t) = k_3_3(60 + 120*(t-1));    
    
    Leftover_E_4(t) = K_4(80 + 160*(t-1));
    Leftover_E_4_1(t) = k_4_1(80 + 160*(t-1));
    Leftover_E_4_2(t) = k_4_2(80 + 160*(t-1));
    Leftover_E_4_3(t) = k_4_3(80 + 160*(t-1)); 
    
    Leftover_T_1(t) = L_T_1(20 + 40*(t-1));
    Leftover_T_2(t) = L_T_2(40 + 80*(t-1));
    Leftover_T_3(t) = L_T_3(60 + 120*(t-1));
    Leftover_T_4(t) = L_T_4(80 + 160*(t-1));
end
distance = 1:1:11;
x3 = [1 2 3 4 5 7 8 9 10 11];

err_k_1 = range([Leftover_E_1_1 ; Leftover_E_1_2; Leftover_E_1_3])./2;
err_k_2 = range([Leftover_E_2_1 ; Leftover_E_2_2; Leftover_E_2_3])./2;
err_k_3 = range([Leftover_E_3_1 ; Leftover_E_3_2; Leftover_E_3_3])./2;
err_k_4 = range([Leftover_E_4_1 ; Leftover_E_4_2; Leftover_E_4_3])./2;

err_k_3(6) = [];
Leftover_E_3(6) =[];

CT_M = cbrewer('qual', 'Paired', 10);


figure

plot(Leftover_T_1,'Linewidth',2 )
hold on
plot(Leftover_T_2,'Linewidth',2)
hold on
plot(Leftover_T_3,'Linewidth',2 )
hold on
plot(Leftover_T_4,'Linewidth',2  )
hold on
he1 = errorbar(distance,Leftover_E_1,err_k_1,'o');
hold on
he2 = errorbar(distance,Leftover_E_2,err_k_2,'o');
hold on
he3 = errorbar(x3,Leftover_E_3,err_k_3,'o', 'Color'     , CT_M(1,:));
hold on
he4 = errorbar(distance,Leftover_E_4,err_k_4,'o', 'Color'     , CT_M(1,:));
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
%title({'$D$ = 0.2 $\mathrm{cm^2/s}$, $u_x$ = 0.155 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$'; '$ M_{k=1}$ = 0.6 $\mathrm{ng}$, $M_{k=2}$ = 0.9 $\mathrm{ng}$'; '$M_{k=3}$ = 1.2 $\mathrm{ng}$, $M_{k=4}$ = 1.9 $\mathrm{ng}$'},'Interpreter','Latex','Fontsize',16)

Background_legend = legend('$k$ = 1','$k$ = 2','$k$ = 3','$k$ = 4');
set(Background_legend,'Interpreter','Latex','Location','Northwest')
xlim([1 10])


% Correlation Values ======================================================

rho_k_1 = corrcoef(L_T_1,k1(1:length(L_T_1)));
rho_k_2 = corrcoef(L_T_2,K_2(1:length(L_T_2)));
rho_k_3 = corrcoef(L_T_3,K_3(1:length(L_T_3)));
rho_k_4 = corrcoef(L_T_4,K_4(1:length(L_T_4)));