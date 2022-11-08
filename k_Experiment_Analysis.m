clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")
addpath("./Experimental Data")

opts = detectImportOptions('kExperimentalData.xlsx');
opts.DataRange = [1 Inf]; %rows;
cols = length(opts.VariableNames); %finding number of columns;
opts.SelectedVariableNames = [1, 1:cols];
A = readtable('kExperimentalData.xlsx',opts);

% Memory Allocation -------------------------------------------------------

leftoverTheoretical = zeros(4,11);
leftoverExperimentalAverage = zeros(4,11);
leftoverExperimental = zeros(12,11);

CmNoisedk = cell(4,1);
kData     = cell(4,1);

% Model Parameters --------------------------------------------------------

diffusion       = 0.15; 
advectiveFlow   = 0.17; 
distance        = 2.5;      % Parameters
mass            = [0.6, 0.9, 1.2, .19];
symbolDuration  = 20; 

% Noise Parameters --------------------------------------------------------

mu_a = 1.21*10^-3;              % Mean
sigma_a = sqrt(0.0960*10^-6);   % Standard deviation

% Experimental Values -----------------------------------------------------

kExperimentalData = table2cell([A(2,15:end) ;  A(3,15:end) ;  A(4,15:end) ; ...
                                 A(6,15:end) ;  A(7,15:end) ;  A(8,15:end) ; ...
                                 A(10,15:end);  A(11,15:end);  A(12,15:end); ...
                                 A(14,15:end);  A(15,15:end);  A(16,15:end)]);

% k = 1 Analysis ----------------------------------------------------------

bitSequence     = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]*mass(1);
bitRepetition   = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

Cm(1,:) = {transmission(diffusion, advectiveFlow, ...
                    distance, ....
                    bitSequence, bitRepetition, ...
                    symbolDuration)};

CmNoisedk(1,:) = {addNoise(cell2mat(Cm), mu_a, sigma_a)};            

%k = cell2mat(kExperimentalData(1,:)) + cell2mat(kExperimentalData(2,:));
p = mean(cell2mat(kExperimentalData(1:3,:)))
figure
plot(cell2mat(k(1,:)),'Linewidth',2)
hold on
plot(cell2mat(CmNoisedk(1,:)),'Linewidth',2)

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

CmNoisedk(2,:) = {addNoise(Cmk2, mu_a, sigma_a)};            

k{2} = sum([kExperimentalData{4:6}])/3;

figure
plot(cell2mat(k(2,:)),'Linewidth',2)
hold on
plot(cell2mat(CmNoisedk(2,:)),'Linewidth',2)

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
xticks(0:40:420*2)

xticklabels({'0','','','120','','','240','','','360','','','480','','','600','','','720','','','840'})

xlim([0 840])
title({'$D$ = 0.15 $\mathrm{cm^2/s}$, $u_x$ = 0.17 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 0.9 $\mathrm{ng}$, $\rho$ = 0.97'},'Interpreter','Latex','Fontsize',16)
ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
k2_legend = legend('Experimental Results','Theoretical Model');
set(k2_legend,'Interpreter','Latex','Location','Northeast')

% K = 3 Analysis ==========================================================

bitSequence = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]*1.2;
bitRepetition  = [0 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];

Cmk3 = transmission(diffusion, advectiveFlow, ...
                    distance, ....
                    bitSequence, bitRepetition, ...
                    symbolDuration);

CmNoisedk(3,:) = {addNoise(Cmk3, mu_a, sigma_a)};            

k{3} = sum([kExperimentalData{7:9}])/3;

figure
plot(cell2mat(k(3,:)),'Linewidth',2)
hold on
plot(cell2mat(CmNoisedk(3,:)),'Linewidth',2)

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
xticks(0:60:420*3)

xticklabels({'0','','','180','','','360','','','540','','','720','','','900','','','1080','','','1260'})

ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.15 $\mathrm{cm^2/s}$, $u_x$ = 0.17 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 1.2 $\mathrm{ng}$, $\rho$ = 0.94'},'Interpreter','Latex','Fontsize',16)

k3_legend = legend('Experimental Results','Theoretical Model');
set(k3_legend,'Interpreter','Latex','Location','Northeast')

% K = 4 Analysis ==========================================================

bitSequence     = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]*1.9;
bitRepetition   = [0 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4];

Cmk4 = transmission(diffusion, advectiveFlow, ...
                    distance, ....
                    bitSequence, bitRepetition, ...
                    symbolDuration);

CmNoisedk(4,:) = {addNoise(Cmk4, mu_a, sigma_a)};            

k{4} = sum([kExperimentalData{10:12}])/3;

figure
plot(cell2mat(k(4,:)),'Linewidth',2)
hold on
plot(cell2mat(CmNoisedk(4,:)),'Linewidth',2)

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
xticks(0:80:420*4)
xticklabels({'0','','','240','','','480','','','720','','','960','','','1200','','','1440','','','1680'})

xlim([0 420*4])

k4_legend = legend('Experimental Results','Theoretical Model');
set(k4_legend,'Interpreter','Latex','Location','Northeast')
ylabel('Signal Current [nA]','Interpreter','Latex')
xlabel('Transmission Time [s]','Interpreter','Latex')
title({'$D$ = 0.15 $\mathrm{cm^2/s}$, $u_x$ = 0.17 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$';' $M$ = 1.9 $\mathrm{ng}$, $\rho$ = 0.94'},'Interpreter','Latex','Fontsize',16)



for t = 1:1:11
    
    counter = 1;

    for i = 1:1:4

        leftoverTheoretical(i,t)            = CmNoisedk{i}(1,20*i + 40*i*(t-1));
        leftoverExperimentalAverage(i,t)    = k{i}(1,20*i + 40*i*(t-1));

        for j = 1:1:3
            
            leftoverExperimental(i*j,t)     = kExperimentalData{counter}(1,20*i + 40*i*(t-1));
            
            counter = counter + 1;
        end
    end

    
    leftoverExperimental(:,t) = [kExperimentalData(1, 20 + 40 *(t-1)); ...
                                 kExperimentalData(2, 20 + 40 *(t-1)); ...
                                 kExperimentalData(3, 20 + 40 *(t-1)); ...
                                 kExperimentalData(4, 40 + 80 *(t-1)); ...
                                 kExperimentalData(5, 40 + 80 *(t-1)); ...
                                 kExperimentalData(6, 40 + 80 *(t-1)); ...
                                 kExperimentalData(7, 60 + 120*(t-1)); ...
                                 kExperimentalData(8, 60 + 120*(t-1)); ...
                                 kExperimentalData(9, 60 + 120*(t-1)); ...
                                 kExperimentalData(10,80 + 160*(t-1)); ...
                                 kExperimentalData(11,80 + 160*(t-1)); ...
                                 kExperimentalData(12,80 + 160*(t-1))];
    

end



distance = 1:1:11;
x3 = [1 2 3 4 5 7 8 9 10 11];

errBar  = [range([leftoverExperimental(1,:);  leftoverExperimental(2,:);  leftoverExperimental(3,:)]); ...
           range([leftoverExperimental(4,:);  leftoverExperimental(5,:);  leftoverExperimental(6,:)]); ...
           range([leftoverExperimental(7,:);  leftoverExperimental(8,:);  leftoverExperimental(9,:)]); ...
           range([leftoverExperimental(10,:); leftoverExperimental(11,:); leftoverExperimental(12,:)])]./2;

%err_k_3(6) = [];
%Leftover_E_3(6) =[];

figure
for i = 1:1:4
plot(leftoverTheoretical(i,:),'Linewidth',2 )
hold on
end
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