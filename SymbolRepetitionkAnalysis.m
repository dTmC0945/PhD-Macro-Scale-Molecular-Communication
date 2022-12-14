% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- k- Repetititon by Daniel T. McGuiness --------------------------------
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

bitSequence                 = zeros(4, 22); bitRepetition   = zeros(4, 22);

CmNoisedk                   = cell(4,1);    kData            = cell(4,1);

% Model Parameters --------------------------------------------------------

diffusion       = 0.15;                 % diffusion coefficient (m2/s)
advectiveFlow   = 0.17;                 % advective flow (m/s)
distance        = 2.5;                  % transmission distance (m)
mass            = [0.6, 0.9, 1.2, 1.9]; % injected mass (ng)
symbolDuration  = 20;                   % symbol duration (s)

for i = 1:1:4
    bitSequence(i,:)    = mass(i).*[0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
    bitRepetition(i,:)  =       i.*[0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

end

mu_a    = 1.21*10^-3;                   % Mean
sigma_a = sqrt(0.0960*10^-6);           % Standard deviation

% Experimental Values -----------------------------------------------------

A = dataRead("kExperimentalData.xlsx"); % reading data using dataread

% Allocating the rows and removing rows and columns with text 

kExperimentalData = table2cell([A(2,15:end) ;  A(3,15:end) ;  A(4,15:end) ; ...
                                A(6,15:end) ;  A(7,15:end) ;  A(8,15:end) ; ...
                                A(10,15:end);  A(11,15:end);  A(12,15:end); ...
                                A(14,15:end);  A(15,15:end);  A(16,15:end)]);

% -------------------------------------------------------------------------

for kValue = 1:1:4      % Comparing theoretical and experimental values ---

    % Generating the theoretical transmission using transmission function

    CmNoisedk(kValue,:) = {addNoise(transmission(diffusion, advectiveFlow, ...
                                            distance, ....
                                            bitSequence(kValue,:), ...
                                            bitRepetition(kValue,:), ...
                                            symbolDuration), ...
                                            mu_a, sigma_a)};            
    
    % Getting the experimental data matching the kValue

    kData(kValue,:) = {mean(cell2mat(kExperimentalData(1 + 3*(kValue-1): ...
                                                       3 + 3*(kValue-1),:)))};

    % Comparing them by plotting

    kPlot(kData, CmNoisedk, ...
          kValue, ...
          diffusion, advectiveFlow, distance)

end

% converting cell Data to array (because the loop is keep giving me errors)

arraykExperimentalData = cell2mat(kExperimentalData);

for t = 1:1:11
    
    counter = 1;

    for i = 1:1:4

        % Generating the leftover values after flushing and its average

        leftoverTheoretical(i,t)            = CmNoisedk{i}(1,20*i + 40*i*(t-1));
        leftoverExperimentalAverage(i,t)    = kData{i}(1,20*i + 40*i*(t-1))*10^9;

        for j = 1:1:3
            
            leftoverExperimental(counter,t) = arraykExperimentalData(i*j,20*i + 40*i*(t-1))*10^9;
            
            counter = counter + 1;

        end
    end   
end

repetitionCount = 1:1:11;               % number of repetitions done

errBar  = [range([leftoverExperimental(1,:);  leftoverExperimental(2,:);  leftoverExperimental(3,:)]); ...
           range([leftoverExperimental(4,:);  leftoverExperimental(5,:);  leftoverExperimental(6,:)]); ...
           range([leftoverExperimental(7,:);  leftoverExperimental(8,:);  leftoverExperimental(9,:)]); ...
           range([leftoverExperimental(10,:); leftoverExperimental(11,:); leftoverExperimental(12,:)])]./2;

figure

for i = 1:1:4
    plot(leftoverTheoretical(i,:),'Linewidth',2 );                             hold on;
    errorbar(repetitionCount,leftoverExperimentalAverage(i,:),errBar(i),'o');  hold on;
end

% Nice plot code ---------------------------------------------------------
    
    set(gca, ...
      'TickLabelInterpreter'    , 'latex'       , ...
      'Fontsize'                , 16            , ...
      'TickDir'                 , 'in'          , ...
      'XGrid'                   , 'on'          , ...
      'YGrid'                   , 'on'          , ...
      'TickLength'              , [.02 .02]     , ...
      'GridLineStyle'           , '--'          , ...
      'YMinorTick'              , 'on'          , ...
      'Box'                     , 'on'          , ...
      'LineWidth'               , 1.5           );

 % ------------------------------------------------------------------------

xlabel('Repetition'             , 'Interpreter', 'Latex')
ylabel('Leftover Signal [nA]'   , 'Interpreter', 'Latex')

%title({'$D$ = 0.2 $\mathrm{cm^2/s}$, $u_x$ = 0.155 $\mathrm{cm/s}$, $x_d$ = 2.5 $\mathrm{cm}$'; '$ M_{k=1}$ = 0.6 $\mathrm{ng}$, $M_{k=2}$ = 0.9 $\mathrm{ng}$'; '$M_{k=3}$ = 1.2 $\mathrm{ng}$, $M_{k=4}$ = 1.9 $\mathrm{ng}$'},'Interpreter','Latex','Fontsize',16)

set(legend('$k$ = 1', '$k$ = 2', '$k$ = 3', '$k$ = 4'), ...
           'Interpreter', 'Latex', ...
           'Location'   , 'Northwest')

xlim([1 10])

% Correlation Values ------------------------------------------------------

for kValue = 1:1:4
    correlationValues(kValue) = findCorrelation(cell2mat(CmNoisedk(kValue,:)), ...
                                               (kData{kValue}(1:420*kValue)));
end

% ============================ FUNCTIONS ==================================

function kPlot(kData, CmNoisedk, kValue, diffusion, advectiveFlow, distance)

    tickValues = zeros(1, 7);

    figure
    plot(cell2mat(kData(kValue,:))*10^9  ,'Linewidth', 2);  hold on;
    plot(cell2mat(CmNoisedk(kValue,:))   ,'Linewidth', 2);
    
    % Nice plot code ------------------------------------------------------
    
    set(gca, ...
      'TickLabelInterpreter'    , 'latex'                       , ...
      'Fontsize'                , 16                            , ...
      'TickDir'                 , 'in'                          , ...
      'XGrid'                   , 'on'                          , ...
      'YGrid'                   , 'on'                          , ...
      'TickLength'              , [.02 .02]                     , ...
      'GridLineStyle'           , '--'                          , ...
      'YMinorTick'              , 'on'                          , ...
      'Box'                     , 'on'                          , ...
      'XTick'                   , 0:20*kValue:420* kValue       , ...
      'YTick'                   , 0:0.4:2                       , ...
      'LineWidth'               , 1.5                           );

    % ---------------------------------------------------------------------
    
    xlim([0 420*kValue]);   ylim([0 2]);
    
    for n = 1:1:7
        tickValues(n) = 60*n*kValue;
    end

    xticklabels({'0', '', '', tickValues(1), '', '', tickValues(2), ...
                      '', '', tickValues(3), '', '', tickValues(4), ...
                      '', '', tickValues(5), '', '', tickValues(6), ...
                      '', '', tickValues(7)})

    ylabel('Signal Current [nA]'    , 'Interpreter', 'Latex')
    xlabel('Transmission Time [s]'  , 'Interpreter', 'Latex')
    
   % title({"$D$ = " + diffusion + " $\mathrm{cm^2/s}$, $u_x$ =  $\mathrm{cm/s}$, $x_d$ = " distance + " $\mathrm{cm}$";" $M$ = 1.2 $\mathrm{ng}$, $\rho$ = 0.94"},'Interpreter','Latex','Fontsize',16)
    
    set(legend('Experimental Results','Theoretical Model'), ...
               'Interpreter','Latex', ...
               'Location','Northeast')
end

% ======================== END OF CODE ====================================