% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- M-ary Analysis by Daniel T. McGuiness --------------------------------
% -------------------------------------------------------------------------

clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")
addpath("./Experimental Data")

% Anaysing experimental data (%) ==========================================

A = maryDataCruncher('MaryExperimentalData.xlsx');

M230 = A(1,:); M430 = A(2,:); M830 = A(3,:);
M260 = A(4,:); M460 = A(5,:); M860 = A(6,:);
M290 = A(7,:); M490 = A(8,:); M890 = A(9,:);

% Model Parameters ========================================================

D   = 0.2; v = 0.2; x = 2.5;            % Parameters

m   = [0.0785, 0.0270, 0.0080; ...      % the mass matrix for the exp.
       0.2800, 0.0800, 0.0160; ...
       0.1300, 0.0390, 0.0120];

Tt  = [30 60 90];                       % time matrix for symbol length
    
% Noise Parameters ========================================================

mu_a = 1.21*10^-3;                      % Mean
    
sigma_a = sqrt(0.0960*10^-6);           % Standard deviation

% M2 - 30 =================================================================

message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0] * m(1,1);   % symbol types
k       = [0 1 1 2 2 1 1 1 1 4 2 1 1 1 2 3];            % repeating symbols

% physical transmission which calls the defined function "transmission".

Cm      = transmission(D, v, x, message, k, Tt(1));   
      
% Adds white gaussian noise to the defined transmission.

Cm230   = addNoise(Cm, mu_a, sigma_a);     

% To find correlation the following self defined function is used.
   
Corr230 = findCorrelation(M230(45:764), Cm230);      
                                            
% Plots the data and compares theoretical model to experimental data.                                    

maryPlotData(2, M230(45:765), Cm230 , ...
         6, Tt(1) , 0.09, ...
         0 ,30 ,800, ...
         length(Cm), 0.12)

xticklabels({'0','','','90' ,'','','180','','','270' ...
                ,'','','360','','','450','','','540' ...
                ,'','','630','','','720'})

% M4 - 30 =================================================================


message = [0 1 0 3 1 0 3 1 2 0] * m(1,2);
k       = [0 1 1 1 2 2 1 2 1 1];

Cm      = transmission(D, v, x, message, k, Tt(1));
        
Cm430   = addNoise(Cm, mu_a, sigma_a);

maryPlotData(4, M430(45:405), Cm430 , ...
         8, Tt(1) , 0.09, ...
         0 ,30 ,800, ...
         length(Cm), 0.12)

Corr430 = findCorrelation(M430(45:404), Cm430);

xticklabels({'0','','60' ,'','120','','180' ...
                ,'','240','','300','','360'})

% M8 - 30 =================================================================

message = [0 2 3 2 4 1 5 3 0] * m(1,3);
k       = [0 1 1 1 1 1 1 1 1];

Cm      = transmission(D, v, x, message, k, Tt(1));

Cm830   = addNoise(Cm, mu_a, sigma_a);
       
maryPlotData(8, M830(45:285), Cm830 , ...
         8, Tt(1) , 0.042, ...
         0 ,30 ,800, ...
         length(Cm), 0.06)

Corr830 = findCorrelation(M830(45:284), Cm830);

% M2 - 60s ================================================================

message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0] * m(2,1);
k       = [0 1 1 2 2 1 1 1 1 4 2 1 1 1 2 3];

Cm      = transmission(D, v, x, message, k, Tt(2));
       
Cm260   = addNoise(Cm, mu_a, sigma_a);
   
Corr260 = findCorrelation(M260(75:1514), Cm260);

maryPlotData(2, M260(75:1515), Cm260 , ...
         15, Tt(2) , 0.32, ...
         0 ,60 ,1440, ...
         length(Cm), 0.45)

xticklabels({'0','','','180' ,'','','360','','','540' , ...
                 '','','720' ,'','','900','','','1080', ...
                 '','','1260','','','1440'})

% M4 - 60s ================================================================

message = [0 1 0 3 1 0 3 1 2 0] * m(2,2);
k       = [0 1 1 1 2 2 1 2 1 1];

Cm      = transmission(D, v, x, message, k, Tt(2));
        
Cm460   = addNoise(Cm, mu_a, sigma_a);
   
Corr460 = findCorrelation(M460(75:794), Cm460);

maryPlotData(4, M460(75:795), Cm460 , ...
         15, Tt(2) , 0.28, ...
         0 ,60 ,720, ...
         length(Cm), 0.4)

xticklabels({'0','','120','','240','','360', ...
                 '','480','','600','','720'})

% M8 - 60s ================================================================

message = [0 2 3 2 4 1 5 3 0] * m(2,3);
k       = [0 1 1 1 1 1 1 1 1];

Cm      = transmission(D, v, x, message, k, Tt(2));

Cm860   = addNoise(Cm, mu_a, sigma_a);
   
Corr860 = findCorrelation(M860(75:554), Cm860);

maryPlotData(8, M860(75:555), Cm860 , ...
         15, Tt(2) , 0.09, ...
         0 ,60 ,480, ...
         length(Cm), 0.12)

% M2 - 90s ================================================================

message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0] * m(3,1);
k       = [0 1 1 2 2 1 1 1 1 4 2 1 1 1 2 3];

Cm      = transmission(D, v, x, message, k, Tt(3));

Cm290   = addNoise(Cm, mu_a, sigma_a);
   
Corr290 = findCorrelation(M290(105:2264), Cm290);

maryPlotData(2, M290(105:2265), Cm290 , ...
         22.5, Tt(3), 0.14, ...
         0 ,90 ,2160, ...
         length(Cm), 0.2)

xticklabels({'0','','','270' ,'','','540' ,'','','810' , ...
                 '','','1080','','','1350','','','1620', ...
                 '','','1890','','','2160'})

% M4 - 90s ================================================================

message = [0 1 0 3 1 0 3 1 2 0] * m(3,2);
k       = [0 1 1 1 2 2 1 2 1 1];

Cm      = transmission(D, v, x, message, k, Tt(3));
        
Cm490   = addNoise(Cm, mu_a, sigma_a);
   
Corr490 = findCorrelation(M490(105:1184), Cm490);

maryPlotData(4, M490(105:1185), Cm490 , ...
         22.5, Tt(3) , 0.13, ...
         0 ,90 ,1080, ...
         length(Cm), 0.2)

xticklabels({'0','','180','','360','','540' , ...
                 '','720','','900','','1080'})

% M8 - 90s ================================================================

message = [0 2 3 2 4 1 5 3 0] * m(3,3);
k       = [0 1 1 1 1 1 1 1 1];

Cm      = transmission(D, v, x, message, k, Tt(3));
        
Cm890   = addNoise(Cm, mu_a, sigma_a);
   
Corr890 = findCorrelation(M890(105:824), Cm890);

maryPlotData(8, M890(105:825), Cm890 , ...
         22.5, Tt(3) , 0.07, ...
         0 ,90 ,1080, ...
         length(Cm), 0.1)

function maryPlotData(mod, expData, Cm, iT, Tt, offset, min, mid, max, xmax, ymax)
    
    % -------------------------------------------------------------------------
    % Function which prints the experimental data and the theoretical data in
    % comparison
    % Parameters are
    % mod: the modulation level of the transmission (2, 4, 8)
    % expData: the experimental data
    % Cm: theoretical data
    % iT: the initial bit delay from experiments (empirical value)
    % offset: the offset for the presented symbols (aestetics)
    % min: low range of the x-axis ticks
    % mid: the increments for the x-axis ticks
    % max: the max value for the x-axis ticks
    % xmax: the max value for the plot x-axis
    % ymax: the max value for the plot y-axis
    % -------------------------------------------------------------------------
    
    figure
    
    plot(expData*10^9,'Linewidth',2)
    hold on
    plot(Cm,'Linewidth',2)
    
    symbol  = length(Cm)/Tt;
    bitData = zeros(1,symbol);
    
    message = [0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
    k       = [0 1 1 2 2 1 1 1 1 4 2 1 1 1 2 3];
    
    op = 1;
    
    % This loop is to extract the transmitted symbols and arrange them based on
    % the modulation level
    
    for j = 2:1:length(message)
        while k(j) > 0
            if message(j) == 0
                bitData(op) = 0; op = op + 1; k(j) = k(j) - 1;
            
            elseif message(j) == 1
                bitData(op) = 1; op = op + 1; k(j) = k(j) - 1;
            end
        end
    end
    
    stringBitData = string(bitData);
    
    % These swicth cases allows to print the symbols correctly for each
    % modulation level
    
    switch mod
        case 2
            for k = 0:1:symbol-1
                text(iT + k*Tt, offset, stringBitData(k+1), ...
                    'Fontsize', 16, 'Interpreter', 'Latex')
                drawnow;
            end
    
        case 4
            for k = 0:1:symbol-1
                text(iT + k*Tt, offset, append(stringBitData(2*k+1),  ...
                                               stringBitData(2*k+2)), ...
                                               'FontSize', 16,        ...
                                               'Interpreter', 'Latex')
                drawnow;
            end
    
        case 8
            for k = 0:1:symbol-1
                text(iT + k*Tt, offset, append(stringBitData(2*k+1),  ...
                                               stringBitData(2*k+2),  ...
                                               stringBitData(2*k+3)), ...
                                               'FontSize', 16,        ... 
                                               'Interpreter', 'Latex')
                drawnow;
            end
    end
    
    % Nice plot code ----------------------------------------------------------
    
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca, ...
      'Fontsize'        ,   18          , ...
      'TickDir'         ,   'out'       , ...
      'YGrid'           ,   'on'        , ...
      'GridLineStyle'   ,   '--'        , ...
      'TickLength'      ,   [.02 .02]   , ...
      'XTick'           ,   min:mid:max , ...
      'Box'             ,   'off'       , ...
      'XGrid'           ,   'on'        , ...
      'LineWidth'       ,   1.2         );
    
    % ------------------------------------------------------------------------- 
    
    xticks(0:Tt:length(Cm))
    
    ylabel('Signal Current [nA]','Interpreter','Latex')
    xlabel('Time [s]','Interpreter','Latex')
    
    xlim([0 xmax]); ylim([0 ymax])
    
end

% ======================== END OF CODE ====================================