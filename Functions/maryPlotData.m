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