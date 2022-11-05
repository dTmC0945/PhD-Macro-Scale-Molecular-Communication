function plotDataChannelSER(Error,I_XY, Channel, bitCounter, min, mid, max)

% -------------------------------------------------------------------------
% First Plot (Symbol Error Rate) ------------------------------------------
% -------------------------------------------------------------------------

dB_axis = (min:mid:max);                        % dB values for the x-axis 

semilogy(dB_axis,Error,'Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'      , ...
  'LineWidth'   , 1.5         );

% Label information -------------------------------------------------------

ylabel('Symbol Error Rate [SER]','Interpreter','Latex')
xlabel('SNR [dB]','Interpreter','Latex')
ylim([1e-4,1e-0]);

% Legend information ------------------------------------------------------

modLegend = cell(1,bitCounter + 1);

for i = 1:1:bitCounter
    modLegend{i} = join(compose("M = %d", 2^ i));

    if i == bitCounter
        modLegend{i+1} = "AWGN Channel";
    end
end

set(legend(modLegend{1:bitCounter}),...
    'Interpreter','Latex','Location','Southwest')

% -------------------------------------------------------------------------
% Second Plot (Mutual Information) ----------------------------------------
% -------------------------------------------------------------------------

figure

for i = 1:1:bitCounter
    plot(dB_axis,I_XY(:,i),'Linewidth',2)
    hold on
end

plot(dB_axis,Channel,'--k','Linewidth',2)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'     , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'Box'         , 'on'     , ...
  'LineWidth'   , 1.5         );

% Label information -------------------------------------------------------

ylabel('Mutual Information [bit/sym]','Interpreter','Latex')
xlabel('SNR [dB]','Interpreter','Latex')

set(legend(modLegend),...
    'Interpreter','Latex','Location','Northwest')

end