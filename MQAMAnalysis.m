% ========================= BEGIN CODE ====================================

% -------------------------------------------------------------------------
% -- MQAM Analysis by Daniel T. McGuiness ---------------------------------
% -------------------------------------------------------------------------

clear variables
clc

% Function Paths ----------------------------------------------------------

addpath("./Functions")

% Variables ---------------------------------------------------------------

distance        = 1;                % Transmission Distance (cm)
advectiveFlow   = 0.1;              % Diffusivity coefficient (cm^2/s)
mass            = 1;                % Injected mass (kg)
radialDiffusion = 1e-5;
symbolDuration  = 20;               % Bit duration (counter) (sec)
detectorRadius  = 0.03;             % Radius of the inlet (detector) (m)
Es_N0           = 25;               % Set value for energy per symbol

bitLength       = 1000;             % Transmitted Random bit length
varIterator     = 1;                % Variable iterator                
modLevelCounter = 2;                % Number of modulation levels of MQAM

% -------------------------------------------------------------------------
    

plotReal        = zeros(modLevelCounter^2, bitLength);
plotImag        = zeros(modLevelCounter^2, bitLength);

% -------------------------------------------------------------------------

for  longitudinalDiffusion = 0.1:0.1:0.4

    
    % Cleainng the input matrix and the probability matrix before each
    % iteration


    randomBits1 = 2*(randperm(length(...
                     randi(modLevelCounter,1,bitLength))));
    randomBits2 = 2*(randperm(length(...
                     randi(modLevelCounter,1,bitLength))));


    % Delegating the random bits to two distinct matrices -------------

    [bitSequence1, ...
     bitRepetition1]     = splitter(randomBits1);

    [bitSequence2, ...
     bitRepetition2]     = splitter(randomBits2);

    CmP = transmission2D(radialDiffusion, longitudinalDiffusion, ...
                         advectiveFlow, ...
                         distance, detectorRadius,...
                         bitSequence1, bitRepetition1,...
                         symbolDuration);
    
    CmQ = transmission2D(radialDiffusion, longitudinalDiffusion, ...
                         advectiveFlow, ...
                         distance, detectorRadius,...
                         bitSequence2, bitRepetition2,...
                         symbolDuration);

    signal2D = complex(CmP - (modLevelCounter + 1),...
                       CmQ - (modLevelCounter + 1));
    
    normal = sqrt(sum(abs(signal2D .^2)) ...
               / (symbolDuration*bitLength));

    noise = 10^(-Es_N0/20)/sqrt(2)...
          * complex(randn(1,bitLength*symbolDuration), ...
                    randn(1,bitLength*symbolDuration)); % white guassian noise, 0dB variance
    
    %signalNoised = addNoise(signal2D, 0, 1./noise);

    noisedSignal = signal2D/normal + noise; % additive white gaussian noise

    K = snr(noisedSignal,noise);
    
    T_Prob = decodingTheory2D(modLevelCounter, Es_N0);

    [inputMatrix, probabilityMatrix, dataPoints] = decodingSimulation2D(modLevelCounter, bitLength, symbolDuration, noisedSignal, normal, randomBits1, randomBits2);

    % Check values for the probability matrix and the input matrix
    % Both of them must equal to 1

    Theo_Probability_Matrix = T_Prob.*bitLength./modLevelCounter^2;
    checkProbability = sum(sum(probabilityMatrix))./bitLength;
    checkTheoryMatrix= sum(sum(Theo_Probability_Matrix))./bitLength;
    checkInputMatrix = sum(sum(inputMatrix))      ./bitLength;

    % Adjusting the probability matrix by dividing the probability matrix
    % with the probability length


    % Plot data for P and Q axis ------------------------------------------
    
    plotReal(varIterator,:) = real(dataPoints);    % P axis 
    plotImag(varIterator,:) = imag(dataPoints);    % Q axis
     
    varIterator = varIterator + 1;                  % variable loop counter
    
   
end

figure

scatter(plotReal, plotImag)


figure

subtightplot(3,3,[4 5 7 8],[0.01,0.01],0.15,0.15)

scatter(plotReal(1,:),plotImag(1,:),'filled','MarkerEdgeColor', [51 77 92]/255,'MarkerFaceColor', [51 77 92]/255)
hold on
scatter(plotReal(2,:),plotImag(2,:),'filled','MarkerEdgeColor', [69 178 157]/255,'MarkerFaceColor', [69 178 157]/255)
hold on
scatter(plotReal(3,:),plotImag(3,:),'filled','MarkerEdgeColor', [239 201 76]/255,'MarkerFaceColor', [239 201 76]/255)
hold on
scatter(plotReal(4,:),plotImag(4,:),'filled','MarkerEdgeColor', [226 122 63]/255,'MarkerFaceColor', [226 122 63]/255)

% Nice plot code ----------------------------------------------------------

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 14        , ...
  'TickDir'     , 'out'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'GridLineStyle','--'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XMinorGrid'  , 'on'      , ...
  'YMinorGrid'  , 'on'      , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.1         );

% -------------------------------------------------------------------------
xlabel('In-phase Chemical ($m_I$)','Interpreter','Latex')
ylabel('Out-phase Chemical ($m_Q$)','Interpreter','Latex')
xticks([1,2,3])
yticks([1,2,3])


MI_legend = legend('$D = 0.1\, \mathrm{cm^2/s}$',...
                   '$D = 0.2\, \mathrm{cm^2/s}$',...
                   '$D = 0.3\, \mathrm{cm^2/s}$',...
                   '$D = 0.4\, \mathrm{cm^2/s}$');
set(MI_legend,'Interpreter','Latex','box','off');
set(MI_legend,'position',[.71 .71 .1 .1]);
ylim([0.5 3.5])
xlim([0.5 3.5])




subtightplot(3,3,[1 2],[0.01,0.01],0.10,0.15)
nbins = 150;
h1 = histogram(plotReal(1,:),nbins, 'EdgeColor', 'none','FaceAlpha',1);
hold on
h2 = histogram(plotReal(2,:),nbins, 'EdgeColor', 'none','FaceAlpha',1);
hold on
h3 = histogram(plotReal(3,:),nbins, 'EdgeColor', 'none','FaceAlpha',1);
hold on
h4 = histogram(plotReal(4,:),nbins, 'EdgeColor', 'none','FaceAlpha',1);

xlim([0.5 3.5])

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   , ...
  'Fontsize'            , 14        , ...
  'TickDir'             , 'out'     , ...
  'YGrid'               , 'on'      , ...
  'XGrid'               , 'on'      , ...
  'GridLineStyle'       , '--'      , ...
  'TickLength'          , [.02 .02] , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'XMinorGrid'          , 'on'      , ...
  'YMinorGrid'          , 'on'      , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       , ...
  'xticklabel'          , {[]}      );

% -------------------------------------------------------------------------

subtightplot(3,3,[6 9],[0.01,0.01],0.15,0.10)

nbins = 150;
h5 = histogram(plotImag(1,:),nbins, ...
               'EdgeColor', 'none', ...
               'FaceAlpha', 1)
hold on
h6 = histogram(plotImag(2,:),nbins,'EdgeColor','none','FaceColor', [69 178 157]/255,'FaceAlpha',1)
hold on
h7 = histogram(plotImag(3,:),nbins,'EdgeColor', 'none','FaceColor', [239 201 76]/255,'FaceAlpha',1)
hold on
h8 = histogram(plotImag(4,:),nbins,'EdgeColor', 'none','FaceColor', [226 122 63]/255,'FaceAlpha',1)
xlim([0.5 3.5])

% Nice plot code ----------------------------------------------------------

set(gca, ...
  'TickLabelInterpreter', 'latex'   , ...
  'Fontsize'            , 14        , ...
  'TickDir'             , 'out'     , ...
  'TickLength'          , [.02 .02] , ...
  'GridLineStyle'       , '--'      , ...
  'XGrid'               , 'on'      , ...
  'YGrid'               , 'on'      , ...
  'XMinorGrid'          , 'on'      , ...
  'YMinorGrid'          , 'on'      , ...
  'XMinorTick'          , 'on'      , ...
  'YMinorTick'          , 'on'      , ...
  'YAxisLocation'       , 'right'   , ...
  'Box'                 , 'off'     , ...
  'LineWidth'           , 1.2       , ...
  'CameraUpVector'      , [1,0,0]   , ...
  'XDir'                , 'reverse' , ...
  'xticklabel'          , {[]}       );

% ======================== END OF CODE ====================================