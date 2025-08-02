clear all
clc

% Closed Distance Transmission ============================================

% Daniel T. McGuiness =====================================================

% XLSREAD -----------------------------------------------------------------

A = xlsread('CLOSED_TEST_FILE_Q.xlsx');
N = xlsread('Background (30 mins)_Q.xlsx');

N_adj = N(1:600);   % Adjusting the time frame of the noise data

% -------------------------------------------------------------------------

x_d = [0.5 1 1.5 2 2.5];    % Transmission distance x-axis values

A_00_1 = A(1,2:601);        % 0 m rep. 1
A_00_2 = A(2,2:601);        % 0 m rep. 1
A_00_3 = A(3,2:601);        % 0 m rep. 1

A_05_1 = A(4,2:601);        % 0.5 m rep. 1
A_05_2 = A(5,2:601);        % 0.5 m rep. 1
A_05_3 = A(6,2:601);        % 0.5 m rep. 1

A_10_1 = A(7,2:601);        % 1 m rep. 1
A_10_2 = A(8,2:601);        % 1 m rep. 1
A_10_3 = A(9,2:601);        % 1 m rep. 1

A_10_1(191) = (1*A_10_1(189) + 2*A_10_1(192))./3;
A_10_1(190) = (2*A_10_1(189) + 1*A_10_1(192))./3;

A_15_1 = A(10,2:601);        % 1.5 m rep. 1
A_15_2 = A(11,2:601);        % 1.5 m rep. 1
A_15_3 = A(12,2:601);        % 1.5 m rep. 1

A_15_1(525) = (1*A_15_1(524) + 1*A_15_1(526))./2;

A_20_1 = A(13,2:601);        % 2.0 m rep. 1
A_20_2 = A(14,2:601);        % 2.0 m rep. 1
A_20_3 = A(15,2:601);        % 2.0 m rep. 1

A_25_1 = A(16,2:601);        % 2.5 m rep. 1
A_25_2 = A(17,2:601);        % 2.5 m rep. 1
A_25_3 = A(18,2:601);        % 2.5 m rep. 1

%find(all(isnan(X),1))

A_25_1(206) = (1*A_25_1(205) + 1*A_25_1(207))./2;
A_25_1(266) = (1*A_25_1(265) + 1*A_25_1(267))./2;
A_25_1(312) = (1*A_25_1(311) + 1*A_25_1(313))./2;
A_25_1(406) = (1*A_25_1(405) + 1*A_25_1(407))./2;
A_25_1(512) = (1*A_25_1(511) + 1*A_25_1(513))./2;
A_25_1(536) = (1*A_25_1(535) + 1*A_25_1(537))./2;

A_30_1 = A(19,2:601);        % 3 m rep. 1
A_30_2 = A(20,2:601);        % 3 m rep. 1
A_30_3 = A(21,2:601);        % 3 m rep. 1

A_00 = (A_00_1(1:600) + A_00_2(1:600) + A_00_3(1:600))./3;
A_05 = (A_05_1(1:600) + A_05_2(1:600) + A_05_3(1:600))./3;
A_10 = (A_10_1(1:600) + A_10_2(1:600) + A_10_3(1:600))./3;
A_15 = (A_15_1(1:600) + A_15_2(1:600) + A_15_3(1:600))./3;
A_20 = (A_20_1(1:600) + A_20_2(1:600) + A_20_3(1:600))./3;
A_25 = (A_25_1(1:600) + A_25_2(1:600) + A_25_3(1:600))./3;
A_30 = (A_30_1(1:600) + A_30_2(1:600) + A_30_3(1:600))./3;

A_10(191) = (1*A_10(189) + 2*A_10(192))./3;
A_10(190) = (2*A_10(189) + 1*A_10(192))./3;

A_15(525) = (1*A_15(524) + 1*A_15(526))./2;

A_25(206) = (1*A_25(205) + 1*A_25(207))./2;
A_25(266) = (1*A_25(265) + 1*A_25(267))./2;
A_25(312) = (1*A_25(311) + 1*A_25(313))./2;
A_25(406) = (1*A_25(405) + 1*A_25(407))./2;
A_25(512) = (1*A_25(511) + 1*A_25(513))./2;
A_25(536) = (1*A_25(535) + 1*A_25(537))./2;

A_30(560) = (1*A_10(559) + 1*A_10(561))./2;

CT_M = cbrewer('qual', 'Set1', 10);


figure % Amplitude signal plot---------------------------------------------

plot(A_05*10^9,'Linewidth',2,  'Color',    CT_M(1,:))
hold on
plot(A_10*10^9,'Linewidth',2,  'Color',    CT_M(2,:))
hold on
plot(A_15*10^9,'Linewidth',2,  'Color',    CT_M(3,:))
hold on
plot(A_20*10^9,'Linewidth',2,  'Color',    CT_M(4,:))
hold on
plot(A_25*10^9,'Linewidth',2,  'Color',    CT_M(5,:))
hold on
plot(A_30*10^9,'Linewidth',2,  'Color',    CT_M(6,:))

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

L = legend('$x$ = 0.5 m', ...                          % Legend data
           '$x$ = 1.0 m', ...
           '$x$ = 1.5 m', ...
           '$x$ = 2.0 m', ...
           '$x$ = 2.5 m', ...
           '$x$ = 3.0 m');                             
set(L,'Interpreter','Latex','Location','Northeast');    % Latex Legend

xlabel('Time [s]','Interpreter','Latex','Fontsize'    , 15)                % x-axis
ylabel('Signal Current [nA]','Interpreter','Latex','Fontsize'    , 15)     % y-axis

% Signal Amplitude Plot ---------------------------------------------------

An = [max(A_00_1) max(A_00_2) max(A_00_3);
      max(A_05_1) max(A_05_2) max(A_05_3);
      max(A_10_1) max(A_10_2) max(A_10_3);
      max(A_15_1) max(A_15_2) max(A_15_3);
      max(A_20_1) max(A_20_2) max(A_20_3);
      max(A_25_1) max(A_25_2) max(A_25_3);
      max(A_30_1)  max(A_30_2) max(A_30_3)];

Anavg = [max(A_00) max(A_05) max(A_10) max(A_15) max(A_20) max(A_25) max(A_30)];
  
An_05percent = 0.05.*Anavg.*10^9;
An_95percent = 0.95.*Anavg.*10^9;


time = [90 103 120 135 149 156] - 60; 

An = An(1:7,:);

AJ = [sum(A_05) sum(A_10) sum(A_15) sum(A_20) sum(A_25)];

Energy = [sum(abs(A_05_1*10^9).^2) sum(abs(A_05_2*10^9).^2) ...
          sum(abs(A_05_3*10^9).^2);
          sum(abs(A_10_1*10^9).^2) sum(abs(A_10_2*10^9).^2) ...
          sum(abs(A_10_3*10^9).^2);
          sum(abs(A_15_1*10^9).^2) sum(abs(A_15_2*10^9).^2) ...
          sum(abs(A_15_3*10^9).^2);
          sum(abs(A_20_1*10^9).^2) sum(abs(A_20_2*10^9).^2) ...
          sum(abs(A_20_3*10^9).^2);
          sum(abs(A_25_1*10^9).^2) sum(abs(A_25_2*10^9).^2) ...
          sum(abs(A_25_3*10^9).^2);sum(abs(A_30_1*10^9).^2) ...
          sum(abs(A_30_2*10^9).^2) sum(abs(A_30_3*10^9).^2)];

E = Energy(:,2);

range_SNR05 = range([snr(A_05_1,N_adj) snr(A_05_2,N_adj) snr(A_05_3,N_adj)]);
range_SNR10 = range([snr(A_10_1,N_adj) snr(A_10_2,N_adj) snr(A_10_3,N_adj)]);
range_SNR15 = range([snr(A_15_1,N_adj) snr(A_15_2,N_adj) snr(A_15_3,N_adj)]);
range_SNR20 = range([snr(A_20_1,N_adj) snr(A_20_2,N_adj) snr(A_20_3,N_adj)]);
range_SNR25 = range([snr(A_25_1,N_adj) snr(A_25_2,N_adj) snr(A_25_3,N_adj)]);

err_SNR = [range_SNR05/2 range_SNR10/2  range_SNR15/2 ...
           range_SNR20/2  range_SNR25/2];

SNR = [([snr(A_05_1,N_adj) snr(A_05_2,N_adj) snr(A_05_3,N_adj)]);
       ([snr(A_10_1,N_adj) snr(A_10_2,N_adj) snr(A_10_3,N_adj)]);
       ([snr(A_15_1,N_adj) snr(A_15_2,N_adj) snr(A_15_3,N_adj)]);
       ([snr(A_20_1,N_adj) snr(A_20_2,N_adj) snr(A_20_3,N_adj)]);
       ([snr(A_25_1,N_adj) snr(A_25_2,N_adj) snr(A_25_3,N_adj)]);
       ([snr(A_30_1,N_adj) snr(A_30_2,N_adj) snr(A_30_3,N_adj)])];

% Modelling of The Transmission ===========================================

% Parameters --------------------------------------------------------------

radius = [0.8 1 1.2];

for l = 1:1:length(radius)

m = 0.29;

v = 1.25*1e-5/(2*pi*(radius(l)/100)^2);
y = 0;
r = radius(l)*10^-2;

mu = 1.66*10^-5;
%mu =  1.81*10^-5
rho = 1.165;
%rho = 1.56;
Re = rho*v*r*2/mu;

f_D = 64/Re;

tau = 1/8*f_D*rho*v^2;

u_shear = sqrt(tau/rho);

r = 100*r;
u_shear = 100*u_shear;

v = 100*v;

d = 5.93*2*r*u_shear;

%d = d + 1/48*(r^2*v*2);

for op = 1:1:500
    
    x = op*1 - 1;
    
    for t = 1:1:4000
        theta(op,t) = m - (m*erf(10^100*r/(d*sqrt(1/(d*t))*t))^2*(erf((x+t*v)/(2*d*sqrt(1/(d*t))*t))+erf((x-t*v)/(2*d*sqrt(1/(d*t))*t))))/2;
        
    end
        
end

for op = 1:1:500
    
    for t = 1:1:2000
        t_1(op,t) = ind2sub(theta(op,t),theta(op,t) >= 1/200*m);
    end
    
    t_spot(op) = find(t_1(op,:),1,'first');
    
    tm(l,op) = t_spot(op);
    
    Imax(l,op) = theta(op,t_spot(op) + 60);    
end

for op = 1:1:500

    x = op*1 - 1 ;
    
    for t = 1:1:2000
        theta0(op,t) = (theta(op,t_spot(op) + 60)*erf(10^100*r/(d*sqrt(1/(d*t))*t))^2*(erf((x+t*v)/(2*d*sqrt(1/(d*t))*t))+erf((x-t*v)/(2*d*sqrt(1/(d*t))*t))))/2;
        
    end
    
end

for op = 1:1:500
    
    for t = 1:1:2000
        t_0(op,t) = ind2sub(theta0(op,t),theta(op,t_spot(op) + 60) ~= theta0(op,t));
    end
    
    t0_spot(op) = find(t_0(op,:),1,'first');
end


for i = 1:1:500
   
    S(i,:) = [theta(i,1:t_spot(i) + 60) theta0(i,t_spot(i):end)];
    
    En(l,i) =  sum(abs(S(i,1:600)).^2);
    
    Snr(l,i) = snr(S(i,1:600),N_adj*10^9);
    
end


end
x_model = 0.01:0.01:5;

x_data  = (0.5:0.5:3) + 0.025;


S2 =[theta(50,1:t_spot(50) + 60) theta0(50,t_spot(50):end)];
S3 =[theta(100,1:t_spot(100) + 60) theta0(100,t_spot(100):end)];
S4 =[theta(150,1:t_spot(150) + 60) theta0(150,t_spot(150):end)];
S5 =[theta(200,1:t_spot(200) + 60) theta0(200,t_spot(200):end)];
S6 =[theta(250,1:t_spot(250) + 60) theta0(250,t_spot(250):end)];
S7 =[theta(300,1:t_spot(300) + 60) theta0(300,t_spot(300):end)];

figure
plot(S2,'Linewidth',2)
hold on
plot(A_05(75:end)*10^9,'--','Linewidth',2)
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

xlim([0 600])
ylim([0 0.4])
xlabel('Time [s]','Interpreter','Latex','Fontsize'    , 16)                % x-axis
ylabel('Signal Current [nA]','Interpreter','Latex','Fontsize'    , 16)     % y-axis

S_label = legend('Experimental Results','Theoretical Model');
set(S_label,'Interpreter','Latex','Orientation','Vertical')


figure
plot(S3,'Linewidth',2)
hold on
plot(A_10(75:end)*10^9,'--','Linewidth',2)
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

xlabel('Time [s]','Interpreter','Latex','Fontsize'    , 16)              % x-axis
ylabel('Signal Current [nA]','Interpreter','Latex','Fontsize'    , 16)     % y-axis

S_label = legend('Experimental Results','Theoretical Model');
set(S_label,'Interpreter','Latex','Orientation','Vertical')

xlim([0 600])
ylim([0 0.4])

figure
plot(S4,'Linewidth',2)
hold on
plot(A_15(75:end)*10^9,'--','Linewidth',2)
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

xlabel('Time [s]','Interpreter','Latex','Fontsize'    , 16)              % x-axis
ylabel('Signal Current [nA]','Interpreter','Latex','Fontsize'    , 16)     % y-axis

S_label = legend('Experimental Results','Theoretical Model');
set(S_label,'Interpreter','Latex','Orientation','Vertical')

xlim([0 600])
ylim([0 0.4])

figure
plot(S5,'Linewidth',2)
hold on
plot(A_20(75:end)*10^9,'--','Linewidth',2)
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

xlabel('Time [s]','Interpreter','Latex','Fontsize'    , 16)              % x-axis
ylabel('Signal Current [nA]','Interpreter','Latex','Fontsize'    , 16)     % y-axis

S_label = legend('Experimental Results','Theoretical Model');
set(S_label,'Interpreter','Latex','Orientation','Vertical')

xlim([0 600])
ylim([0 0.4])

figure
plot(S6,'Linewidth',2)
hold on
plot(A_25(75:end)*10^9,'--','Linewidth',2)
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

S_label = legend('Experimental Results','Theoretical Model');
set(S_label,'Interpreter','Latex','Orientation','Vertical')


xlim([0 600])
ylim([0 0.4])

xlabel('Time [s]','Interpreter','Latex','Fontsize'    , 16)              % x-axis
ylabel('Signal Current [nA]','Interpreter','Latex','Fontsize'    , 16)     % y-axis



figure
plot(A_30(75:end)*10^9,'Linewidth',2)
hold on
plot(S7,'--','Linewidth',2)
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

xlabel('Time [s]','Interpreter','Latex','Fontsize'    , 16)              % x-axis
ylabel('Signal Current [nA]','Interpreter','Latex','Fontsize'    , 16)     % y-axis

S_label = legend('Experimental Results','Theoretical Model');
set(S_label,'Interpreter','Latex','Orientation','Vertical')


xlim([0 600])
ylim([0 0.4])


A7 = S7(1:525);
B7 = A_30(76:600)*10^9;
Co7 = corrcoef(A7,B7);

A6 = S6(1:525);
B6 = A_25(76:600)*10^9;
Co6 = corrcoef(A6,B6);

A5 = S5(1:525);
B5 = A_20(76:600)*10^9;
Co5 = corrcoef(A5,B5);

A4 = S4(1:525);
B4 = A_15(76:600)*10^9;
Co4 = corrcoef(A4,B4);

A3 = S3(1:525);
B3 = A_10(76:600)*10^9;
Co3 = corrcoef(A3,B3);

A2 = S2(1:525);
B2 = A_05(76:600)*10^9;
Co2 = corrcoef(A2,B2);

mean_en = [mean(Energy(1,:)), ...
    mean(Energy(2,:)), ...
    mean(Energy(3,:)), ...
    mean(Energy(4,:)), ...
    mean(Energy(5,:)), ...
    mean(Energy(6,:))];

range_en = [range(Energy(1,:)), ...
    range(Energy(2,:)), ...
    range(Energy(3,:)), ...
    range(Energy(4,:)), ...
    range(Energy(5,:)), ...
    range(Energy(6,:))];

E = [sum(abs(S2.^2)) sum(abs(S3.^2)) sum(abs(S4.^2)) sum(abs(S5.^2)) sum(abs(S6.^2)) sum(abs(S7.^2))];

corrcoef(mean_en,E)


figure
semilogx(x_model,En(1,:),'Linewidth',2,  'Color',    CT_M(1,:))
hold on
semilogx(x_model,En(2,:),'Linewidth',2,  'Color',    CT_M(2,:))
hold on
semilogx(x_model,En(3,:),'Linewidth',2,  'Color',    CT_M(3,:))
hold on
he1=errorbar(x_data,mean_en,range_en/2,'o','Linewidth',1.2,'MarkerFaceColor',CT_M(2,:),'MarkerEdgeColor',CT_M(2,:))
set(he1                           , ...
  'LineWidth'       , 1.2           , ...
  'Marker'          , 'o'         , ...
  'Color'           , 'k'           ,...
  'MarkerSize'      , 8           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , CT_M(2,:)  );


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
  'Yscale'      ,'log'      , ...
  'Xscale'      ,'log'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Distance [m]','Interpreter','Latex')
ylabel('Signal Energy [nW]','Interpreter','Latex')

xlim([0 5])
ylim([0.01 10])
En_label = legend('$R$ = 0.8 cm','$R$ = 1 cm','$R$ = 1.2 cm','Experimental Results');
set(En_label,'Interpreter','Latex','Location','Southwest')



mean_snr = [mean(SNR(1,:)), ...
    mean(SNR(2,:)), ...
    mean(SNR(3,:)), ...
    mean(SNR(4,:)), ...
    mean(SNR(5,:)), ...
    mean(SNR(6,:))];

range_snr = [range(SNR(1,:)), ...
    range(SNR(2,:)), ...
    range(SNR(3,:)), ...
    range(SNR(4,:)), ...
    range(SNR(5,:)), ...
    range(SNR(6,:))];

corrcoef(mean_snr,10*log10(E./0.096));


figure

semilogx(x_model,Snr(1,:),'Linewidth',2,  'Color',    CT_M(1,:))
hold on
semilogx(x_model,Snr(2,:),'Linewidth',2,  'Color',    CT_M(2,:))
hold on
semilogx(x_model,Snr(3,:),'Linewidth',2,  'Color',    CT_M(3,:))
hold on
he1=errorbar(x_data,mean_snr,range_snr/2,'o','Linewidth',1.2,'MarkerFaceColor',CT_M(2,:),'MarkerEdgeColor',CT_M(2,:))
set(he1                           , ...
  'LineWidth'       , 1.2           , ...
  'Marker'          , 'o'         , ...
  'Color'           , 'k'           ,...
  'MarkerSize'      , 8           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , CT_M(2,:)  );

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
      'Xscale'      ,'log'      , ...
  'Box'         , 'off'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Distance [m]','Interpreter','Latex')
ylabel('Signal-to-Noise Ratio [dB]','Interpreter','Latex')

xlim([0 5])

Snr_label = legend('$R$ = 0.8 cm','$R$ = 1 cm','$R$ = 1.2 cm','Experimental Results');
set(Snr_label,'Interpreter','Latex','Location','Southwest')


time_delay_T = [tm(1,50) + 20 tm(1,100) + 20 tm(1,150) + 20 tm(1,200) + 20 tm(1,250) + 20 tm(1,300) + 20];
time_delay_E = time;

T_corr = corrcoef(time_delay_E, time_delay_T)

figure
semilogx(x_model,tm(1,:) + 20,'Linewidth',2,  'Color',    CT_M(1,:))
hold on
semilogx(x_model,tm(2,:) + 20,'Linewidth',2,  'Color',    CT_M(2,:))
hold on
semilogx(x_model,tm(3,:) + 20,'Linewidth',2,  'Color',    CT_M(3,:))
hold on
he1=semilogx(x_data,time,'o','Linewidth',1.2,'MarkerFaceColor',CT_M(2,:),'MarkerEdgeColor',CT_M(2,:))
set(he1                           , ...
  'LineWidth'       , 1.2           , ...
  'Marker'          , 'o'         , ...
  'Color'           , 'k'           ,...
  'MarkerSize'      , 8           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , CT_M(2,:)  );

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
    'Xscale'      ,'log'      , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Distance [m]','Interpreter','Latex')
ylabel('Transmission Delay [s]','Interpreter','Latex')

xlim([0 5])

Snr_label = legend('$R$ = 0.8 cm','$R$ = 1 cm','$R$ = 1.2 cm','Experimental Results');
set(Snr_label,'Interpreter','Latex','Location','Northwest')


% Amplitude Plot ----------------------------------------------------------

mean_amp = [mean(An(2,:)*10^9), ...
    mean(An(3,:)*10^9), ...
    mean(An(4,:)*10^9), ...
    mean(An(5,:)*10^9), ...
    mean(An(6,:)*10^9), ...
    mean(An(7,:)*10^9)];

range_amp = [range(An(2,:)*10^9), ...
    range(An(3,:)*10^9), ...
    range(An(4,:)*10^9), ...
    range(An(5,:)*10^9), ...
    range(An(6,:)*10^9), ...
    range(An(7,:)*10^9)];

figure

amp_T = [Imax(1,50) Imax(1,100) Imax(1,150) Imax(1,200) Imax(1,250) Imax(1,300)];

corrcoef(amp_T,mean_amp)

semilogx(x_model,Imax(1,:),'Linewidth',2,  'Color',    CT_M(1,:))
hold on
semilogx(x_model,Imax(2,:),'Linewidth',2,  'Color',    CT_M(2,:))
hold on
semilogx(x_model,Imax(3,:),'Linewidth',2,  'Color',    CT_M(3,:))
hold on
he1=errorbar(x_data,mean_amp,range_amp/2,'o','Linewidth',1.2,'MarkerFaceColor',CT_M(2,:),'MarkerEdgeColor',CT_M(2,:))
set(he1                           , ...
  'LineWidth'       , 1.2           , ...
  'Marker'          , 'o'         , ...
  'Color'           , 'k'           ,...
  'MarkerSize'      , 8           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , CT_M(2,:)  );

Amp_label = legend('$R$ = 0.8 cm','$R$ = 1 cm','$R$ = 1.2 cm','Experimental Results');
set(Amp_label,'Interpreter','Latex','Location','Southwest')

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
  'Box'         , 'off'     , ...
  'XScale'      , 'log'     , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

xlabel('Distance [m]','Interpreter','Latex')
ylabel('Signal Amplitude [nA]','Interpreter','Latex')

xlim([0 5])

