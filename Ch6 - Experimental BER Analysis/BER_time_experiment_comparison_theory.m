% MOLECULAR SINGLE VELOCITY ANALYSIS -----------------------------------------
% Daniel McGuiness 06.07.2018 ---------------------------------------------
% -------------------------------------------------------------------------

clear all
clc

% Parameters ==============================================================

x = 2.5;                                    % Transmission Distance(cm)
v = 0.10;                                    % Advective Flow       (cm/s)
r = 0.3;                                    % Detector Radius      (cm)

D_L = 0.124;                                  % Longitudional coefficient of 
                                            % Diffusivity          (cm^2/s)

M = 1;                                      % Injected Mass        (ng)
D_R = 0.001;                                % Radial coefficient of 
                                            % Diffusivity          (cm^2/s)

%x = 2.5;                                % Transmission Distance (cm)
%d = 0.124;                              % Advective Flow (cm)
%M = 0.9;                                % Injected mass (kg)
%v = 0.135;
bit_length = 10000;                    % Transmitted Random bit length

CBER_Error  = [26 22 28 29; 20 17 11 10 ; 12 1 4 1; 11 5 4 0; 2 3 2 0]./100;


for bit_counter = 1
    
    level_counter = 2^(bit_counter);    % Number of levels in the channel
    
    % Simulation Parameters
    
    counter = 1;
    
for Tt = [5 10 15 20 25]                 % Counter operator (dB)
    dB = 25;
    SNR = 10^(dB/10);
    
    
    H_YX            = zeros(level_counter,level_counter); % Contional
    H_Y             = zeros(1,level_counter);               % Ent. (Y)
    Enput           = zeros(level_counter,level_counter); % Probability 
    Error_Matrix    = zeros(level_counter,level_counter); % Error
    Decoded         = zeros(1,bit_length);                    % Decoded
    Cecoded         = zeros(1,bit_length);                    % Cecoded
    
    Adj_Probability_Matrix = zeros(level_counter,level_counter);
    

    % Cleainng the input matrix and the probability matrix before each
    % iteration
    
    Probability_Matrix  = zeros(level_counter,level_counter);
    Input_Matrix        = zeros(1,level_counter);
    
    %Tt      = 60;                  % Bit duration (counter) (sec)
    
    Cm1     = 0;                        % First value of transmission
                                        % (keep it 0)

    M_old   = 0;                        % Leftover chemicals in the system  
                                        % (keep it 0)
    Ts      = 0;                        % Transmission time                 
                                        % (keep it 0)
    Q_neg   = 0;                        % Flushed chemicals                 
                                        % (keep it 0)

    p       = randi([0 level_counter-1],1,bit_length);  % Random number generator
    
    % Randomize the generated number for both message1 (Prime_seed_1) and

    Prime_seed_1 = p(randperm(length(p)));
    
    % Assigning the generated values into the message parameters
    
    message1    = Prime_seed_1;

    % Creating a D parameter for use in for loops
    
    D           = message1;

    Re1 = [];
    Pe1 = [];
    
    i = 0;
    
    w = 1;
    
    for n = 1:1:length(message1)

        if n == 1

            i = i + 1;

            Re1(w) = message1(n);

            Pe1(w) = i;

        elseif message1(n-1) == message1(n)

            i = i + 1;

            Pe1(w) = i;

            Re1(w) = message1(n);

        elseif message1(n-1) ~= message1(n)

            Pe1(w) = i;

            w = w + 1;

            Re1(w) = message1(n);

            i = 1;

            if n == length(message1)

                Pe1(w) = 1;

                w = w + 1;

             end
        end
    end

    sum_P1 = sum(Pe1);      % Sum check (must be equal to bit_length)

    message1 = [0 Re1];     % reassign message1 (must add 0)

    k1  = [0 Pe1];          % Reassign bit frequency for message1

    for j = 2:1:length(message1)

            if message1(j) > message1(j-1)

                M_old = Cm1(end);

                for t = 1:1:k1(j)*Tt

                    Q(t) = M_old + (message1(j) - message1(j-1))- ...
                                       (message1(j) - message1(j-1))* ...
                                       erf(r/sqrt(4*D_R*t))^2*...
                                      (erf((x - v*t)/(2*sqrt(D_L*t))) + ...
                                       erf((x + v*t)/(2*sqrt(D_L*t))))/2 ;

                    Cm1(t + Ts) = Q(t); 

                end

            else
                    M_r = Cm1(end); 

                for t = 1:1:k1(j)*Tt

                    Q(t) = (abs(M_r - message1(j)))* ...
                                  erf(r/sqrt(4*D_R*t))^2* ...
                               (erf((x - v*t)/(2*sqrt(D_L*t))) + ...
                                erf((x + v*t)/(2*sqrt(D_L*t))))/2 +...
                                message1(j);

                    Cm1(t + Ts) = Q(t); 

                end

            end

            Ts = k1(j)*Tt + Ts;

    end

    Cm1 = awgn(Cm1,dB,'measured','dB');     % Add AWGN to the transmitted 
                                                 % signal (message1)
                                              
    t = 1:1:length(Cm1);    % Transmission time x-axis value

    % Generating discrete received values from the transmission
    
    for i = 2:1:length(D)

               R1(i) = Cm1(Tt*i);   % message1 received values

    end
    
    % Generation of the symbol values for use in communications
    
    for a = 1:1:level_counter % Number of symbol values
        
        Alphabet(a) = a - 1; 
        
    end
    
    % Naming of the symbol for each constellation values
    
    alpha_counter = 1;
    
    for a = 1:1:level_counter % In-phase value
            
            Input_Alphabet(a) = Alphabet(alpha_counter);
            
            alpha_counter = alpha_counter + 1;
    end

    
    % Decoded Values
    
    for i = 2:1:length(D)
        
        for a = 1:1:level_counter
                
                if R1(i) < Input_Alphabet(1,1) - 0.5 
                    
                    Decoded(i) = 0;
                    
                elseif R1(i) > Input_Alphabet(1,level_counter) + 0.5
                    
                    Decoded(i) = level_counter - 1;
                    
                elseif  R1(i) >= (a - 1 - 1/2) ...
                     && R1(i) <= (a - 1 + 1/2)
                    Decoded(i) = Input_Alphabet(1,a);
                    
                end
                
        end

        
        for a = 1:1:level_counter

                if  Prime_seed_1(i) >= (a  - 1 ...
                        - 1/2) ...
                 && Prime_seed_1(i) <= (a  - 1 ...
                        + 1/2) 
                    
                        Cecoded(i) = Input_Alphabet(1,a);
                end
                
       
        end
        
    end

    for c = 1:1:length(Decoded)
        
        for a = 1:1:level_counter + 1
            
            if Cecoded(c) == (a - 1)
                
                Input_Matrix(a) = Input_Matrix(a) + 1; 
            end
            
            for b = 1:1:level_counter + 1
             
                if Decoded(c) == (a - 1) && Cecoded(c) == (b - 1)
                    
                    Probability_Matrix(a,b) = Probability_Matrix(a,b) + 1;

                end
            end
        end
    end
    
    % Check values for the probability matrix and the input matrix
    % Both of them must equal to 1
    
    Probability_check = sum(sum(Probability_Matrix))./length(D);
    Input_check       = sum(sum(Input_Matrix))      ./length(D);
    
    % Generation the amount of input a single output received.
    
    Probability_length = sum(Probability_Matrix,1);

    % Adjusting the probability matrix by dividing the probability matrix
    % with the probability length
    
    for a = 1:1:level_counter
            
        for b = 1:1:level_counter
            
            Adj_Probability_Matrix(a,b) = ...
                Probability_Matrix(a,b)./Probability_length(b);
            
        end
    end
        
    % Error Count ---------------------------------------------------------
    
    Error_Input_Matrix = (Adj_Probability_Matrix - eye(level_counter)...
                        .*Adj_Probability_Matrix);
    
    for a = 1:1:level_counter
            
        for b = 1:1:level_counter
            
            Error_Matrix(a,b) = ...
                Error_Input_Matrix(b,a)*Input_Matrix(a)/length(D) ;
            
        end
    end
               
    Error(counter) = sum(sum(Error_Matrix));
    
    % ---------------------------------------------------------------------
    
    for c = 1:1:length(Decoded)
        
        for a = 1:1:level_counter
            
            for b = 1:1:level_counter
                
                Enput(a,b) =  Input_Matrix(a)...
                    /length(D)*Adj_Probability_Matrix(b,a);
                         
            end
        end
    end
    
    Enput_check = sum(sum(Enput));      % P_Y Check the sum 
                                        % all probabilities equal to 1
                                        
    P_Y_adj = sum(Enput);
    
    for a = 1:1:level_counter
        
            H_Y(a) = - P_Y_adj(a)*mylog2(P_Y_adj(a));
            
    end
    
    Adj_Input = Input_Matrix./length(D);
    
    Adj_H_Y = sum(H_Y);
    
    for a = 1:1:level_counter
        for b = 1:1:level_counter
            
            H_YX(a,b) = Enput(a,b)*mylog2((Adj_Input(1,a))/Enput(a,b));

        end
    end
    
    H_YX(isnan(H_YX)) = 0;
    
    Adj_H_YX = sum(sum(H_YX));
    
    I_XY(counter) = Adj_H_Y - Adj_H_YX
        
    counter = counter + 1;
    
end
end

y_pos = [3.33333333333334,3.66666666666667,1,0.666666666666667,0]./100;
y_neg = [4.66666666666666,2.33333333333333,1,0.333333333333333,0]./100;

BERE =[38.6666666666667,19.3333333333333,6,4.33333333333333,0.333333333333333]./100;

dB_axis = [5 10 15 20 25];
exp_axis = [5 10 15 20 25];


hold on
he1 = errorbar(exp_axis, BERE,y_pos,y_neg,'o','Linewidth',1.2)
hold on
semilogy(dB_axis,Error,'s','Linewidth',2,...
                'MarkerFaceColor',[0.8500    0.3250    0.0980]   , ...
                'MarkerEdgeColor',[0.8500    0.3250    0.0980]   , ...
                'MarkerSize', 8)
            
set(he1                        , ...
  'LineWidth'       , 1.5      , ...
  'Color'           , [0    0.4470    0.7410]   , ...
  'MarkerSize'      , 8        , ...
  'MarkerEdgeColor' , [0    0.4470    0.7410]  , ...
  'MarkerFaceColor' , [0    0.4470    0.7410]   );

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
  'YScale'      , 'log'     , ...
  'LineWidth'   , 1.2         );

% -------------------------------------------------------------------------

ylim([1e-4,1e-0]);

ylabel('Symbol Error Rate [SER]','Interpreter','Latex')
xlabel('Symbol Duration [s]','Interpreter','Latex')
%title({'$D$ = 0.124 $\mathrm{cm^2/s}$, $u_x$ = 0.136 $\mathrm{cm/s}$'; '$x$ = 2.5 $\mathrm{cm}$, $\rho$ = 0.94'},'Interpreter','Latex','Fontsize',16)

label = legend('Experimental Results','Theoretical Model');
set(label,'Interpreter','Latex','Orientation','Vertical','Location','Southwest')


figure

plot(dB_axis,I_XY,'-s','Linewidth',2,...
                'MarkerFaceColor',[0    0.4470    0.7410]   , ...
                'MarkerEdgeColor',[0    0.4470    0.7410])
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

ylabel('Mutual Information [bit/sym]','Interpreter','Latex')
xlabel('Symbol Duration [s]','Interpreter','Latex')

BER_Error = [Error(1) Error(2) Error(3) Error(4) Error(5)];

corrcoef(BER_Error,BERE)