clear all
clc

filename = 'ber25sec (Quantised).xlsx';

A = xlsread(filename);

i = 1;

m = 1;

T = zeros(1,200);
B = zeros(1,200);
C = zeros(1,20);
D = zeros(1,20);
error_array = zeros(1,20);
error_array01 = zeros(1,20);
error_array10 = zeros(1,20);
H_b_e0 = zeros(1,20);
H_b_e1 = zeros(1,20);

figure
time = 0:1:5056;
plot(time,A*10^9,'Linewidth',2)
ylabel('Signal Current (nA)','Interpreter','latex')
xlabel('Transmission Time (sec)','Interpreter','latex')
title('Bit-Error-Rate Transmission of Acetone ($T_{bit} = 25s$)','Interpreter','latex')

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XTick'       , 0:1000:6000   , ...
  'YTick'       , 0:5:60 , ...
  'LineWidth'   , 1.5         );

xlim([0 5056])

while i <= 200
    
    if i == 1
    
        B(m) = A(45);
        
        T(m) = 45;
    
    else
        
        B(m) = A(20 + 25*(i));
        
        T(m) = 20 + 25*(i);
        
    end
        
    m = m + 1;
    
    i = i + 1;
        
end

n = 1;

while n <= length(B)
    
    if n == 1
        
        if B(n) < 1 * 10^-10
        
        D(n) = 0;
        
        end
        
    elseif B(n) > B(n-1)
            
        if abs(B(n-1) - B(n)) < 0.1 * B(n)
            
            D(n) = D(n-1);
        
        else
            
            D(n) = 1;
        
        end
    
    elseif B(n) < B(n-1)
        
        if abs(B(n-1) - B(n)) < 0.1 * B(n)
            
            D(n) = D(n-1);
        
        else
        
            D(n) =  0;
        
        end
        

    

    end
    
    n = n + 1;
    
end

        
data = [0	1	0	1	1	1	0	1	0	1	0	1	0	0	0	0	0	...
        0	1	1	0	0	1	0	1	1	0	0	1	1	0	0	1	0	...
        1	0	1	0	0	1	0	0	0	1	0	1	1	0	0	0	0	...
        0	1	1	1	1	1	0	0	0	1	1	1	1	0	0	0	1	...
        1	0	1	1	0	1	0	0	0	1	1	0	1	0	0	0	0	...
        1	0	0	1	0	1	0	0	0	1	1	1	0	0	1	0	1	...
        1	0	0	0	0	1	0	1	0	1	1	0	1	1	0	0	1	...
        1	0	1	0	0	1	1	0	0	0	1	1	1	1	1	0	0	...
        1	1	1	1	0	1	0	0	1	1	1	1	1	1	0	0	1	...
        0	0	1	0	0	1	0	1	1	1	0	0	0	1	0	1	1	...
        1	0	0	0	1	0	0	1	0	0	1	1	1	1	0	1	1	...
        0	0	1	0	1	1	1	1	1	0	1	1	0];
    
u = 1;

e = 0;

for k  = 1:1:20

    while u <= 10*k
    
    if D(u) == data(u)
        
        fprintf('%d is correct \n',u)
        
    else
        
        fprintf('%d is wrong \n',u)
        
        e = e + 1;
        
    end
    
    u = u + 1;
    
    end
    error_array(k) = e/(10*k);

end

u = 1;

e01     = 1e-11;
e10     = 1e-11;
e01m    = 1e-11;
e10m    = 1e-11;
c1 = 0;
c0 = 0;
for k  = 1:1:20

    while u <= 10*k
        
        if D(u) == 1
        
            c0 = c0 + 1;
        
        else
        
            c1 = c1 + 1;
        end
    
    if D(u) == 0 && data(u) == 1
        
        e01 = e01 + 1;
        e01m(k) = e01/(c0);
        
    elseif D(u) == 1 && data(u) == 0
        
        e10 = e10 + 1;
        e10m(k) = e10/c1;
        
    end
    
    u = u + 1;
    
    end
    error_array01(k) = e01/(c0);
    error_array10(k) = e10/(c1);
    
    H_b_e0(k) = -error_array01(k) * log2(error_array01(k)) - (1 - error_array01(k)) * log2(1 - error_array01(k));
    H_b_e1(k) = -error_array10(k) * log2(error_array10(k)) - (1 - error_array10(k)) * log2(1 - error_array10(k));

% Binary Asymmetric Channel (BAC) main equation

    C(k) = - (1 - error_array10(k)) / (1 - error_array01(k) - error_array10(k)) * H_b_e0(k) + ...
                (error_array01(k)) / (1 - error_array01(k) - error_array10(k)) * H_b_e1(k) + ...
                log2(1 + 2^((H_b_e0(k) - H_b_e1(k))/(1-error_array01(k)-error_array10(k))));
    
end

packet_array = 10:10:200;

figure
plot(packet_array,error_array01,'-o','Linewidth',2)
hold on
plot(packet_array,error_array10,'-o','Linewidth',2)
hold on
plot(packet_array,error_array,'-o','Linewidth',2)
hold on

ylabel('Error Probability ($\%$)','Interpreter','latex')
xlabel('Packet Size (bits)','Interpreter','latex')
title('Error Probability of Acetone ($T_{bit} = 25s$)','Interpreter','latex')

L = legend('$e_{P(1|0)}$','$e_{P(0|1)}$','$e_{P(1|0)}$ + $e_{P(0|1)}$');

set(L,'Interpreter','latex','Orientation','horizontal')

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XTick'       , 0:20:200   , ...
  'YTick'       , 0:0.1:0.5    , ...
  'LineWidth'   , 1.5         );

ylim([0 0.5])

figure

plot(packet_array,C,'-o','Linewidth',2)

ylabel('Channel Capacity (bits/sec)','Interpreter','latex')
xlabel('Packet Size (bits)','Interpreter','latex')
title('Channel Capacity of Acetone ($T_{bit} = 25s$)','Interpreter','latex')

set(L,'Interpreter','latex','Orientation','horizontal')

set(gca,'TickLabelInterpreter', 'latex');
set(gca, ...
  'Fontsize'    , 16        , ...
  'TickDir'     , 'in'      , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XTick'       , 0:20:200   , ...
  'YTick'       , 0.1:0.1:0.90 , ...
  'LineWidth'   , 1.5         );


m = 1;
j = 1;

for n = 1:1:200
    
    if D(n) == 0
        
        B_0(m) = B(n)
        
        m = m + 1;
        
    else
        
        B_1(j) = B(n)
        
        j = j + 1;
        
    end
end
