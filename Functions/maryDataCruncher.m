function Output = maryDataCruncher(filename)

% This function reads the experimental data saved in an excel function.

A = xlsread(filename);

M2_ary_30_1 = A(1,:); M2_ary_30_2 = A(2,:); M2_ary_30_3 = A(3,:);
M4_ary_30_1 = A(5,:); M4_ary_30_2 = A(6,:); M4_ary_30_3 = A(7,:);
M8_ary_30_1 = A(9,:); M8_ary_30_2 = A(10,:); M8_ary_30_3 = A(11,:);

M2_ary_60_1 = A(14,:); M2_ary_60_2 = A(15,:); M2_ary_60_3 = A(16,:);
M4_ary_60_1 = A(18,:); M4_ary_60_2 = A(19,:); M4_ary_60_3 = A(20,:);
M8_ary_60_1 = A(22,:); M8_ary_60_2 = A(23,:); M8_ary_60_3 = A(24,:);


M2_ary_90_1 = A(27,:); M2_ary_90_2 = A(28,:); M2_ary_90_3 = A(29,:);
M4_ary_90_1 = A(31,:); M4_ary_90_2 = A(32,:); M4_ary_90_3 = A(33,:);
M8_ary_90_1 = A(35,:); M8_ary_90_2 = A(36,:); M8_ary_90_3 = A(37,:);

M230 = (M2_ary_30_1 + M2_ary_30_2 + M2_ary_30_3)./3;
M430 = (M4_ary_30_1 + M4_ary_30_2 + M4_ary_30_3)./3;
M830 = (M8_ary_30_1 + M8_ary_30_2 + M8_ary_30_3)./3;

M260 = (M2_ary_60_1 + M2_ary_60_2 + M2_ary_60_3)./3;
M460 = (M4_ary_60_1 + M4_ary_60_2 + M4_ary_60_3)./3;
M860 = (M8_ary_60_1 + M8_ary_60_2 + M8_ary_60_3)./3;

M290 = (M2_ary_90_1 + M2_ary_90_2 + M2_ary_90_3)./3;
M490 = (M4_ary_90_1 + M4_ary_90_2 + M4_ary_90_3)./3;
M890 = (M8_ary_90_1 + M8_ary_90_2 + M8_ary_90_3)./3;

Output = [M230;M430;M830;
          M260;M460;M860;
          M290;M490;M890];

end