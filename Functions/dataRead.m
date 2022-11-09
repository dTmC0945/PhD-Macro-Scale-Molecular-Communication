function A = dataRead(filename)
    
    opts = detectImportOptions(filename);
    opts.DataRange = [1 Inf]; %rows;
    cols = length(opts.VariableNames); %finding number of columns;
    opts.SelectedVariableNames = [1, 1:cols];
    A = readtable('kExperimentalData.xlsx',opts);
    
end

