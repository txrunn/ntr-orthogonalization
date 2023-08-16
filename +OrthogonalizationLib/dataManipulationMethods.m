function matchInputsToParameters(obj)
    % Implementation depends on what inputs need to be matched to what parameters
end

% Placeholder

function formatTables(obj)

    % Read masterlist into table 'parameters' from 'ntr_masterlist_gp.xlsx'
    parameters = readtable('data/ntr_masterlist_gp.xlsx');

    % Read data into table 'jpglove' from 'jpglove.csv'
    jpglove = readtable('data/jpglove.csv');

    % Read onset-rime input into table 'ORTable' from 'ntr_masterlist_onset_rimes.xlsx'
    ORTable = readtable('data/ntr_masterlist_onset_rimes.xlsx');

end
