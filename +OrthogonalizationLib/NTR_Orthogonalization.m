classdef NTR_Orthogonalization < handle

    properties
        numOfWords;
        numOfIterations;
        increment;
        output_path;
        parameters;
        jpglove;
        ORTable;
        scope;
        ELP;
        biphone;
    end

    methods
        % Constructor for the NTR_Orthogonalization class
        function obj = NTR_Orthogonalization(wordlistinput_file)
            % Load the word list
            wordlistinput = readtable(wordlistinput_file, 'ReadVariableNames', false);
            obj.numOfWords = size(wordlistinput, 2);
            obj.numOfIterations = size(wordlistinput, 1);
            obj.increment = obj.numOfIterations;
            obj.output_path = pwd;

            % Load the data files
            obj.parameters = readtable('ntr_masterlist_gp.xlsx');
            obj.jpglove = readtable('jpglove.csv');
            obj.ORTable = readtable('ntr_masterlist_onset_rimes.xlsx');
            obj.ORTable = table2array(obj.ORTable);
            obj.scope = readtable('ntr_masterlist_scope_upd.csv');
            obj.ELP = readtable('ntr_masterlist_elp_with_values_upd.xlsx', 'Sheet', 'Values');
            obj.biphone = readtable('iphod_wohoms_phonprob_edit.csv');
        end

        % Method to perform all the filter operations
        function filterData(obj)
            obj.filterHomographs();
            obj.filterWordLength();
            obj.filterGLOVEParameters();
            obj.filterELPbyMorphemesAndPOS();
            obj.filterScopeByMissingValues();
            obj.filterBiphoneProbabilityByMissingValues();
            obj.matchInputsToParameters();
            obj.formatTables();
        end

        % Method to filter out homographs from the parameters table
        function filterHomographs(obj)
            [list_size, ~] = size(obj.parameters);
            [~, ia, ic] = unique(obj.parameters.string); 
            v = accumarray(ic, 1); 
            obj.parameters = obj.parameters(ia(v == 1), :);
            obj.parameters = sortrows(obj.parameters, 'Item');
            disp([num2str(list_size - size(obj.parameters, 1)), ' words filtered (homographs)']);
            disp([num2str(size(obj.parameters, 1)), ' words in remaining list']);
        end

        % Method to filter out words by their length from the parameters table
        function filterWordLength(obj)
            [list_size, ~] = size(obj.parameters);
            mask2 = (strlength(obj.parameters.string) >= 3) & (strlength(obj.parameters.string) <= 9);
            obj.parameters = obj.parameters(mask2, :);
            disp([num2str(list_size - size(obj.parameters, 1)), ' words filtered (word length)']);
            disp([num2str(size(obj.parameters, 1)), ' words in remaining list']);
        end

        % Method to filter the GLOVE parameters table
        function filterGLOVEParameters(obj)
            [list_size, ~] = size(obj.parameters);
            obj.jpglove(~ismember(obj.jpglove.Var1, obj.parameters.string), :) = [];
            disp([num2str(list_size - size(obj.jpglove, 1)), ' words filtered (glove)']);
            obj.parameters(~ismember(obj.parameters.string, obj.jpglove.Var1), :) = [];
            disp([num2str(size(obj.parameters, 1)), ' words in remaining list']);
        end

        % Method to filter the ELP table by Morphemes and POS
        function filterELPbyMorphemesAndPOS(obj)
            [list_size, ~] = size(obj.parameters);
            mask = obj.ELP.MorphemeCount > 0 & ~isnan(obj.ELP.POS);
            obj.ELP = obj.ELP(mask, :);
            disp([num2str(list_size - size(obj.ELP, 1)), ' words filtered (ELP Morphemes and POS)']);
            disp([num2str(size(obj.ELP, 1)), ' words in remaining list']);
        end

        % Method to filter the Scope table by missing values
        function filterScopeByMissingValues(obj)
            [list_size, ~] = size(obj.parameters);
            obj.scope = rmmissing(obj.scope);
            disp([num2str(list_size - size(obj.scope, 1)), ' words filtered (missing scope values)']);
            disp([num2str(size(obj.scope, 1)), ' words in remaining list']);
        end

        % Method to filter the Biphone table by missing probabilities
        function filterBiphoneProbabilityByMissingValues(obj)
            [list_size, ~] = size(obj.parameters);
            obj.biphone = rmmissing(obj.biphone);
            disp([num2str(list_size - size(obj.biphone, 1)), ' words filtered (missing biphone probabilities)']);
            disp([num2str(size(obj.biphone, 1)), ' words in remaining list']);
        end

        % Placeholder function to match inputs to parameters
        function matchInputsToParameters(obj)
            % Implementation depends on what inputs need to be matched to what parameters
        end

        % Placeholder function to format tables
        function formatTables(obj)
            % Implementation depends on the desired format of the tables
        end

    end

end
