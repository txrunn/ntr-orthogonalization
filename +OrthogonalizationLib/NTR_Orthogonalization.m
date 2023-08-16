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

            % Read masterlist into table 'parameters' from 'ntr_masterlist_gp.xlsx'
            parameters = readtable('data/ntr_masterlist_gp.xlsx');

            % Read data into table 'jpglove' from 'jpglove.csv'
            jpglove = readtable('data/jpglove.csv');

            % Read onset-rime input into table 'ORTable' from 'ntr_masterlist_onset_rimes.xlsx'
            ORTable = readtable('data/ntr_masterlist_onset_rimes.xlsx');

        end

    end

end

function obj = NTR_Orthogonalization(obj)
    % % Constructor for the NTR_Orthogonalization class

    % Read onset-rime input into table 'ORTable' from 'onset_rimes.xlsx'
    ORTable = readtable('data/onset_rimes.xlsx');
end

function obj = size(obj)
    % % Method to perform all the filter operations

    % Determine the size of relevant data structures or arrays here.

    function filterHomographs(obj)
        % Filters homographs from the parameters
        homographs_filtered = [num2str(obj.list_size - size(obj.parameters, 1)), ' words filtered (homographs)'];
        disp(homographs_filtered);
    end

    function filterWordLength(obj)
        % Filters words based on length criteria
        word_length_filtered = [num2str(obj.list_size - size(obj.parameters, 1)), ' words filtered (word length)'];
        disp(word_length_filtered);
    end

    function filterGLOVEParameters(obj)
        % Placeholder for filtering words based on GLOVE parameters

        % Read GLOVE data
        glove = readtable('data/jpglove.csv');

        % Filter out rows in the glove table that do not match with the parameters string
        glove(~ismember(glove.Var1, parameters.string), :) = [];

        % Remove rows from the glove table where the Var2 column contains NaN values
        gloveindex = find(isnan(glove.Var2));
        glove(gloveindex, :) = [];

        % Generate a message indicating how many words were filtered out due to GLOVE filtering
        glove_filtered = [num2str(list_size - size(glove, 1)), ' words filtered (glove)'];

        disp('Filtering based on GLOVE parameters not yet implemented.');
    end

    function filterELPData(obj)
        % Placeholder for filtering and operations based on ELP data

        % Read ELP data
        ELP = readtable('data/ntr_masterlist_elp_with_values_upd.xlsx', 'Sheet', 'Values');

        % Filter out rows in the ELP table where the Occurrences column contains NaN values
        ELP((isnan(ELP.Occurrences)), :) = [];

        % Remove rows from the ELP table that do not match with the parameters string
        ELP(~ismember(ELP.Word, parameters.string), :) = [];

        % Further filtering based on NMorph, POS, and MorphSp columns
        % This is a simplified representation of the filters; more specific filtering criteria can be added based on original script
        ELP(ELP.NMorph > 3, :) = [];
        ELP(ismember(ELP.POS, {'minor', 'encl', '#'}), :) = [];
        ELP(ismember(ELP.MorphSp, {'>s', '>ed', '>ing'}), :) = [];

        % Generate a message indicating how many words were filtered out due to ELP filtering
        ELP_filtered = [num2str(list_size - size(ELP, 1)), ' words filtered (ELP)'];

        disp('Operations based on ELP data not yet implemented.');
    end

    function filterScopeData(obj)
        % Placeholder for filtering and operations based on scope data

        % Read scope data
        scope = readtable('data/ntr_masterlist_scope_upd.csv');

        % Filter out rows in the scope table that do not match with the parameters string
        scope(~ismember(scope.Word, parameters.string), :) = [];

        % Remove rows from the scope table based on NaN value checks in multiple columns
        nan_columns = {'Freq_SUBTLEXUS', 'BigramF_Avg_C_Log', 'TrigramF_Avg_C_Log', 'OLD20',
                       'PLD20', 'Sem_N_D', 'Phonographic_N'};

        for col = nan_columns
            col_name = col{1};
            scope(isnan(scope.(col_name)), :) = [];
        end

        % Generate a message indicating how many words were filtered out due to scope filtering
        scope_filtered = [num2str(list_size - size(scope, 1)), ' words filtered (scope)'];

        disp('Operations based on scope data not yet implemented.');
    end

    function filterBiphoneProbabilityData(obj)
        % Placeholder for filtering and operations based on biphone probability data

        % Read biphone probability data
        biphone = readtable('data/iphod_wohoms_phonprob_edit.csv');

        % Filter out rows in the biphone table that do not match with the parameters string
        biphone(~ismember(biphone.Word, parameters.string), :) = [];

        % Remove rows from the biphone table where the BiphonProb column contains NaN values
        biphone(isnan(biphone.BiphonProb), :) = [];

        % Generate a message indicating how many words were filtered out due to biphone probability filtering
        biphone_filtered = [num2str(list_size - size(biphone, 1)), ' words filtered (biphone prob)'];

        disp('Operations based on biphone probability data not yet implemented.');
    end

    function manipulateStrings(obj, inputString)
        % Placeholder for string manipulations

        % Various string-based filtering conditions and manipulations
        parameters = stringManipulations.performStringFilters(parameters, list_size);
        parameters = stringManipulations.performStringConversions(parameters);

        % Tokenized document processing
        parameters = stringManipulations.tokenizeStrings(parameters);

        % File naming and saving operations
        filenames = stringManipulations.generateFilenames(parameters, output_code);

        trimmedString = strtrim(inputString); % Removes leading and trailing white space from strings
        stringWithoutPlus = strrep(trimmedString, '+', ''); % Replaces '+' with empty string
        disp(['Original String: ', inputString]);
        disp(['Modified String: ', stringWithoutPlus]);
    end

    function matchInputsToParameters(obj)
        % Placeholder for matching other inputs to parameters

        % Read the input data
        input = readtable("data/wordinput_1.csv", 'ReadVariableNames', false);

        % Filter out rows in the input table that do not match with the parameters string
        input(~ismember(input.Word, parameters.string), :) = [];

        % Further transformations or processing on the input list as required

        disp('Matching inputs to parameters not yet implemented.');
    end

    function formatTables(obj)
        % Placeholder for formatting tables related to gp, bigp, and onset-rimes

        % Reading data from various files into tables
        inputTable = readtable("data/wordinput_1.csv", 'ReadVariableNames', false);

        % Filtering rows in the table based on certain criteria
        inputTable(~ismember(inputTable.Word, parameters.string), :) = [];

        % Converting tables to arrays and vice versa
        tableArray = table2array(inputTable);
        newArrayTable = array2table(tableArray);

        % Accessing specific rows or columns in a table
        specificRow = inputTable(5, :);
        specificColumn = inputTable.Word;

        % Further table manipulations and processing as required

        disp('Formatting tables not yet implemented.');
    end

    function runIterations(obj)
        % Placeholder for running iterations and computations

        % Sample logic for running iterations and computations:

        % Iterating over a specific size
        for g = 1:list_size
            % Relevant logic or computations for each iteration
        end

        % Running a specific number of iterations
        for a = 1:numOfIterations
            % Relevant logic or computations for each iteration
        end

        % Iterating over the number of words
        for s = 1:numOfWords
            % Relevant logic or computations for each iteration
        end

        % Iterating through specific tables and their rows and columns
        for t = 1:size(wp_tally, 1)

            for j = 1:numOfWords
                % Relevant logic or computations for each iteration
            end

        end

        % Further iterations and computations as required

        disp('Running iterations and computations not yet implemented.');
    end

    function processNextSegment(obj)
        % Placeholder for operations related to the next segment of the script

        % Matching word indices from the input list based on the number of words
        for s = 1:numOfWords
            r(1, s) = find(strcmp(wordlistinput{a, s}, string(parameters.string)) == 1);
        end

        % Initializing wordkey from the parameters
        wordkey = parameters.string;

        % Creating a new wordlist based on the found indices
        wordlist = strings(numOfWords, 1);

        for i = 1:numOfWords
            wordlist(i) = wordkey(r(i));
        end

        % Converting the word input list into a vertical table format
        % Additional logic as needed

        disp('Operations related to the next segment of the script not yet implemented.');
    end

    function processFurtherSegment(obj)
        % Placeholder for operations related to the further segment of the script

        % Takes wordinput list and makes it to a vertical table for parameter output
        % Generates table with words and its parameters
        r2 = rows2vars(array2table(r));
        r2(:, 1) = [];
        r2 = renamevars(r2, 'Var1', 'Index');
        wordlistOutput = [r2, array2table(wordlist)];

        % Additional logic as needed based on the original script

        disp('Operations related to the further segment of the script not yet implemented.');
    end

    function processSubsequentSegment(obj)
        % Placeholder for operations related to the subsequent segment of the script

        % Indexing parameters table and adding it to final output table
        for z = 1:numOfWords
            wordlistOutput(z, 1:width(parameters)) = parameters(r(z), :);
            % Ignore warning: The new variables being added to the table have fewer rows than the table.
        end

        % Save words and parameters as a .csv file
        save(strcat('allWordParameters_', num2str(a), '.csv'), 'wordlistOutput');

        % Indexing the words
        % Additional logic for word indexing based on the original script

        disp('Operations related to the subsequent segment of the script not yet implemented.');
    end

    function processNextSegment2(obj)
        % Placeholder for operations related to the next segment of the script

        % Matching word indices from the input list based on the number of words
        for s = 1:numOfWords
            r(1, s) = find(strcmp(wordlistinput{a, s}, string(parameters.string)) == 1);
        end

        % Initializing wordkey from the parameters
        wordkey = parameters.string;

        % Creating a new wordlist based on the found indices
        wordlist = strings(numOfWords, 1);

        for i = 1:numOfWords
            wordlist(i) = wordkey(r(i));
        end

        % Converting the word input list into a vertical table format
        % Additional logic as needed

        disp('Operations related to the next segment of the script (2) not yet implemented.');
    end

    function processUpcomingSegment(obj)
        % Placeholder for operations related to the upcoming segment of the script
        % Logic for the upcoming segment of the script has been addressed or is empty
        disp('Operations related to the upcoming segment of the script not yet implemented.');
    end

    function prepareDataForOrthogonalization(obj)
        % Placeholder for operations related to the subsequent segment of the script

        % Indexing parameters table and adding it to final output table
        for z = 1:numOfWords
            wordlistOutput(z, 1:width(parameters)) = parameters(r(z), :);
            % Ignore warning: The new variables being added to the table have fewer rows than the table.
        end

        % Save words and parameters as a .csv file
        save(strcat('allWordParameters_', num2str(a), '.csv'), 'wordlistOutput');

        % Indexing the words
        % Additional logic for word indexing based on the original script

        disp('Operations related to the subsequent segment of the script (2) not yet implemented.');
    end

    function processFollowingSegment(obj)
        % Placeholder for operations related to the following segment of the script
        % Logic for the following segment of the script has been addressed or is empty
        disp('Operations related to the following segment of the script not yet implemented.');
    end

    function processNextSegment3(obj)
        % Placeholder for operations related to the next segment of the script

        % Matching word indices from the input list based on the number of words
        for s = 1:numOfWords
            r(1, s) = find(strcmp(wordlistinput{a, s}, string(parameters.string)) == 1);
        end

        % Initializing wordkey from the parameters
        wordkey = parameters.string;

        % Creating a new wordlist based on the found indices
        wordlist = strings(numOfWords, 1);

        for i = 1:numOfWords
            wordlist(i) = wordkey(r(i));
        end

        % Converting the word input list into a vertical table format
        % Additional logic as needed

        disp('Operations related to the next segment of the script (3) not yet implemented.');
    end

    function performOrthogonalization(obj)
        % Placeholder for operations related to the subsequent segment of the script

        % Indexing parameters table and adding it to final output table
        for z = 1:numOfWords
            wordlistOutput(z, 1:width(parameters)) = parameters(r(z), :);
            % Ignore warning: The new variables being added to the table have fewer rows than the table.
        end

        % Save words and parameters as a .csv file
        save(strcat('allWordParameters_', num2str(a), '.csv'), 'wordlistOutput');

        % Indexing the words
        % Additional logic for word indexing based on the original script

        disp('Operations related to the subsequent segment of the script (3) not yet implemented.');
    end

    function processNextSegment4(obj)
        % Placeholder for operations related to the next segment of the script

        % Matching word indices from the input list based on the number of words
        for s = 1:numOfWords
            r(1, s) = find(strcmp(wordlistinput{a, s}, string(parameters.string)) == 1);
        end

        % Initializing wordkey from the parameters
        wordkey = parameters.string;

        % Creating a new wordlist based on the found indices
        wordlist = strings(numOfWords, 1);

        for i = 1:numOfWords
            wordlist(i) = wordkey(r(i));
        end

        % Converting the word input list into a vertical table format
        % Additional logic as needed

        disp('Operations related to the next segment of the script (4) not yet implemented.');
    end

    function postProcessResults(obj)
        % Placeholder for operations related to the subsequent segment of the script

        % Indexing parameters table and adding it to final output table
        for z = 1:numOfWords
            wordlistOutput(z, 1:width(parameters)) = parameters(r(z), :);
            % Ignore warning: The new variables being added to the table have fewer rows than the table.
        end

        % Save words and parameters as a .csv file
        save(strcat('allWordParameters_', num2str(a), '.csv'), 'wordlistOutput');

        % Indexing the words
        % Additional logic for word indexing based on the original script

        disp('Operations related to the subsequent segment of the script (4) not yet implemented.');
    end

    function processFollowingSegment2(obj)
        % Placeholder for operations related to the following segment of the script
        % Logic for the following segment of the script has been addressed or is empty
        disp('Operations related to the following segment of the script (2) not yet implemented.');
    end

    function generateOutputFiles(obj)
        % Placeholder for operations related to the subsequent segment of the script

        % Indexing parameters table and adding it to final output table
        for z = 1:numOfWords
            wordlistOutput(z, 1:width(parameters)) = parameters(r(z), :);
            % Ignore warning: The new variables being added to the table have fewer rows than the table.
        end

        % Save words and parameters as a .csv file
        save(strcat('allWordParameters_', num2str(a), '.csv'), 'wordlistOutput');

        % Indexing the words
        % Additional logic for word indexing based on the original script

        disp('Operations related to the subsequent segment of the script (5) not yet implemented.');
    end

end
