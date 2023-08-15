%function NTR_Orthogonalization_v13_wordinput(output_code) % censor if running as SCRIPT
% example for FUNCTION: NTR_Orthogonalization_v11('test')
%% Creates orthogonalized matrix comparing phoneme units, grapheme units,
% phoneme + grapheme (GP) units, letter edit distance (Levenstein), and semantic distance
%% Set Parameters
output_code = 'finalist_lists'; % uncommment for SCRIPT
%numOfIterations = 10000;
wordlistinput = readtable("wordinput_1.csv", 'ReadVariableNames', false);
numOfWords = size(wordlistinput, 2);
numOfIterations = size(wordlistinput, 1);
increment = numOfIterations;
%increment= 30;

output_path = pwd; % navigate to the program folder and uncenser this for running as SCRIPT
%output_path=fileparts(mfilename('fullpath'));% censor if running as SCRIPT. Get output path from the m-file location.

%% Read in info for SCRIPT (uncensor if running as script)
parameters = readtable('ntr_masterlist_gp.xlsx'); % creates table from masterlist input (10701 words)
jpglove = readtable('jpglove.csv'); %creates table from jpglove input (10701 words)
ORTable = readtable("ntr_masterlist_onset_rimes.xlsx"); % creates table from onset-rime input (10701 words)
ORTable = table2array(ORTable); %formats table from onset-rime input
scope = readtable("ntr_masterlist_scope_upd.csv");
ELP = readtable('ntr_masterlist_elp_with_values_upd.xlsx', 'Sheet', 'Values'); % removed words repated with ed/s/ing morphemes
biphone = readtable('iphod_wohoms_phonprob_edit.csv');

%% Read in info for function (censor if running as SCRIPT)
% parameters = readtable(fullfile(output_path,'ntr_masterlist_gp.xlsx')); % creates table from masterlist input (10701 words)
% jpglove=readtable('jpglove.csv'); %creates table from jpglove input (10701 words)
% ORTable = readtable(fullfile(output_path,'ntr_masterlist_onset_rimes.xlsx')); % creates table from onset-rime input (10701 words)
% ORTable = table2array(ORTable); %formats table from onset-rime input
% scope = readtable(fullfile(output_path,'ntr_masterlist_scope_upd.csv'));
% ELP=readtable(fullfile(output_path,'ntr_masterlist_elp_with_values_upd.xlsx'),'Sheet','Values');% removed words repated with ed/s/ing morphemes
% biphone=readtable(fullfile(output_path,'iphod_wohoms_phonprob_edit.csv'));
% GPfreq = readtable(fullfile(output_path, 'Lookup_Toolkit'));

%% filter homographs in parameters
[list_size, ~] = size(parameters);
[~, ia, ic] = unique(parameters.string); % Unique Elements in string column
v = accumarray(ic, 1); % Tally Occurrences Of Rows in string column
parameters = parameters(ia(v == 1), :); % Keep Rows for string that Only Appear Once
parameters = sortrows(parameters, 'Item'); %Re-sort parameters by Item #
homographs_filtered = [num2str(list_size - size(parameters, 1)), ' words filtered (homographs)'];
disp(homographs_filtered)
list_size = [num2str(size(parameters, 1)), ' words in remaining list'];
disp(list_size)

% %% filters homophones in parameters
% [list_size,~] = size(parameters);
% [~,ia,ic] = unique(parameters.IPA); % Unique Elements in IPA column
% v = accumarray(ic, 1); % Tally Occurrences Of Rows in IPA column
% parameters = parameters(ia(v==1),:); % Keep Rows for IPAs that Only Appear Once
% parameters=sortrows(parameters,'Item'); %Re-sort parameters by Item #
% homophones_filtered = [num2str(list_size-size(parameters,1)),' words filtered (homophones)'];
% disp(homophones_filtered)
% list_size = [num2str(size(parameters,1)),' words in remaining list'];
% disp(list_size)

%% filter by wordlength
[list_size, ~] = size(parameters);
mask2 = false;

for g = 1:list_size
    wordlength = strlength(parameters.string(g));

    if (9 >= wordlength) && (wordlength >= 3)
        mask2(g, 1) = 1;
    else
        mask2(g, 1) = 0;
    end

end

parameters = parameters(mask2, :);
wordlength_filtered = [num2str(list_size - size(parameters, 1)), ' words filtered (word length)'];
disp(wordlength_filtered)
list_size = [num2str(size(parameters, 1)), ' words in remaining list'];
disp(list_size)
%% Filter words missing GLOVE parameters
[list_size, ~] = size(parameters);
jpglove(~ismember(jpglove.Var1, parameters.string), :) = [];
[numrowsgloveorig, ~] = size(jpglove);
jpgloveindex = find(isnan(jpglove.Var2));
jpglove(jpgloveindex, :) = [];
glove_filtered = [num2str(list_size - size(jpglove, 1)), ' words filtered (glove)'];
disp(glove_filtered)

parameters(~ismember(parameters.string, jpglove.Var1), :) = [];
[numrowsnew, ~] = size(parameters);
list_size = [num2str(size(parameters, 1)), ' words in remaining list'];
disp(list_size)

%% Filter ELP by # Morphemes and POS
[list_size, ~] = size(parameters);
ELP((isnan(ELP.Occurrences)), :) = [];
ELP(~ismember(ELP.Word, parameters.string), :) = [];
[numrowsELPorig, ~] = size(ELP);
%nmorphindex = find(isnan(ELP.NMorph));
%ELP(nmorphindex,:) = [];
%nmorphindex2 = find(ELP.NMorph>3);
%ELP(nmorphindex2,:) = [];
% index = strfind(ELP.POS,'minor');
% index = ~cellfun(@isempty, index);
% ELP(index,:) = [];
index = strfind(ELP.POS, 'encl');
index = ~cellfun(@isempty, index);
ELP(index, :) = [];
% index = strfind(ELP.POS,'#');
% index = ~cellfun(@isempty, index);
% ELP(index,:) = [];
% index = strfind(ELP.MorphSp,'>s');
% index = ~cellfun(@isempty, index);
% ELP(index,:) = [];
% index = strfind(ELP.MorphSp,'>ed');
% index = ~cellfun(@isempty, index);
% ELP(index,:) = [];
% index = strfind(ELP.MorphSp,'>ing');
% index = ~cellfun(@isempty, index);
%ELP(index,:) = [];
ELP_filtered = [num2str(list_size - size(ELP, 1)), ' words filtered (ELP)'];
disp(ELP_filtered)

parameters(~ismember(parameters.string, ELP.Word), :) = [];
[numrowsnew, ~] = size(parameters);
list_size = [num2str(size(parameters, 1)), ' words in remaining list'];
disp(list_size)
%% Filter Scope by missing values
[list_size, ~] = size(parameters);
scope(~ismember(scope.Word, parameters.string), :) = [];
[numrowsscopeorig, ~] = size(scope);
scopeindex = find(isnan(scope.Freq_SUBTLEXUS));
scope(scopeindex, :) = [];
scopeindex2 = find(isnan(scope.BigramF_Avg_C_Log));
scope(scopeindex2, :) = [];
scopeindex3 = find(isnan(scope.TrigramF_Avg_C_Log));
scope(scopeindex3, :) = [];
scopeindex4 = find(isnan(scope.OLD20));
scope(scopeindex4, :) = [];
scopeindex5 = find(isnan(scope.PLD20));
scope(scopeindex5, :) = [];
scopeindex6 = find(isnan(scope.Sem_N_D));
scope(scopeindex6, :) = [];
scopeindex7 = find(isnan(scope.Phonographic_N));
scope(scopeindex7, :) = [];
scope_filtered = [num2str(list_size - size(scope, 1)), ' words filtered (scope)'];
disp(scope_filtered)

parameters(~ismember(parameters.string, scope.Word), :) = [];
[numrowsnew, ~] = size(parameters);
list_size = [num2str(size(parameters, 1)), ' words in remaining list'];
disp(list_size)

%% Filter biphone probability by missing values
[list_size, ~] = size(parameters);
biphone(~ismember(biphone.Word, parameters.string), :) = [];
[numrowsbiphoneorig, ~] = size(biphone);
biphoneindex = find(isnan(biphone.BiphonProb));
biphone(biphoneindex, :) = [];
biphone_filtered = [num2str(list_size - size(biphone, 1)), ' words filtered (biphone prob)'];
disp(biphone_filtered)

parameters(~ismember(parameters.string, biphone.Word), :) = [];
[numrowsnew, ~] = size(parameters);
list_size = [num2str(size(parameters, 1)), ' words in remaining list'];
disp(list_size)

disp('Filtering complete');

%% Match other inputs to parameters
% parameters(~ismember(parameters.string,jpglove.Var1),:)=[];
% parameters(~ismember(parameters.string,ELP.Word),:)=[];
% parameters(~ismember(parameters.string,scope.Word),:)=[];
% parameters(~ismember(parameters.string,biphone.Word),:)=[];
jpglove(~ismember(jpglove.Var1, parameters.string), :) = [];
scope(~ismember(scope.Word, parameters.string), :) = [];
ELP(~ismember(ELP.Word, parameters.string), :) = [];
biphone(~ismember(biphone.Word, parameters.string), :) = [];
%bigpTable(~ismember(bigpTable.bigpTable1,parameters.string),:)=[];
ORTable = array2table(ORTable);
ORTable(~ismember(ORTable.ORTable1, parameters.string), :) = [];

%% Format gp, bigp, onset-rimes tables
gpTablewords = parameters{:, {'string'}}; % formats empty table for stimulus words, finds every pg* column, makes table of just gp
mask = startsWith(parameters.Properties.VariableNames, 'pg'); % determines what columns in ntr_word_masterlist begin with lowercase "pg"
Fcols = parameters(:, mask); % omits columns that do not begin with "pg"
gpTable = table2cell([gpTablewords Fcols]); % combines gpTablewords and "pg" columns and converts table to cell
%bigpTable = table2cell(bigpTable);
ORTable = table2cell(ORTable);
[numrowsgp, numcolsgp] = size(gpTable); %Get # of pg columns to use later in script
%[numrowsbigp, numcolsbigp] = size(bigpTable); %Get bipg columns#, needed?
[numrowsOR, numcolsOR] = size(ORTable);
slash = "/"; %retain slash for bigp edit distance
%%
tic
it = 0;
tally = 0;

for a = 1:numOfIterations
    fprintf("iteration %d\n", a); % printing a start for visual cue; fprintf writes data to text file; %d start\n,a prints row for each n of iteration followed by "start"

    for s = 1:numOfWords
        r(1, s) = find(strcmp(wordlistinput{a, s}, string(parameters.string)) == 1); % JP sitched a and s variables
    end

    wordkey = parameters.string;
    wordlist = strings(numOfWords, 1);

    for i = 1:numOfWords %loop respective rows from word key for first N (number of words) numbers from r
        wordlist(i) = wordkey(r(i)); % %wordlist and randWords are the same
    end

    % Takes wordinput list and makes it to a vertical table for parameter output
    % Generates table with words and its parameters
    r2 = rows2vars(array2table(r));
    r2(:, 1) = [];
    r2 = renamevars(r2, 'Var1', 'Index');
    wordlistOutput = [r2, array2table(wordlist)];

    %Indexing parameters table and adding it to final output table
    for z = 1:numOfWords
        wordlistOutput(z, 1:width(parameters)) = parameters(r(z), :); % Ignore warning: The new variables being added to the table have fewer rows than the table.
    end

    save(strcat('allWordParameters_', num2str(a), '.csv'), 'wordlistOutput') % saves words and paramaters as .csv file

    wordindex(a, 2:[numOfWords + 1]) = transpose(wordlist(1:numOfWords)); %ERROR
    wordindex(a, 1) = sprintf('%d', a); %iteration
    %     % creates triangles for pairwise distances
    %     w = pdist(randWords); % determines pairwise distance
    %     E=squareform(w); % formats distance matrix
    %     PTriangles=triu(E,-1); % grabs upper triangular part of matrix from one diagonal below the center

    % Get POS, NMorph, NSyll Counts for Word List
    countELP = ELP(:, {'Word', 'POS', 'NMorph', 'NSyll'});
    countELPlist = countELP(r, :);

    nouncount(a, 1) = sum(strcmp(countELPlist.POS, 'NN'));
    verbcount(a, 1) = sum(strcmp(countELPlist.POS, 'VB'));
    adjcount(a, 1) = sum(strcmp(countELPlist.POS, 'JJ'));
    advcount(a, 1) = sum(strcmp(countELPlist.POS, 'RB'));
    multcount(a, 1) = size(find(contains(countELPlist.POS, '|')), 1);
    othercount(a, 1) = size(wordlist, 1) - nouncount(a, 1) - verbcount(a, 1) - adjcount(a, 1) - advcount(a, 1) - multcount(a, 1);

    onemorphcount(a, 1) = length(find(countELPlist.NMorph == 1));
    twomorphcount(a, 1) = length(find(countELPlist.NMorph == 2));
    threemorphcount(a, 1) = length(find(countELPlist.NMorph == 3));
    fourmorphcount(a, 1) = length(find(countELPlist.NMorph == 4));
    othermorphcount(a, 1) = size(wordlist, 1) - onemorphcount(a, 1) - twomorphcount(a, 1) - threemorphcount(a, 1) - fourmorphcount(a, 1);

    onesyllcount(a, 1) = length(find(countELPlist.NSyll == 1));
    twosyllcount(a, 1) = length(find(countELPlist.NSyll == 2));
    threesyllcount(a, 1) = length(find(countELPlist.NSyll == 3));
    foursyllcount(a, 1) = length(find(countELPlist.NSyll == 4));
    fivesyllcount(a, 1) = length(find(countELPlist.NSyll == 5));
    sixsyllcount(a, 1) = length(find(countELPlist.NSyll == 6));
    othersyllcount(a, 1) = size(wordlist, 1) - onesyllcount(a, 1) - twosyllcount(a, 1) - threesyllcount(a, 1) - foursyllcount(a, 1) - fivesyllcount(a, 1) - sixsyllcount(a, 1);
    %poslabels={'Noun';'Verb';'Adjective';'Adverb';'minor';'encl';'#';'Multiple'};
    %     poslabels={'Noun';'Verb';'Adjective';'Adverb';'Multiple';'Other'};
    %     POScount=table(poslabels);
    %     POScount{1,2}=num2cell(sum(strcmp(countELPlist.POS,'NN')));
    %     POScount{2,2}=num2cell(sum(strcmp(countELPlist.POS,'VB')));
    %     POScount{3,2}=num2cell(sum(strcmp(countELPlist.POS,'JJ')));
    %     POScount{4,2}=num2cell(sum(strcmp(countELPlist.POS,'RB')));
    %     POScount{5,2}=num2cell(size(find(contains(countELPlist.POS,'|')),1));
    %     POScount{6,2}=num2cell(size(wordlist,1)-sum(cell2mat(POScount{1:5,2}));
    %     POScount{5,2}=num2cell(sum(strcmp(countELPlist.POS,'minor')));
    %     POScount{6,2}=num2cell(sum(strcmp(countELPlist.POS,'encl')));
    %     POScount{7,2}=num2cell(sum(strcmp(countELPlist.POS,'#')));
    %     NMorphcount=table(num2cell(unique(ELP.NMorph)));
    %     for p=1:4
    %         NMorphcount{p,2}=length(find(countELPlist.NMorph==p));
    %     end
    %     NSyllcount=table(num2cell(unique(ELP.NSyll)));
    %     for w=1:4
    %         NSyllcount{w,2}=length(find(countELPlist.NSyll==w));
    %     end

    tri_mask = triu(ones(numOfWords) > 0, 1);

    % Setting up the results table
    letterTri = zeros(numOfWords);
    semanticTri = zeros(numOfWords);
    gpTri = zeros(numOfWords);
    % bigpTri = array2table(zeros(numOfWords));
    ORTri = zeros(numOfWords);
    phonemeTri = zeros(numOfWords);
    graphemeTri = zeros(numOfWords);
    Freq_SUBTLEXUS = zeros(numOfWords);
    OLD20sum = zeros(numOfWords);
    PLD20sum = zeros(numOfWords);
    Phonographic_N = zeros(numOfWords);
    Sem_N_D = zeros(numOfWords);
    BigramF_Avg_C_Log = zeros(numOfWords);
    GPfreq = zeros(numOfWords);
    Biphon_Prob = zeros(numOfWords);
    Length = zeros(numOfWords);
    NSyll = zeros(numOfWords);
    GPprob = zeros(numOfWords);
    ORprob = zeros(numOfWords);
    Gfreq = zeros(numOfWords);
    Pfreq = zeros(numOfWords);

    % get uniqe list
    wp_tally = nchoosek(1:numOfWords, 2);

    %finding the word pairings and their distances
    for t = 1:size(wp_tally, 1)
        % for j = 1:numOfWords
        %  for k = 1:numOfWords
        word1 = wordlist(wp_tally(t, 1), 1); % grabs word 1 from wordlist
        word2 = wordlist(wp_tally(t, 2), 1); % grabs word 2 from wordlist

        if (word1 < word2) %finding which word comes first in the alphabet
            lessAlpha = word1;
            moreAlpha = word2;
        else
            lessAlpha = word2;
            moreAlpha = word1;
        end

        % ORTHOGRAPHIC EDIT DISTANCE
        wLength(1) = strlength(word1);
        wLength(2) = strlength(word2);
        letterTri(wp_tally(t, 1), wp_tally(t, 2)) = editDistance(word1, word2, 'SwapCost', 1) / max(wLength);

        % GP EDIT DISTANCE - make function
        index1 = find(strcmp(word1, string(parameters.string)) == 1); %finding where word1 is in the unique word list
        index2 = find(strcmp(word2, string(parameters.string)) == 1); %find word2 in unique word list

        if isempty(index1) == 0 && isempty(index2) == 0 %check if string was found
            string1gp = strjoin(gpTable(index1, 2:end)); % condenses all 29 pg columns from gpTable into one string for word1
            string1gp(string1gp == '+') = []; % removes all plus signs
            word1GP = tokenizedDocument(string1gp); % converts to tokenized document (splits each unit into separate cells)

            string2gp = strjoin(gpTable(index2, 2:end)); % condenses all 29 pg columns from gpTable into one string for word1
            string2gp(string2gp == '+') = []; % removes all plus signs
            word2GP = tokenizedDocument(string2gp); % converts to tokenized document (splits each unit into separate cells)

            gpDistance = editDistance(word1GP, word2GP, 'SwapCost', 1); % gives edit distance for words 1 and 2 gps

            % string1gp = strtrim(string1gp); % removes leading and trailing white space from strings
            % string2gp = strtrim(string2gp);

            gpLength(1) = size(tokenDetails(word1GP), 1);
            gpLength(2) = size(tokenDetails(word2GP), 1);

            % wordLength(1) = count(string1gp,' ') + 1;
            % wordLength(2) = count(string2gp,' ') + 1;
            gpEditDistance = gpDistance(1, 1) / max(gpLength);

            gpTri(wp_tally(t, 1), wp_tally(t, 2)) = gpEditDistance;
        end

        % Onset rimes EDIT DISTANCE
        index1 = find(strcmp(word1, string(parameters.string)) == 1); %finding where word1 is in the unique word list
        index2 = find(strcmp(word2, string(parameters.string)) == 1); %find word2 in unique word list

        if isempty(index1) == 0 && isempty(index2) == 0 %check if string was found
            %string1OR = strjoin(ORTable(index1,2:end)); % condenses all 29 pg columns from gpTable into one string for word1
            for n = 2:numcolsOR
                string1OR = string(ORTable(index1, n));
                string1OR = strrep(string1OR, '+', '');
                string1ORfull(1, n) = string1OR;
                string2OR = string(ORTable(index2, n));
                string2OR = strrep(string2OR, '+', '');
                string2ORfull(1, n) = string2OR;
            end

            string1ORtok = tokenizedDocument(string1ORfull, 'TokenizeMethod', 'none');
            string2ORtok = tokenizedDocument(string2ORfull, 'TokenizeMethod', 'none');
            ORDistance = editDistance(string1ORtok, string2ORtok, 'SwapCost', 1);

            % string1OR = strtrim(string1OR); % removes leading and trailing white space from strings
            % string2OR = strtrim(string2OR);
            ORLength(1) = size(tokenDetails(string1ORtok), 1);
            ORLength(2) = size(tokenDetails(string2ORtok), 1);
            OREditDistance = ORDistance(1, 1) / max(ORLength);

            ORTri(wp_tally(t, 1), wp_tally(t, 2)) = OREditDistance;
        end

        % SEMANTIC EDIT DISTANCE - change to use file, not variables from workspace
        word1Vars = jpglove(r(1, wp_tally(t, 1)), 2:301); % grabs all glove variables for word 1
        word2Vars = jpglove(r(1, wp_tally(t, 2)), 2:301); % grabs all glove variables for word 2
        semanticCorr = corrcoef(table2array(word1Vars), table2array(word2Vars)); %grabs correlation coefficient
        semanticTri(wp_tally(t, 1), wp_tally(t, 2)) = 1 - semanticCorr(1, 2); % records correlation

        % PHONEME UNIT AND GRAPHEME UNIT EDIT DISTANCE (separated, not P+G)
        if isempty(index1) == 0 && isempty(index2) == 0 %check if string was found
            phonemeWord1 = string; % sets up four variables with empty strings
            graphemeWord1 = string;
            phonemeWord2 = string;
            graphemeWord2 = string;

            for m = 2:numcolsgp %Read the number of columns in gpTable (don't include words column)
                pString1 = strjoin(gpTable(index1, 2:end)); % check if p correct

                pgWord1 = gpTable{index1, m};
                pWord1 = extractBefore(pgWord1, "+"); %grabs just phoneme (before plus)
                phonemeWord1 = append(phonemeWord1, pWord1, ' '); % puts phonemes into one string

                pgWord1 = gpTable{index1, m};
                gWord1 = extractAfter(pgWord1, "+"); %grabs just grapheme (after plus)
                graphemeWord1 = append(graphemeWord1, gWord1, ' '); % puts graphemes into one string

                pString2 = strjoin(gpTable(index2, 2:end)); % check if p correct
                pgWord2 = gpTable{index2, m};
                pWord2 = extractBefore(pgWord2, "+");
                phonemeWord2 = append(phonemeWord2, pWord2, ' ');

                pgWord2 = gpTable{index2, m};
                gWord2 = extractAfter(pgWord2, "+");
                graphemeWord2 = append(graphemeWord2, gWord2, ' ');
            end

            % creates tokenized documents
            phonemeWord1 = tokenizedDocument(phonemeWord1);
            graphemeWord1 = tokenizedDocument(graphemeWord1);
            phonemeWord2 = tokenizedDocument(phonemeWord2);
            graphemeWord2 = tokenizedDocument(graphemeWord2);

            % calculates edit distances for phonemes and graphemes
            phonemeLength(1) = size(tokenDetails(phonemeWord1), 1);
            phonemeLength(2) = size(tokenDetails(phonemeWord2), 1);

            graphemeLength(1) = size(tokenDetails(graphemeWord1), 1);
            graphemeLength(2) = size(tokenDetails(graphemeWord2), 1);

            pEditDistance = editDistance(phonemeWord1, phonemeWord2, 'SwapCost', 1);
            gEditDistance = editDistance(graphemeWord1, graphemeWord2, 'SwapCost', 1);
            phonemeEditDistance = pEditDistance(1, 1) / max(phonemeLength);
            graphemeEditDistance = gEditDistance(1, 1) / max(graphemeLength);

            phonemeTri(wp_tally(t, 1), wp_tally(t, 2)) = phonemeEditDistance; % records phonemes and graphemes edit distances
            graphemeTri(wp_tally(t, 1), wp_tally(t, 2)) = graphemeEditDistance;
        end

        %% SUM DISTANCES
        index1 = find(strcmp(word1, string(scope.Word)) == 1); %finding where word1 is in the unique word list
        index2 = find(strcmp(word2, string(scope.Word)) == 1);

        Freq_SUBTLEXUS(wp_tally(t, 1), wp_tally(t, 2)) = scope.Freq_SUBTLEXUS(index1) + scope.Freq_SUBTLEXUS(index2);
        OLD20sum(wp_tally(t, 1), wp_tally(t, 2)) = scope.OLD20(index1) + scope.OLD20(index2);
        %PLD20sum(wp_tally(t,1),wp_tally(t,2)) = scope.PLD20(index1) + scope.PLD20(index2);
        Phonographic_N(wp_tally(t, 1), wp_tally(t, 2)) = scope.Phonographic_N(index1) + scope.Phonographic_N(index2);
        Sem_N_D(wp_tally(t, 1), wp_tally(t, 2)) = scope.Sem_N_D(index1) + scope.Sem_N_D(index2);
        BigramF_Avg_C_Log(wp_tally(t, 1), wp_tally(t, 2)) = scope.BigramF_Avg_C_Log(index1) + scope.BigramF_Avg_C_Log(index2);
        Biphon_Prob(wp_tally(t, 1), wp_tally(t, 2)) = biphone.BiphonProb(index1) + biphone.BiphonProb(index2);

        %GPprob(wp_tally(t,1),wp_tally(t,2)) = parameters.GP_mean(index1) + parameters.GP_mean(index2);
        % ORprob(wp_tally(t,1),wp_tally(t,2)) = parameters.ORGP_mean(index1) + parameters.ORGP_mean(index2);
        % Gfreq(wp_tally(t,1),wp_tally(t,2)) = parameters.G_freq_mean(index1) + parameters.G_freq_mean(index2);
        %Pfreq(wp_tally(t,1),wp_tally(t,2)) = parameters.P_freq_mean(index1) + parameters.P_freq_mean(index2);

        GPprob(wp_tally(t, 1), wp_tally(t, 2)) = parameters.GPprob_Min(index1) + parameters.GPprob_Min(index2);
        ORprob(wp_tally(t, 1), wp_tally(t, 2)) = parameters.ORprob_Min(index1) + parameters.ORprob_Min(index2);
        Gfreq(wp_tally(t, 1), wp_tally(t, 2)) = parameters.Gfreq_Mean(index1) + parameters.Gfreq_Mean(index2);
        Pfreq(wp_tally(t, 1), wp_tally(t, 2)) = parameters.Pfreq_Mean(index1) + parameters.Pfreq_Mean(index2);
        GPfreq(wp_tally(t, 1), wp_tally(t, 2)) = parameters.GPfreq_Min(index1) + parameters.GPfreq_Min(index2);

        %% SUBTRACTION DISTANCES
        index1 = find(strcmp(word1, string(ELP.Word)) == 1); %finding where word1 is in the unique word list
        index2 = find(strcmp(word2, string(ELP.Word)) == 1);

        Length(wp_tally(t, 1), wp_tally(t, 2)) = ELP.Length(index1) - ELP.Length(index2);
        NSyll(wp_tally(t, 1), wp_tally(t, 2)) = ELP.NSyll(index1) - ELP.NSyll(index2);
    end

    %fprintf("%d after triangle\n", a);

    % Spearman correlation comparing GP edit distance to phonological and orthographic
    GPprobToPhono(a, 1) = corr(phonemeTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToletter(a, 1) = corr(letterTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToSem(a, 1) = corr(semanticTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToGraph(a, 1) = corr(graphemeTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToFreq_SUBTLEXUS(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToOLD20sum(a, 1) = corr(OLD20sum(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToGPfreq(a, 1) = corr(GPfreq(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToPhonographic_N(a, 1) = corr(Phonographic_N(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToSem_N_D(a, 1) = corr(Sem_N_D(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToBigramF_Avg_C_Log(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToLength(a, 1) = corr(Length(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToNSyll(a, 1) = corr(NSyll(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToBiphon_Prob(a, 1) = corr(Biphon_Prob(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToOR(a, 1) = corr(ORTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobTogp(a, 1) = corr(gpTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToORprob(a, 1) = corr(ORprob(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToGfreq(a, 1) = corr (Gfreq(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToPfreq(a, 1) = corr (Pfreq(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');

    ORprobToPhono(a, 1) = corr(phonemeTri(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToletter(a, 1) = corr(letterTri(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToSem(a, 1) = corr(semanticTri(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToGraph(a, 1) = corr(graphemeTri(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToFreq_SUBTLEXUS(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToOLD20sum(a, 1) = corr(OLD20sum(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToGPfreq(a, 1) = corr(GPfreq(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToPhonographic_N(a, 1) = corr(Phonographic_N(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToSem_N_D(a, 1) = corr(Sem_N_D(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToBigramF_Avg_C_Log(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToLength(a, 1) = corr(Length(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToNSyll(a, 1) = corr(NSyll(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToBiphon_Prob(a, 1) = corr(Biphon_Prob(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToOR(a, 1) = corr(ORTri(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobTogp(a, 1) = corr(gpTri(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToGfreq(a, 1) = corr(Gfreq(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');
    ORprobToPfreq(a, 1) = corr(Pfreq(tri_mask), ORprob(tri_mask), 'Type', 'Spearman');

    gpToPhono(a, 1) = corr(phonemeTri(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToletter(a, 1) = corr(letterTri(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToSem(a, 1) = corr(semanticTri(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToGraph(a, 1) = corr(graphemeTri(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToFreq_SUBTLEXUS(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToOLD20sum(a, 1) = corr(OLD20sum(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToGPfreq(a, 1) = corr(GPfreq(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToPhonographic_N(a, 1) = corr(Phonographic_N(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToSem_N_D(a, 1) = corr(Sem_N_D(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToBigramF_Avg_C_Log(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToLength(a, 1) = corr(Length(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToNSyll(a, 1) = corr(NSyll(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToBiphon_Prob(a, 1) = corr(Biphon_Prob(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToOR(a, 1) = corr(ORTri(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToGfreq(a, 1) = corr(Gfreq(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');
    gpToPfreq(a, 1) = corr(Pfreq(tri_mask), gpTri(tri_mask), 'Type', 'Spearman');

    ORToPhono(a, 1) = corr(phonemeTri(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToletter(a, 1) = corr(letterTri(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToSem(a, 1) = corr(semanticTri(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToGraph(a, 1) = corr(graphemeTri(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToFreq_SUBTLEXUS(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToOLD20sum(a, 1) = corr(OLD20sum(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToGPfreq(a, 1) = corr(GPfreq(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToPhonographic_N(a, 1) = corr(Phonographic_N(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToSem_N_D(a, 1) = corr(Sem_N_D(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToBigramF_Avg_C_Log(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToLength(a, 1) = corr(Length(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToNSyll(a, 1) = corr(NSyll(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToBiphon_Prob(a, 1) = corr(Biphon_Prob(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToGfreq(a, 1) = corr(Gfreq(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');
    ORToPfreq(a, 1) = corr(Pfreq(tri_mask), ORTri(tri_mask), 'Type', 'Spearman');

    Freq_SUBTLEXUSToPhono(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), phonemeTri(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToLetter(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), letterTri(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToSem(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToGraph(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), graphemeTri(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToOLD20sum(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), OLD20sum(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToGPfreq(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), GPfreq(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToPhonographic_N(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), Phonographic_N(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToSem_N_D(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), Sem_N_D(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToBigramF_Avg_C_Log(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), BigramF_Avg_C_Log(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToLength(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToNSyll(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToBiphon_Prob(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToGfreq(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    Freq_SUBTLEXUSToPfreq(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    OLD20sumToGPfreq(a, 1) = corr(OLD20sum(tri_mask), GPfreq(tri_mask), 'Type', 'Spearman');
    OLD20sumToPhonographic_N(a, 1) = corr(OLD20sum(tri_mask), Phonographic_N(tri_mask), 'Type', 'Spearman');
    OLD20sumToSem_N_D(a, 1) = corr(OLD20sum(tri_mask), Sem_N_D(tri_mask), 'Type', 'Spearman');
    OLD20sumToBigramF_Avg_C_Log(a, 1) = corr(OLD20sum(tri_mask), BigramF_Avg_C_Log(tri_mask), 'Type', 'Spearman');
    OLD20sumToLength(a, 1) = corr(OLD20sum(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    OLD20sumToNSyll(a, 1) = corr(OLD20sum(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    OLD20sumToBiphon_Prob(a, 1) = corr(OLD20sum(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    OLD20sumToPhono(a, 1) = corr(OLD20sum(tri_mask), phonemeTri(tri_mask), 'Type', 'Spearman');
    OLD20sumToLetter(a, 1) = corr(OLD20sum(tri_mask), letterTri(tri_mask), 'Type', 'Spearman');
    OLD20sumToGraph(a, 1) = corr(OLD20sum(tri_mask), graphemeTri(tri_mask), 'Type', 'Spearman');
    OLD20sumToSem(a, 1) = corr(OLD20sum(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    OLD20sumToGfreq(a, 1) = corr(OLD20sum(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    OLD20sumToPfreq(a, 1) = corr(OLD20sum(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    GPfreqToPhonographic_N(a, 1) = corr(GPfreq(tri_mask), Phonographic_N(tri_mask), 'Type', 'Spearman');
    GPfreqToSem_N_D(a, 1) = corr(GPfreq(tri_mask), Sem_N_D(tri_mask), 'Type', 'Spearman');
    GPfreqToBigramF_Avg_C_Log(a, 1) = corr(GPfreq(tri_mask), BigramF_Avg_C_Log(tri_mask), 'Type', 'Spearman');
    GPfreqToLength(a, 1) = corr(GPfreq(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    GPfreqToNSyll(a, 1) = corr(GPfreq(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    GPfreqToBiphon_Prob(a, 1) = corr(GPfreq(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    GPfreqToSem(a, 1) = corr(GPfreq(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    GPfreqToPhono(a, 1) = corr(GPfreq(tri_mask), phonemeTri(tri_mask), 'Type', 'Spearman');
    GPfreqToLetter(a, 1) = corr(GPfreq(tri_mask), letterTri(tri_mask), 'Type', 'Spearman');
    GPfreqToGraph(a, 1) = corr(GPfreq(tri_mask), graphemeTri(tri_mask), 'Type', 'Spearman');
    GPfreqToGfreq(a, 1) = corr(GPfreq(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    GPfreqToPfreq(a, 1) = corr(GPfreq(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    Phonographic_NToSem_N_D(a, 1) = corr(Phonographic_N(tri_mask), Sem_N_D(tri_mask), 'Type', 'Spearman');
    Phonographic_NToBigramF_Avg_C_Log(a, 1) = corr(Phonographic_N(tri_mask), BigramF_Avg_C_Log(tri_mask), 'Type', 'Spearman');
    Phonographic_NToLength(a, 1) = corr(Phonographic_N(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    Phonographic_NToNSyll(a, 1) = corr(Phonographic_N(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    Phonographic_NToBiphon_Prob(a, 1) = corr(Phonographic_N(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    Phonographic_ToSem(a, 1) = corr(Phonographic_N(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    Phonographic_NToPhono(a, 1) = corr(Phonographic_N(tri_mask), phonemeTri(tri_mask), 'Type', 'Spearman');
    Phonographic_NToLetter(a, 1) = corr(Phonographic_N(tri_mask), letterTri(tri_mask), 'Type', 'Spearman');
    Phonographic_NTo_Graph(a, 1) = corr(Phonographic_N(tri_mask), graphemeTri(tri_mask), 'Type', 'Spearman');
    Phonographic_NToGfreq(a, 1) = corr(Phonographic_N(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    Phonographic_NToPfreq(a, 1) = corr(Phonographic_N(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    Sem_N_DToBigramF_Avg_C_Log(a, 1) = corr(Sem_N_D(tri_mask), BigramF_Avg_C_Log(tri_mask), 'Type', 'Spearman');
    Sem_N_DToLength(a, 1) = corr(Sem_N_D(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    Sem_N_DToNSyll(a, 1) = corr(Sem_N_D(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    Sem_N_DToBiphon_Prob(a, 1) = corr(Sem_N_D(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    Sem_N_DToSem(a, 1) = corr(Sem_N_D(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    Sem_N_DToPhono(a, 1) = corr(Sem_N_D(tri_mask), phonemeTri(tri_mask), 'Type', 'Spearman');
    Sem_N_DToLetter(a, 1) = corr(Sem_N_D(tri_mask), letterTri(tri_mask), 'Type', 'Spearman');
    Sem_N_DToGraph(a, 1) = corr(Sem_N_D(tri_mask), graphemeTri(tri_mask), 'Type', 'Spearman');
    Sem_N_DToGfreq(a, 1) = corr(Sem_N_D(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    Sem_N_DToPfreq(a, 1) = corr(Sem_N_D(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    BigramF_Avg_C_LogToLength(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    BigramF_Avg_C_LogToNSyll(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    BigramF_Avg_C_LogToBiphon_Prob(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    BigramF_Avg_C_LogToSem(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    BigramF_Avg_C_LogToPhono(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), phonemeTri(tri_mask), 'Type', 'Spearman');
    BigramF_Avg_C_LogToLetter(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), letterTri(tri_mask), 'Type', 'Spearman');
    BigramF_Avg_C_LogToGraph(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), graphemeTri(tri_mask), 'Type', 'Spearman');
    BigramF_Avg_C_LogToGfreq(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    BigramF_Avg_C_LogToPfreq(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    graphToletter(a, 1) = corr(graphemeTri(tri_mask), letterTri(tri_mask), 'Type', 'Spearman');
    graphToPhono(a, 1) = corr(graphemeTri(tri_mask), phonemeTri(tri_mask), 'Type', 'Spearman');
    graphToLength(a, 1) = corr(graphemeTri(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    graphToNSyll(a, 1) = corr(graphemeTri(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    graphToBiphon_Prob(a, 1) = corr(graphemeTri(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    graphToSem(a, 1) = corr(graphemeTri(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    graphToGfreq(a, 1) = corr(graphemeTri(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    graphToPfreq(a, 1) = corr(graphemeTri(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    phonoToletter(a, 1) = corr(phonemeTri(tri_mask), letterTri(tri_mask), 'Type', 'Spearman');
    phonoToLength(a, 1) = corr(phonemeTri(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    phonoToNSyll(a, 1) = corr(phonemeTri(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    phonoToBiphon_Prob(a, 1) = corr(phonemeTri(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    phonoToSem(a, 1) = corr(phonemeTri(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    phonoToGfreq(a, 1) = corr(phonemeTri(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    phonoToPfreq(a, 1) = corr(phonemeTri(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    LetterToLength(a, 1) = corr(letterTri(tri_mask), Length(tri_mask), 'Type', 'Spearman');
    LetterToNSyll(a, 1) = corr(letterTri(tri_mask), NSyll(tri_mask), 'Type', 'Spearman');
    LetterToBiphon_Prob(a, 1) = corr(letterTri(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    LetterToSem(a, 1) = corr(letterTri(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    LetterToGfreq(a, 1) = corr(letterTri(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    LetterToPfreq(a, 1) = corr(letterTri(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    LengthToBiphon_Prob(a, 1) = corr(Length(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    LengthToSem(a, 1) = corr(Length(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    LengthToGfreq(a, 1) = corr(Length(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    LengthToPfreq(a, 1) = corr(Length(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    NSyllToBiphon_Prob(a, 1) = corr(NSyll(tri_mask), Biphon_Prob(tri_mask), 'Type', 'Spearman');
    NSyllToSem(a, 1) = corr(NSyll(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    NSyllToGfreq(a, 1) = corr(NSyll(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    NSyllToPfreq(a, 1) = corr(NSyll(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    Biphon_ProbToSem_Corr(a, 1) = corr(Biphon_Prob(tri_mask), semanticTri(tri_mask), 'Type', 'Spearman');
    Biphon_ProbToGfreq(a, 1) = corr(Biphon_Prob(tri_mask), Gfreq(tri_mask), 'Type', 'Spearman');
    Biphon_ProbToPfreq(a, 1) = corr(Biphon_Prob(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    NSyllToLen(a, 1) = corr(NSyll(tri_mask), Length(tri_mask), 'Type', 'Spearman');

    GfreqToPfreq(a, 1) = corr (Gfreq(tri_mask), Pfreq(tri_mask), 'Type', 'Spearman');

    %Example of getting 1 stdev for the whole correlation matrix
    w = 0; %(w=0 is n-1, w=1 is n for stdev)
    std_phonemeTri(a, 1) = std(phonemeTri(tri_mask), w, "all");
    std_gpTri(a, 1) = std(gpTri(tri_mask), w, "all");
    std_letterTri(a, 1) = std(letterTri(tri_mask), w, "all");
    std_semanticTri(a, 1) = std(semanticTri(tri_mask), w, "all");
    std_graphemeTri(a, 1) = std(graphemeTri(tri_mask), w, "all");
    std_ORTri(a, 1) = std(ORTri(tri_mask), w, "all");
    std_Freq_SUBTLEXUS(a, 1) = std(Freq_SUBTLEXUS(tri_mask), w, "all");
    std_OLD20sum(a, 1) = std(OLD20sum(tri_mask), w, "all");
    std_GPfreq(a, 1) = std(GPfreq(tri_mask), w, "all");
    std_Phonographic_N(a, 1) = std(Phonographic_N(tri_mask), w, "all");
    std_Sem_N_D(a, 1) = std(Sem_N_D(tri_mask), w, "all");
    std_BigramF_Avg_C_Log(a, 1) = std(BigramF_Avg_C_Log(tri_mask), w, "all");
    std_Length(a, 1) = std(Length(tri_mask), w, "all");
    std_NSyll(a, 1) = std(NSyll(tri_mask), w, "all");
    std_Biphon_Prob(a, 1) = std(Biphon_Prob(tri_mask), w, "all");
    std_Gfreq(a, 1) = std(Gfreq(tri_mask), w, "all");
    std_Pfreq(a, 1) = std(Pfreq(tri_mask), w, "all");

    % Saves distance matrices into singular structure array and .mat file
    letterTri = fliplr(letterTri);
    semanticTri = fliplr(semanticTri);
    gpTri = fliplr(gpTri);
    ORTri = fliplr(ORTri);
    phonemeTri = fliplr(phonemeTri);
    graphemeTri = fliplr(graphemeTri);
    Freq_SUBTLEXUS = fliplr(Freq_SUBTLEXUS);
    OLD20sum = fliplr(OLD20sum);
    PLD20sum = fliplr(PLD20sum);
    Phonographic_N = fliplr(Phonographic_N);
    Sem_N_D = fliplr(Sem_N_D);
    BigramF_Avg_C_Log = fliplr(BigramF_Avg_C_Log);
    GPfreq = fliplr(GPfreq);
    Biphon_Prob = fliplr(Biphon_Prob);
    Length = fliplr(Length);
    NSyll = fliplr(NSyll);
    GPprob = fliplr(GPprob);
    ORprob = fliplr(ORprob);
    Gfreq = fliplr(Gfreq);
    Pfreq = fliplr(Pfreq);

    allMatrices = struct('letterTri', letterTri, 'semanticTri', semanticTri, 'gpTri', gpTri, ...
        'ORTri', ORTri, 'phonemeTri', phonemeTri, 'graphemeTri', graphemeTri, 'Freq_SUBTLEXUS', Freq_SUBTLEXUS, ...
        'OLD20sum', OLD20sum, 'PLD20sum', PLD20sum, 'Phonographic_N', Phonographic_N, ...
        'Sem_N_D', Sem_N_D, 'BigramF_Avg_C_Log', BigramF_Avg_C_Log, 'GPfreq', GPfreq, 'Biphon_Prob', Biphon_Prob, ...
        'Length', Length, 'NSyll', NSyll, 'GPprob', GPprob, 'ORprob', ORprob, 'Gfreq', Gfreq, 'Pfreq', Pfreq);

    disMatName = strcat('allDistanceMatrices_', num2str(a), '.mat');
    save(disMatName, 'allMatrices')

    % print out results every 30 iterations
    it = it + 1;
    tally = tally + 1;

    if tally == increment
        %% writes everything in excel spreadsheets
        iteration(1:it, 1) = (1:it)';
        %wordindex(:,1) = iteration; %Number rows 1 through number of iterations for word index
        % set up table based on variables

        finalResults = table(iteration, nouncount, verbcount, adjcount, advcount, multcount, othercount, onemorphcount, twomorphcount, threemorphcount, fourmorphcount, othermorphcount, ...
            onesyllcount, twosyllcount, threesyllcount, foursyllcount, fivesyllcount, sixsyllcount, othersyllcount, NSyllToLen, gpToPhono, gpToletter, gpToSem, gpToGraph, gpToFreq_SUBTLEXUS, ...
            gpToOLD20sum, gpToGPfreq, gpToPhonographic_N, gpToLength, gpToNSyll, gpToOR, ORToPhono, ORToletter, ORToSem, ORToGraph, ORToFreq_SUBTLEXUS, ORToLength, ORToNSyll, Freq_SUBTLEXUSToPhono, ...
            Freq_SUBTLEXUSToLetter, Freq_SUBTLEXUSToSem, Freq_SUBTLEXUSToGraph, Freq_SUBTLEXUSToOLD20sum, Freq_SUBTLEXUSToGPfreq, Freq_SUBTLEXUSToPhonographic_N, ...
            Freq_SUBTLEXUSToSem_N_D, Freq_SUBTLEXUSToBigramF_Avg_C_Log, Freq_SUBTLEXUSToLength, Freq_SUBTLEXUSToNSyll, ...
            Freq_SUBTLEXUSToBiphon_Prob, Freq_SUBTLEXUSToGfreq, Freq_SUBTLEXUSToPfreq, OLD20sumToGPfreq, OLD20sumToPhonographic_N, OLD20sumToSem_N_D, OLD20sumToBigramF_Avg_C_Log, OLD20sumToLength, OLD20sumToNSyll, ...
            OLD20sumToBiphon_Prob, OLD20sumToPhono, OLD20sumToLetter, OLD20sumToGraph, OLD20sumToSem, OLD20sumToGfreq, OLD20sumToPfreq, GPfreqToPhonographic_N, GPfreqToSem_N_D, ...
            GPfreqToBigramF_Avg_C_Log, GPfreqToLength, GPfreqToNSyll, GPfreqToBiphon_Prob, GPfreqToSem, GPfreqToGfreq, GPfreqToPfreq, ...
            Phonographic_NToSem_N_D, Phonographic_NToBigramF_Avg_C_Log, Phonographic_NToLength, Phonographic_NToNSyll, Phonographic_NToBiphon_Prob, ...
            Phonographic_NToGfreq, Phonographic_NToPfreq, Sem_N_DToBigramF_Avg_C_Log, Sem_N_DToLength, Sem_N_DToNSyll, ...
            Sem_N_DToBiphon_Prob, Sem_N_DToSem, Sem_N_DToPhono, Sem_N_DToLetter, Sem_N_DToGraph, Sem_N_DToGfreq, Sem_N_DToPfreq, BigramF_Avg_C_LogToLength, ...
            BigramF_Avg_C_LogToNSyll, BigramF_Avg_C_LogToBiphon_Prob, BigramF_Avg_C_LogToSem, ...
            BigramF_Avg_C_LogToGfreq, BigramF_Avg_C_LogToPfreq, graphToletter, graphToPhono, graphToLength, graphToNSyll, ...
            graphToBiphon_Prob, graphToSem, graphToGfreq, graphToPfreq, phonoToletter, phonoToLength, phonoToNSyll, phonoToBiphon_Prob, phonoToSem, phonoToGfreq, phonoToPfreq, LetterToLength, ...
            LetterToNSyll, LetterToSem, LengthToBiphon_Prob, LengthToSem, LengthToGfreq, LengthToPfreq, NSyllToBiphon_Prob, NSyllToGfreq, NSyllToPfreq, Biphon_ProbToGfreq, Biphon_ProbToPfreq, ...
            GPprobToGraph, GPprobToFreq_SUBTLEXUS, GPprobToOLD20sum, GPprobToGPfreq, ...
            GPprobToPhonographic_N, GPprobToSem_N_D, GPprobToBigramF_Avg_C_Log, GPprobToLength, GPprobToNSyll, GPprobToBiphon_Prob, ...
            GPprobToOR, GPprobTogp, GPprobToORprob, GPprobToGfreq, GPprobToPfreq, ORprobToPhono, ORprobToletter, ORprobToFreq_SUBTLEXUS, ...
            ORprobToOLD20sum, ORprobToGPfreq, ORprobToPhonographic_N, ORprobToSem_N_D, ORprobToBigramF_Avg_C_Log, ...
            ORprobToLength, ORprobToNSyll, ORprobToBiphon_Prob, ORprobToOR, ORprobTogp, ORprobToGfreq, ORprobToPfreq, ...
            std_phonemeTri, std_gpTri, std_letterTri, std_semanticTri, std_graphemeTri, std_ORTri, std_Freq_SUBTLEXUS, std_OLD20sum, ...
            std_GPfreq, std_Phonographic_N, std_Sem_N_D, std_BigramF_Avg_C_Log, std_Length, std_NSyll, std_Biphon_Prob, std_Gfreq, std_Pfreq);

        %         %% output for SCRIPT uncensor if running as SCRIPT

        writematrix(wordindex, fullfile('output', strcat(output_code, sprintf('_Orthogonal_WordIndex_iterations_v1_%d', it), '.csv')), 'WriteMode', 'Append'); %write word index
        writetable(finalResults, fullfile('output', strcat(output_code, sprintf('_Orthogonal_Corrs_iterations_v1_%d', it), '.csv')), 'WriteMode', 'Append'); %write orthog correlations
        %         %% output for FUNCTION (censor if running as SCRIPT)
        %         writematrix(wordindex, fullfile(output_path,'output',strcat(output_code,sprintf('_Orthogonal_WordIndex_iterations_%d',it),'.csv')), 'WriteMode', 'Append'); %write word index
        %         writetable(finalResults, fullfile(output_path,'output',strcat(output_code,sprintf('_Orthogonal_Corrs_iterations_%d',it),'.csv')),'WriteMode', 'Append'); %write orthog correlations
        tally = 0;
    end

end

toc
