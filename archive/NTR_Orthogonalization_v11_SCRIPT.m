% Ensure the root folder of this directory is added to the MATLAB path.

%function NTR_Orthogonalization_v11(output_code) % censor if running as SCRIPT
% example for FUNCTION: NTR_Orthogonalization_v11('test')
%% Creates orthogonalized matrix comparing phoneme units, grapheme units,
% phoneme + grapheme (GP) units, letter edit distance (Levenstein), and semantic distance
%% Set Parameters
output_code='test'; % uncommment for SCRIPT
numOfIterations = 1000;
numOfWords =210;
increment=50;

output_path=pwd; % navigate to the program folder and uncenser this for running as SCRIPT

%output_path=fileparts(mfilename('fullpath'));% censor if running as SCRIPT. Get output path from the m-file location.

%% Read in info for SCRIPT (uncensor if running as script)
parameters = readtable('data/ntr_masterlist_gp.xlsx'); % creates table from masterlist input (10701 words)
jpglove=readtable('data/jpglove.csv'); %creates table from jpglove input (10701 words)
ORTable = readtable("data/ntr_masterlist_onset_rimes.xlsx"); % creates table from onset-rime input (10701 words)
ORTable = table2array(ORTable); %formats table from onset-rime input
scope = readtable("data/ntr_masterlist_scope_upd.csv");
ELP=readtable('data/ntr_masterlist_elp_with_values_upd.xlsx','Sheet','Values');% removed words repated with ed/s/ing morphemes
biphone=readtable('data/iphod_wohoms_phonprob_edit.csv');

%% Read in info for function (censor if running as SCRIPT)
% parameters = readtable(fullfile(output_path,'ntr_masterlist_gp.xlsx')); % creates table from masterlist input (10701 words)
% jpglove=readtable('jpglove.csv'); %creates table from jpglove input (10701 words)
% ORTable = readtable(fullfile(output_path,'ntr_masterlist_onset_rimes.xlsx')); % creates table from onset-rime input (10701 words)
% ORTable = table2array(ORTable); %formats table from onset-rime input
% scope = readtable(fullfile(output_path,'ntr_masterlist_scope_upd.csv'));
% ELP=readtable(fullfile(output_path,'ntr_masterlist_elp_with_values_upd.xlsx'),'Sheet','Values');% removed words repated with ed/s/ing morphemes
% biphone=readtable(fullfile(output_path,'iphod_wohoms_phonprob_edit.csv'));

%% filter homographs in parameters
[list_size,~] = size(parameters);
[~,ia,ic] = unique(parameters.string); % Unique Elements in string column
v = accumarray(ic, 1); % Tally Occurrences Of Rows in string column
parameters = parameters(ia(v==1),:); % Keep Rows for string that Only Appear Once
parameters=sortrows(parameters,'Item'); %Re-sort parameters by Item #
homographs_filtered = [num2str(list_size-size(parameters,1)),' words filtered (homographs)'];
disp(homographs_filtered)
list_size = [num2str(size(parameters,1)),' words in remaining list'];
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
[list_size,~] = size(parameters);
mask2=false;
for g = 1:list_size
    wordlength = strlength(parameters.string(g));
    if (9 >= wordlength) && (wordlength >= 3)
        mask2(g,1) = 1;
    else
        mask2(g,1) = 0;
    end
end
parameters = parameters(mask2,:);
wordlength_filtered = [num2str(list_size-size(parameters,1)),' words filtered (word length)'];
disp(wordlength_filtered)
list_size = [num2str(size(parameters,1)),' words in remaining list'];
disp(list_size)
%% Filter words missing GLOVE parameters
[list_size,~] = size(parameters);
jpglove(~ismember(jpglove.Var1,parameters.string),:)=[];
[numrowsgloveorig,~] = size(jpglove);
jpgloveindex = find(isnan(jpglove.Var2));
jpglove(jpgloveindex,:) = [];
glove_filtered = [num2str(list_size-size(jpglove,1)),' words filtered (glove)'];
disp(glove_filtered)

parameters(~ismember(parameters.string,jpglove.Var1),:)=[];
[numrowsnew, ~] = size(parameters);
list_size = [num2str(size(parameters,1)),' words in remaining list'];
disp(list_size)

%% Filter ELP by # Morphemes and POS
[list_size,~] = size(parameters);
ELP((isnan(ELP.Occurrences)),:)=[];
ELP(~ismember(ELP.Word,parameters.string),:)=[];
[numrowsELPorig,~] = size(ELP);
%nmorphindex = find(isnan(ELP.NMorph));
%ELP(nmorphindex,:) = [];
%nmorphindex2 = find(ELP.NMorph>3);
%ELP(nmorphindex2,:) = [];
% index = strfind(ELP.POS,'minor');
% index = ~cellfun(@isempty, index);
% ELP(index,:) = [];
index = strfind(ELP.POS,'encl');
index = ~cellfun(@isempty, index);
ELP(index,:) = [];
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
ELP_filtered = [num2str(list_size-size(ELP,1)),' words filtered (ELP)'];
disp(ELP_filtered)

parameters(~ismember(parameters.string,ELP.Word),:)=[];
[numrowsnew, ~] = size(parameters);
list_size = [num2str(size(parameters,1)),' words in remaining list'];
disp(list_size)
%% Filter Scope by missing values
[list_size,~] = size(parameters);
scope(~ismember(scope.Word,parameters.string),:)=[];
[numrowsscopeorig,~] = size(scope);
scopeindex = find(isnan(scope.Freq_SUBTLEXUS));
scope(scopeindex,:) = [];
scopeindex2 = find(isnan(scope.BigramF_Avg_C_Log));
scope(scopeindex2,:) = [];
scopeindex3 = find(isnan(scope.TrigramF_Avg_C_Log));
scope(scopeindex3,:) = [];
scopeindex4 = find(isnan(scope.OLD20));
scope(scopeindex4,:) = [];
scopeindex5 = find(isnan(scope.PLD20));
scope(scopeindex5,:) = [];
scopeindex6 = find(isnan(scope.Sem_N_D));
scope(scopeindex6,:) = [];
scopeindex7 = find(isnan(scope.Phonographic_N));
scope(scopeindex7,:) = [];
scope_filtered = [num2str(list_size-size(scope,1)),' words filtered (scope)'];
disp(scope_filtered)

parameters(~ismember(parameters.string,scope.Word),:)=[];
[numrowsnew, ~] = size(parameters);
list_size = [num2str(size(parameters,1)),' words in remaining list'];
disp(list_size)

%% Filter biphone probability by missing values
[list_size,~] = size(parameters);
biphone(~ismember(biphone.Word,parameters.string),:)=[];
[numrowsbiphoneorig,~] = size(biphone);
biphoneindex = find(isnan(biphone.BiphonProb));
biphone(biphoneindex,:) = [];
biphone_filtered = [num2str(list_size-size(biphone,1)),' words filtered (biphone prob)'];
disp(biphone_filtered)

parameters(~ismember(parameters.string,biphone.Word),:)=[];
[numrowsnew, ~] = size(parameters);
list_size = [num2str(size(parameters,1)),' words in remaining list'];
disp(list_size)

disp('Filtering complete');

%% Match other inputs to parameters
% parameters(~ismember(parameters.string,jpglove.Var1),:)=[];
% parameters(~ismember(parameters.string,ELP.Word),:)=[];
% parameters(~ismember(parameters.string,scope.Word),:)=[];
% parameters(~ismember(parameters.string,biphone.Word),:)=[];
jpglove(~ismember(jpglove.Var1,parameters.string),:)=[];
scope(~ismember(scope.Word,parameters.string),:)=[];
ELP(~ismember(ELP.Word,parameters.string),:)=[];
biphone(~ismember(biphone.Word,parameters.string),:)=[];
%bigpTable(~ismember(bigpTable.bigpTable1,parameters.string),:)=[];
ORTable = array2table(ORTable);
ORTable(~ismember(ORTable.ORTable1,parameters.string),:)=[];

%% Format gp, bigp, onset-rimes tables
gpTablewords = parameters{:,{'string'}}; % formats empty table for stimulus words, finds every pg* column, makes table of just gp
mask = startsWith(parameters.Properties.VariableNames, 'pg'); % determines what columns in ntr_word_masterlist begin with lowercase "pg"
Fcols = parameters(:,mask); % omits columns that do not begin with "pg"
gpTable = table2cell([gpTablewords Fcols]);% combines gpTablewords and "pg" columns and converts table to cell
%bigpTable = table2cell(bigpTable);
ORTable = table2cell(ORTable);
[numrowsgp, numcolsgp] = size(gpTable); %Get # of pg columns to use later in script
%[numrowsbigp, numcolsbigp] = size(bigpTable); %Get bipg columns#, needed?
[numrowsOR, numcolsOR] = size(ORTable);
slash = "/"; %retain slash for bigp edit distance

%%
tic
it=0;
tally=0;
for a=1:numOfIterations
    fprintf("iteration %d\n", a); % printing a start for visual cue; fprintf writes data to text file; %d start\n,a prints row for each n of iteration followed by "start"

    % Perhaps consolidate later
    r = randperm(size(parameters,1)); %randomly generating numbers within how many words there are
    newparameters = parameters(r(1:numOfWords),:); %takes N (number of words) numbers from r and grabs respective row from parameters
    randWords = newparameters.string; %grabs the selected words

    wordkey = parameters.string; %grabs all words

    wordlist= strings(numOfWords, 1); %creating the list of 200 words used in the triangles
    for i=1:numOfWords %loop respective rows from word key for first N (number of words) numbers from r
        wordlist(i)=wordkey(r(i)); %%wordlist and randWords are the same
    end
    
    wordindex(a,2:[numOfWords+1]) = transpose(wordlist(1:numOfWords)); %ERROR
    wordindex(a,1) = sprintf('%d',a); %iteration
    %     % creates triangles for pairwise distances
    %     w = pdist(randWords); % determines pairwise distance
    %     E=squareform(w); % formats distance matrix
    %     PTriangles=triu(E,-1); % grabs upper triangular part of matrix from one diagonal below the center
    tri_mask=triu(ones(numOfWords)>0,1);

    % Setting up the results table
    letterTri = zeros(numOfWords);
    semanticTri = zeros(numOfWords);
    gpTri = zeros(numOfWords);
    % bigpTri = array2table(zeros(numOfWords));
    ORTri = zeros(numOfWords);
    phonemeTri = zeros(numOfWords);
    graphemeTri = zeros(numOfWords);
    WordFreq = zeros(numOfWords);
    OrthoNeigh = zeros(numOfWords);
    PhonoNeigh = zeros(numOfWords);
    PhonoGraphNeigh = zeros(numOfWords);
    SemNeighDen = zeros(numOfWords);
    BigramFreq = zeros(numOfWords);
    TrigramF_Avg_C_Log = zeros(numOfWords);
    Biphon_Prob = zeros(numOfWords);
    Length = zeros(numOfWords);
    NSyll = zeros(numOfWords);
    GPProb= zeros(numOfWords);
    ORProb= zeros(numOfWords);

    % get uniqe list
    wp_tally=nchoosek(1:numOfWords,2);

    %finding the word pairings and their distances
    for t=1:size(wp_tally,1)
        fprintf("Number of runs: %d of %d\n", t, size(wp_tally, 1))
        % for j = 1:numOfWords
        %  for k = 1:numOfWords
        word1 = wordlist(wp_tally(t,1), 1); % grabs word 1 from wordlist
        word2 = wordlist(wp_tally(t,2), 1); % grabs word 2 from wordlist

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
        letterTri(wp_tally(t,1), wp_tally(t,2)) = editDistance(word1, word2,'SwapCost',1)/max(wLength);

        % GP EDIT DISTANCE - make function
        index1 = find(strcmp(word1, string(parameters.string)) == 1); %finding where word1 is in the unique word list
        index2 = find(strcmp(word2, string(parameters.string)) == 1); %find word2 in unique word list
        if isempty(index1) == 0 && isempty(index2) == 0 %check if string was found
            string1gp = strjoin(gpTable(index1,2:end)); % condenses all 29 pg columns from gpTable into one string for word1
            string1gp(string1gp=='+')=[]; % removes all plus signs
            word1GP = tokenizedDocument(string1gp); % converts to tokenized document (splits each unit into separate cells)

            string2gp = strjoin(gpTable(index2,2:end)); % condenses all 29 pg columns from gpTable into one string for word1
            string2gp(string2gp=='+')=[]; % removes all plus signs
            word2GP = tokenizedDocument(string2gp); % converts to tokenized document (splits each unit into separate cells)

            gpDistance = editDistance(word1GP, word2GP,'SwapCost',1); % gives edit distance for words 1 and 2 gps

           % string1gp = strtrim(string1gp); % removes leading and trailing white space from strings
           % string2gp = strtrim(string2gp);

            gpLength(1) = size(tokenDetails(word1GP),1);
            gpLength(2) = size(tokenDetails(word2GP),1);

           % wordLength(1) = count(string1gp,' ') + 1;
           % wordLength(2) = count(string2gp,' ') + 1;
            gpEditDistance = gpDistance(1,1) / max(gpLength);

            gpTri(wp_tally(t,1),  wp_tally(t,2)) = gpEditDistance;
        end

        % Onset rimes EDIT DISTANCE
        index1 = find(strcmp(word1, string(parameters.string)) == 1); %finding where word1 is in the unique word list
        index2 = find(strcmp(word2, string(parameters.string)) == 1); %find word2 in unique word list
        if isempty(index1) == 0 && isempty(index2) == 0 %check if string was found
            %string1OR = strjoin(ORTable(index1,2:end)); % condenses all 29 pg columns from gpTable into one string for word1
            for n=2:numcolsOR
                string1OR = string(ORTable(index1,n));
                string1OR = strrep(string1OR,'+','');
                string1ORfull(1,n)=string1OR;
                string2OR = string(ORTable(index2,n));
                string2OR = strrep(string2OR,'+','');
                string2ORfull(1,n)=string2OR;
            end
            string1ORtok=tokenizedDocument(string1ORfull,'TokenizeMethod','none');
            string2ORtok=tokenizedDocument(string2ORfull, 'TokenizeMethod','none');
            ORDistance = editDistance(string1ORtok, string2ORtok,'SwapCost',1);

           % string1OR = strtrim(string1OR); % removes leading and trailing white space from strings
           % string2OR = strtrim(string2OR);
            ORLength(1) = size(tokenDetails(string1ORtok),1);
            ORLength(2) = size(tokenDetails(string2ORtok),1);
            OREditDistance = ORDistance(1,1) / max(ORLength);

            ORTri(wp_tally(t,1), wp_tally(t,2)) = OREditDistance;
        end

        % SEMANTIC EDIT DISTANCE - change to use file, not variables from workspace
        word1Vars = jpglove(r(1, wp_tally(t,1)), 2:301); % grabs all glove variables for word 1
        word2Vars = jpglove(r(1, wp_tally(t,2)), 2:301); % grabs all glove variables for word 2
        semanticCorr = corrcoef(table2array(word1Vars), table2array(word2Vars)); %grabs correlation coefficient
        semanticTri(wp_tally(t,1), wp_tally(t,2)) = 1 - semanticCorr(1, 2); % records correlation

        % PHONEME UNIT AND GRAPHEME UNIT EDIT DISTANCE (separated, not P+G)
        if isempty(index1) == 0 && isempty(index2) == 0 %check if string was found
            phonemeWord1 = string; % sets up four variables with empty strings
            graphemeWord1 = string;
            phonemeWord2 = string;
            graphemeWord2 = string;
            for m = 2:numcolsgp %Read the number of columns in gpTable (don't include words column)
                pString1 = strjoin(gpTable(index1,2:end)); % check if p correct

                pgWord1 = gpTable{index1,m};
                pWord1 = extractBefore(pgWord1,"+"); %grabs just phoneme (before plus)
                phonemeWord1 = append(phonemeWord1, pWord1, ' '); % puts phonemes into one string

                pgWord1 = gpTable{index1,m};
                gWord1 = extractAfter(pgWord1,"+"); %grabs just grapheme (after plus)
                graphemeWord1 = append(graphemeWord1, gWord1, ' ');% puts graphemes into one string

                pString2 = strjoin(gpTable(index2,2:end)); % check if p correct
                pgWord2 = gpTable{index2,m};
                pWord2 = extractBefore(pgWord2,"+");
                phonemeWord2 = append(phonemeWord2, pWord2, ' ');

                pgWord2 = gpTable{index2,m};
                gWord2 = extractAfter(pgWord2,"+");
                graphemeWord2 = append(graphemeWord2, gWord2, ' ');
            end

            % creates tokenized documents
            phonemeWord1 = tokenizedDocument(phonemeWord1);
            graphemeWord1 = tokenizedDocument(graphemeWord1);
            phonemeWord2 = tokenizedDocument(phonemeWord2);
            graphemeWord2 = tokenizedDocument(graphemeWord2);

            % calculates edit distances for phonemes and graphemes
            phonemeLength(1) = size(tokenDetails(phonemeWord1),1);
            phonemeLength(2) = size(tokenDetails(phonemeWord2),1);

            graphemeLength(1) = size(tokenDetails(graphemeWord1),1);
            graphemeLength(2) = size(tokenDetails(graphemeWord2),1);

            pEditDistance = editDistance(phonemeWord1,phonemeWord2,'SwapCost',1);
            gEditDistance = editDistance(graphemeWord1,graphemeWord2,'SwapCost',1);
            phonemeEditDistance = pEditDistance(1,1) / max(phonemeLength);
            graphemeEditDistance = gEditDistance(1,1) / max(graphemeLength);

            phonemeTri(wp_tally(t,1), wp_tally(t,2)) = phonemeEditDistance; % records phonemes and graphemes edit distances
            graphemeTri(wp_tally(t,1),wp_tally(t,2)) = graphemeEditDistance;
        end

        %% SUM DISTANCES
        index1 = find(strcmp(word1, string(scope.Word)) == 1); %finding where word1 is in the unique word list
        index2 = find(strcmp(word2, string(scope.Word)) == 1);

        WordFreq(wp_tally(t,1),wp_tally(t,2)) = scope.Freq_SUBTLEXUS(index1) + scope.Freq_SUBTLEXUS(index2);
        OrthoNeigh(wp_tally(t,1),wp_tally(t,2)) = scope.OLD20(index1) + scope.OLD20(index2);
        PhonoNeigh(wp_tally(t,1),wp_tally(t,2)) = scope.PLD20(index1) + scope.PLD20(index2);
        PhonoGraphNeigh(wp_tally(t,1),wp_tally(t,2)) = scope.Phonographic_N(index1) + scope.Phonographic_N(index2);
        SemNeighDen(wp_tally(t,1),wp_tally(t,2)) = scope.Sem_N_D(index1) + scope.Sem_N_D(index2);
        BigramFreq(wp_tally(t,1),wp_tally(t,2)) = scope.BigramF_Avg_C_Log(index1) + scope.BigramF_Avg_C_Log(index2);
        TrigramF_Avg_C_Log(wp_tally(t,1),wp_tally(t,2)) = scope.TrigramF_Avg_C_Log(index1) + scope.TrigramF_Avg_C_Log(index2);
        Biphon_Prob(wp_tally(t,1),wp_tally(t,2)) = biphone.BiphonProb(index1) + biphone.BiphonProb(index2);

        GPProb(wp_tally(t,1),wp_tally(t,2)) = parameters.GP_mean(index1) + parameters.GP_mean(index2);
        ORProb(wp_tally(t,1),wp_tally(t,2)) = parameters.ORGP_mean(index1) + parameters.ORGP_mean(index2);

        %% SUBTRACTION DISTANCES
        index1 = find(strcmp(word1, string(ELP.Word)) == 1); %finding where word1 is in the unique word list
        index2 = find(strcmp(word2, string(ELP.Word)) == 1);

        Length(wp_tally(t,1),wp_tally(t,2)) = ELP.Length(index1) - ELP.Length(index2);
        NSyll(wp_tally(t,1),wp_tally(t,2)) = ELP.NSyll(index1) - ELP.NSyll(index2);
    end
    %fprintf("%d after triangle\n", a);

    % Spearman correlation comparing GP edit distance to phonological and orthographic
    GPProbToPhono(a,1) = corr(phonemeTri(tri_mask), GPProb(tri_mask),'Type','Spearman');
    GPProbToletter(a,1) = corr(letterTri(tri_mask), GPProb(tri_mask),'Type','Spearman' );
    GPProbToSem(a,1) = corr(semanticTri(tri_mask), GPProb(tri_mask),'Type','Spearman' );
    GPProbToGraph(a,1) = corr(graphemeTri(tri_mask), GPProb(tri_mask),'Type','Spearman' );
    GPProbToWordFreq(a,1) = corr( WordFreq(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToOrthoNeigh(a,1) = corr( OrthoNeigh(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToPhonoNeigh(a,1) = corr( PhonoNeigh(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToPhonoGraphNeigh(a,1) = corr( PhonoGraphNeigh(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToSemNeighDen(a,1) = corr( SemNeighDen(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToBigramFreq(a,1) = corr( BigramFreq(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToLength(a,1) = corr( Length(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToNSyll(a,1) = corr( NSyll(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToBiphon_Prob(a,1) = corr( Biphon_Prob(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToOR(a,1) = corr( ORTri(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbTogp(a,1) = corr( gpTri(tri_mask),  GPProb(tri_mask),'Type','Spearman' );
    GPProbToORProb(a,1) = corr( ORProb(tri_mask),  GPProb(tri_mask),'Type','Spearman' );

    ORProbToPhono(a,1) = corr(phonemeTri(tri_mask), ORProb(tri_mask),'Type','Spearman');
    ORProbToletter(a,1) = corr(letterTri(tri_mask), ORProb(tri_mask),'Type','Spearman' );
    ORProbToSem(a,1) = corr(semanticTri(tri_mask), ORProb(tri_mask),'Type','Spearman' );
    ORProbToGraph(a,1) = corr(graphemeTri(tri_mask), ORProb(tri_mask),'Type','Spearman' );
    ORProbToWordFreq(a,1) = corr( WordFreq(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToOrthoNeigh(a,1) = corr( OrthoNeigh(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToPhonoNeigh(a,1) = corr( PhonoNeigh(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToPhonoGraphNeigh(a,1) = corr( PhonoGraphNeigh(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToSemNeighDen(a,1) = corr( SemNeighDen(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToBigramFreq(a,1) = corr( BigramFreq(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToLength(a,1) = corr( Length(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToNSyll(a,1) = corr( NSyll(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToBiphon_Prob(a,1) = corr( Biphon_Prob(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbToOR(a,1) = corr( ORTri(tri_mask),  ORProb(tri_mask),'Type','Spearman' );
    ORProbTogp(a,1) = corr( gpTri(tri_mask),  ORProb(tri_mask),'Type','Spearman' );

    gpToPhono(a,1) = corr(phonemeTri(tri_mask), gpTri(tri_mask),'Type','Spearman');
    gpToletter(a,1) = corr(letterTri(tri_mask), gpTri(tri_mask),'Type','Spearman' );
    gpToSem(a,1) = corr(semanticTri(tri_mask), gpTri(tri_mask),'Type','Spearman' );
    gpToGraph(a,1) = corr(graphemeTri(tri_mask), gpTri(tri_mask),'Type','Spearman' );
    gpToWordFreq(a,1) = corr( WordFreq(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToOrthoNeigh(a,1) = corr( OrthoNeigh(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToPhonoNeigh(a,1) = corr( PhonoNeigh(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToPhonoGraphNeigh(a,1) = corr( PhonoGraphNeigh(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToSemNeighDen(a,1) = corr( SemNeighDen(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToBigramFreq(a,1) = corr( BigramFreq(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToLength(a,1) = corr( Length(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToNSyll(a,1) = corr( NSyll(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToBiphon_Prob(a,1) = corr( Biphon_Prob(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
    gpToOR(a,1) = corr( ORTri(tri_mask),  gpTri(tri_mask),'Type','Spearman' );

    ORToPhono(a,1) = corr(phonemeTri(tri_mask), ORTri(tri_mask),'Type','Spearman' );
    ORToletter(a,1) = corr(letterTri(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToSem(a,1) = corr(semanticTri(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToGraph(a,1) = corr(graphemeTri(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToWordFreq(a,1) = corr( WordFreq(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToOrthoNeigh(a,1) = corr( OrthoNeigh(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToPhonoNeigh(a,1) = corr( PhonoNeigh(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToPhonoGraphNeigh(a,1) = corr( PhonoGraphNeigh(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToSemNeighDen(a,1) = corr( SemNeighDen(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToBigramFreq(a,1) = corr( BigramFreq(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToLength(a,1) = corr( Length(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToNSyll(a,1) = corr( NSyll(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
    ORToBiphon_Prob(a,1) = corr( Biphon_Prob(tri_mask),  ORTri(tri_mask),'Type','Spearman' );

    WordFreqToPhono(a,1) = corr( WordFreq(tri_mask),  phonemeTri(tri_mask),'Type','Spearman' );
    WordFreqToLetter(a,1) = corr( WordFreq(tri_mask),  letterTri(tri_mask),'Type','Spearman' );
    WordFreqToSem(a,1) = corr( WordFreq(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );
    WordFreqToGraph(a,1) = corr( WordFreq(tri_mask),  graphemeTri(tri_mask),'Type','Spearman' );
    WordFreqToOrthoNeigh(a,1) = corr( WordFreq(tri_mask),  OrthoNeigh(tri_mask),'Type','Spearman' );
    WordFreqToPhonoNeigh(a,1) = corr( WordFreq(tri_mask),  PhonoNeigh(tri_mask),'Type','Spearman' );
    WordFreqToPhonoGraphNeigh(a,1) = corr( WordFreq(tri_mask),  PhonoGraphNeigh(tri_mask),'Type','Spearman' );
    WordFreqToSemNeighDen(a,1) = corr( WordFreq(tri_mask),  SemNeighDen(tri_mask),'Type','Spearman' );
    WordFreqToBigramFreq(a,1) = corr( WordFreq(tri_mask),  BigramFreq(tri_mask),'Type','Spearman' );
    WordFreqToLength(a,1) = corr( WordFreq(tri_mask),  Length(tri_mask),'Type','Spearman' );
    WordFreqToNSyll(a,1) = corr( WordFreq(tri_mask),  NSyll(tri_mask),'Type','Spearman' );
    WordFreqToBiphon_Prob(a,1) = corr( WordFreq(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );

    OrthoNeighToPhonoNeigh(a,1) = corr( OrthoNeigh(tri_mask),  PhonoNeigh(tri_mask), 'Type','Spearman' );
    OrthoNeighToPhonoGraphNeigh(a,1) = corr( OrthoNeigh(tri_mask),  PhonoGraphNeigh(tri_mask), 'Type','Spearman' );
    OrthoNeighToSemNeighDen(a,1) = corr( OrthoNeigh(tri_mask),  SemNeighDen(tri_mask), 'Type','Spearman' );
    OrthoNeighToBigramFreq(a,1) = corr( OrthoNeigh(tri_mask),  BigramFreq(tri_mask), 'Type','Spearman' );
    OrthoNeighToLength(a,1) = corr( OrthoNeigh(tri_mask),  Length(tri_mask), 'Type','Spearman' );
    OrthoNeighToNSyll(a,1) = corr( OrthoNeigh(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
    OrthoNeighToBiphon_Prob(a,1) = corr( OrthoNeigh(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
    OrthoNeighToPhono(a,1)= corr( OrthoNeigh(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
    OrthoNeighToLetter(a,1)= corr( OrthoNeigh(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
    OrthoNeighToGraph(a,1)= corr( OrthoNeigh(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );
    OrthoNeighToSem(a,1)= corr( OrthoNeigh(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' ); 

    PhonoNeighToPhonoGraphNeigh(a,1) = corr( PhonoNeigh(tri_mask),  PhonoGraphNeigh(tri_mask), 'Type','Spearman' );
    PhonoNeighToSemNeighDen(a,1) = corr( PhonoNeigh(tri_mask),  SemNeighDen(tri_mask), 'Type','Spearman' );
    PhonoNeighToBigramFreq(a,1) = corr( PhonoNeigh(tri_mask),  BigramFreq(tri_mask), 'Type','Spearman' );
    PhonoNeighToLength(a,1) = corr( PhonoNeigh(tri_mask),  Length(tri_mask), 'Type','Spearman' );
    PhonoNeighToNSyll(a,1) = corr( PhonoNeigh(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
    PhonoNeighToBiphon_Prob(a,1) = corr( PhonoNeigh(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
    PhonoNeighToSem(a,1)= corr( PhonoNeigh(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );
    PhonoNeighToPhono(a,1)= corr( PhonoNeigh(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
    PhonoNeighToLetter(a,1)= corr( PhonoNeigh(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
    PhonoNeighToGraph(a,1)= corr( PhonoNeigh(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );

    PhonoGraphNeighToSemNeighDen(a,1) = corr( PhonoGraphNeigh(tri_mask),  SemNeighDen(tri_mask), 'Type','Spearman' );
    PhonoGraphNeighToBigramFreq(a,1) = corr( PhonoGraphNeigh(tri_mask),  BigramFreq(tri_mask), 'Type','Spearman' );
    PhonoGraphNeighToLength(a,1) = corr( PhonoGraphNeigh(tri_mask),  Length(tri_mask), 'Type','Spearman' );
    PhonoGraphNeighToNSyll(a,1) = corr( PhonoGraphNeigh(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
    PhonoGraphNeighToBiphon_Prob(a,1) = corr( PhonoGraphNeigh(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
    Phonographic_ToSem(a,1)= corr( PhonoGraphNeigh(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );
    PhonoGraphNeighToPhono(a,1)= corr( PhonoGraphNeigh(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
    PhonoGraphNeighToLetter(a,1)= corr( PhonoGraphNeigh(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
    PhonoGraphNeighTo_Graph(a,1)= corr( PhonoGraphNeigh(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );

    SemNeighDenToBigramFreq(a,1) = corr( SemNeighDen(tri_mask),  BigramFreq(tri_mask), 'Type','Spearman' );
    SemNeighDenToLength(a,1) = corr( SemNeighDen(tri_mask),  Length(tri_mask), 'Type','Spearman' );
    SemNeighDenToNSyll(a,1) = corr( SemNeighDen(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
    SemNeighDenToBiphon_Prob(a,1) = corr( SemNeighDen(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
    SemNeighDenToSem(a,1)= corr( SemNeighDen(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );
    SemNeighDenToPhono(a,1)= corr( SemNeighDen(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
    SemNeighDenToLetter(a,1)= corr( SemNeighDen(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
    SemNeighDenToGraph(a,1)= corr( SemNeighDen(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );

    BigramFreqToLength(a,1) = corr( BigramFreq(tri_mask),  Length(tri_mask), 'Type','Spearman' );
    BigramFreqToNSyll(a,1) = corr( BigramFreq(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
    BigramFreqToBiphon_Prob(a,1) = corr( BigramFreq(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
    BigramFreqToSem(a,1)= corr( BigramFreq(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );
    BigramFreqToPhono(a,1)= corr( BigramFreq(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
    BigramFreqToLetter(a,1)= corr( BigramFreq(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
    BigramFreqToGraph(a,1)= corr( BigramFreq(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );

    graphToletter(a,1) = corr( graphemeTri(tri_mask),  letterTri(tri_mask),'Type','Spearman' );
    graphToPhono(a,1) = corr( graphemeTri(tri_mask),  phonemeTri(tri_mask),'Type','Spearman' );
    graphToLength(a,1)= corr( graphemeTri(tri_mask),  Length(tri_mask),'Type','Spearman' );
    graphToNSyll(a,1)= corr( graphemeTri(tri_mask),  NSyll(tri_mask),'Type','Spearman' );
    graphToBiphon_Prob(a,1)= corr( graphemeTri(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );
    graphToSem(a,1)= corr( graphemeTri(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );

    phonoToletter(a,1) = corr( phonemeTri(tri_mask),  letterTri(tri_mask),'Type','Spearman' );
    phonoToLength(a,1)= corr( phonemeTri(tri_mask),  Length(tri_mask),'Type','Spearman' );
    phonoToNSyll(a,1)= corr( phonemeTri(tri_mask),  NSyll(tri_mask),'Type','Spearman' );
    phonoToBiphon_Prob(a,1)= corr( phonemeTri(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );
    phonoToSem(a,1)= corr( phonemeTri(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );

    LetterToLength(a,1)= corr( letterTri(tri_mask),  Length(tri_mask),'Type','Spearman' );
    LetterToNSyll(a,1)= corr( letterTri(tri_mask),  NSyll(tri_mask),'Type','Spearman' );
    LetterToBiphon_Prob(a,1)= corr( letterTri(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );
    LetterToSem(a,1)= corr( letterTri(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );

    LengthToBiphon_Prob(a,1)= corr( Length(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );
    LengthToSem(a,1)= corr( Length(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );

    NSyllToBiphon_Prob(a,1) = corr( NSyll(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
    NSyllToSem(a,1)= corr( NSyll(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );

    Biphon_ProbToSem_Corr(a,1) = corr( Biphon_Prob(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );

    NSyllToLen(a,1) = corr( NSyll(tri_mask),  Length(tri_mask), 'Type','Spearman' );

    %Example of getting 1 stdev for the whole correlation matrix
    w=0; %(w=0 is n-1, w=1 is n for stdev)
    std_phonemeTri(a,1) =std(phonemeTri(tri_mask),w,"all");
    std_gpTri(a,1) =std(gpTri(tri_mask),w,"all");
    std_letterTri(a,1) =std(letterTri(tri_mask),w,"all");
    std_semanticTri(a,1) =std(semanticTri(tri_mask),w,"all");
    std_graphemeTri(a,1) =std(graphemeTri(tri_mask),w,"all");
    std_ORTri(a,1) =std(ORTri(tri_mask),w,"all");
    std_WordFreq(a,1) =std(WordFreq(tri_mask),w,"all");
    std_OrthoNeigh(a,1) =std(OrthoNeigh(tri_mask),w,"all");
    std_PhonoNeigh(a,1) =std(PhonoNeigh(tri_mask),w,"all");
    std_PhonoGraphNeigh(a,1) =std(PhonoGraphNeigh(tri_mask),w,"all");
    std_SemNeighDen(a,1) =std(SemNeighDen(tri_mask),w,"all");
    std_BigramFreq(a,1) =std(BigramFreq(tri_mask),w,"all");
    std_Length(a,1) =std(Length(tri_mask),w,"all");
    std_NSyll(a,1) =std(NSyll(tri_mask),w,"all");
    std_Biphon_Prob(a,1) =std(Biphon_Prob(tri_mask),w,"all");
  
    % print out results every 1000 iterations
    it=it+1;
    tally=tally+1;
    if tally==increment
        %% writes everything in excel spreadsheets
        iteration(1:it,1)=(1:it)';
        %wordindex(:,1) = iteration; %Number rows 1 through number of iterations for word index
        % set up table based on variables

        finalResults=table(iteration,NSyllToLen,gpToPhono	,gpToletter	,gpToSem,gpToGraph,gpToWordFreq,gpToOrthoNeigh,gpToPhonoNeigh,gpToPhonoGraphNeigh,gpToSemNeighDen,  ...
            gpToBigramFreq,gpToLength,gpToNSyll,gpToBiphon_Prob,gpToOR,ORToPhono,ORToletter,ORToSem,ORToGraph,ORToWordFreq,  ...
            ORToOrthoNeigh,ORToPhonoNeigh,ORToPhonoGraphNeigh,ORToSemNeighDen,ORToBigramFreq,ORToLength,ORToNSyll,ORToBiphon_Prob,WordFreqToPhono,  ...
            WordFreqToLetter,WordFreqToSem,WordFreqToGraph,WordFreqToOrthoNeigh,WordFreqToPhonoNeigh,WordFreqToPhonoGraphNeigh,  ...
            WordFreqToSemNeighDen,WordFreqToBigramFreq,WordFreqToLength,WordFreqToNSyll,  ...
            WordFreqToBiphon_Prob,OrthoNeighToPhonoNeigh,OrthoNeighToPhonoGraphNeigh,OrthoNeighToSemNeighDen,OrthoNeighToBigramFreq,OrthoNeighToLength,OrthoNeighToNSyll,  ...
            OrthoNeighToBiphon_Prob,OrthoNeighToPhono,OrthoNeighToLetter,OrthoNeighToGraph,OrthoNeighToSem,PhonoNeighToPhonoGraphNeigh,PhonoNeighToSemNeighDen,  ...
            PhonoNeighToBigramFreq,PhonoNeighToLength,PhonoNeighToNSyll,PhonoNeighToBiphon_Prob,PhonoNeighToSem,PhonoNeighToPhono,PhonoNeighToLetter,PhonoNeighToGraph,  ...
            PhonoGraphNeighToSemNeighDen,PhonoGraphNeighToBigramFreq,PhonoGraphNeighToLength,PhonoGraphNeighToNSyll,PhonoGraphNeighToBiphon_Prob,  ...
            Phonographic_ToSem,PhonoGraphNeighToPhono,PhonoGraphNeighToLetter,PhonoGraphNeighTo_Graph,SemNeighDenToBigramFreq,SemNeighDenToLength,SemNeighDenToNSyll,  ...
            SemNeighDenToBiphon_Prob,SemNeighDenToSem,SemNeighDenToPhono,SemNeighDenToLetter,SemNeighDenToGraph,BigramFreqToLength,  ...
            BigramFreqToNSyll,BigramFreqToBiphon_Prob,BigramFreqToSem,  ...
            BigramFreqToPhono,BigramFreqToLetter,BigramFreqToGraph,graphToletter,graphToPhono,graphToLength,graphToNSyll,  ...
            graphToBiphon_Prob,graphToSem,phonoToletter,phonoToLength,phonoToNSyll,phonoToBiphon_Prob,phonoToSem,LetterToLength,  ...
            LetterToNSyll,LetterToBiphon_Prob,LetterToSem,LengthToBiphon_Prob,LengthToSem,NSyllToBiphon_Prob,NSyllToSem,Biphon_ProbToSem_Corr, ...
            GPProbToPhono,GPProbToletter,GPProbToSem,GPProbToGraph,GPProbToWordFreq	,GPProbToOrthoNeigh,GPProbToPhonoNeigh, ...
            GPProbToPhonoGraphNeigh	,GPProbToSemNeighDen	,GPProbToBigramFreq	,GPProbToLength	,GPProbToNSyll	,GPProbToBiphon_Prob, ...
            GPProbToOR	,GPProbTogp	,GPProbToORProb	,ORProbToPhono	,ORProbToletter	,ORProbToSem	,ORProbToGraph	,ORProbToWordFreq	, ...
            ORProbToOrthoNeigh	,ORProbToPhonoNeigh	,ORProbToPhonoGraphNeigh	,ORProbToSemNeighDen	,ORProbToBigramFreq	, ...
            ORProbToLength	,ORProbToNSyll	,ORProbToBiphon_Prob	,ORProbToOR	,ORProbTogp, ...
            std_phonemeTri,std_gpTri,std_letterTri, std_semanticTri,std_graphemeTri,std_ORTri,std_WordFreq, std_OrthoNeigh, ...
            std_PhonoNeigh,std_PhonoGraphNeigh,std_SemNeighDen, std_BigramFreq,std_Length,std_NSyll, std_Biphon_Prob);

        %% output for SCRIPT uncensor if running as SCRIPT
         writematrix(wordindex, fullfile('output',strcat(output_code,sprintf('_Orthogonal_WordIndex_iterations_%d',it),'.csv')), 'WriteMode', 'Append'); %write word index
         writetable(finalResults, fullfile('output',strcat(output_code,sprintf('_Orthogonal_Corrs_iterations_%d',it),'.csv')),'WriteMode', 'Append'); %write orthog correlations
        %% output for FUNCTION (censor if running as SCRIPT)
        %writematrix(wordindex, fullfile(output_path,'output',strcat(output_code,sprintf('_Orthogonal_WordIndex_iterations_%d',it),'.csv')), 'WriteMode', 'Append'); %write word index
        %writetable(finalResults, fullfile(output_path,'output',strcat(output_code,sprintf('_Orthogonal_Corrs_iterations_%d',it),'.csv')),'WriteMode', 'Append'); %write orthog correlations
        tally=0;
    end
end
toc
