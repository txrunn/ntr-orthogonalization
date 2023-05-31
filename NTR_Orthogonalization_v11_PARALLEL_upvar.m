function NTR_Orthogonalization_v11_PARALLEL(groups,core_num) % censor if running as SCRIPT
% example for FUNCTION: NTR_Orthogonalization_v11(groups)
%% Creates orthogonalized matrix comparing phoneme units, grapheme units,
% phoneme + grapheme (GP) units, letter edit distance (Levenstein), and semantic distance
%parpool;
disp(datetime)

%% Set Parameters
unique_code='jp_zaratan_set4';
numOfIterations = 20;
numOfWords =160;

%% set up parfor loop
parfor pa=1:groups% use up to 100 cores
    %maxNumCompThreads(40);
    %fprintf("iteration %d\n", a); % printing a start for visual cue;
    %fprintf writes data to text file; %d start\n,a prints row for each n
    %of iteration followed by "start"
    %output_path=pwd; % navigate to the program folder and uncenser this for running as SCRIPT
    output_code=sprintf('%s_core%d_group%d',unique_code,core_num,pa); % uncommment for SCRIPT
    output_path=fileparts(mfilename('fullpath'));% censor if running as SCRIPT. Get output path from the m-file location.

    %% Read in info for SCRIPT (uncensor if running as script)
    % parameters = readtable('ntr_masterlist_gp.xlsx'); % creates table from masterlist input (10701 words)
    % jpglove=readtable('jpglove.csv'); %creates table from jpglove input (10701 words)
    % ORTable = readtable("ntr_masterlist_onset_rimes.xlsx"); % creates table from onset-rime input (10701 words)
    % ORTable = table2array(ORTable); %formats table from onset-rime input
    % scope = readtable("ntr_masterlist_scope_upd.csv");
    % ELP=readtable('ntr_masterlist_elp_with_values_upd.xlsx','Sheet','Values');% removed words repated with ed/s/ing morphemes
    % biphone=readtable('iphod_wohoms_phonprob_edit.csv');

    %% Read in info for function (censor if running as SCRIPT)
    parameters = readtable(fullfile(output_path,'ntr_masterlist_gp.xlsx')); % creates table from masterlist input (10701 words)
    jpglove=readtable('jpglove.csv'); %creates table from jpglove input (10701 words)
    ORTable = readtable(fullfile(output_path,'ntr_masterlist_onset_rimes.xlsx')); % creates table from onset-rime input (10701 words)
    ORTable = table2array(ORTable); %formats table from onset-rime input
    scope = readtable(fullfile(output_path,'ntr_masterlist_scope_upd.csv'));
    ELP=readtable(fullfile(output_path,'ntr_masterlist_elp_with_values_upd.xlsx'),'Sheet','Values');% removed words repated with ed/s/ing morphemes
    biphone=readtable(fullfile(output_path,'iphod_wohoms_phonprob_edit.csv'));

    %% filter homographs in parameters
    [list_size,~] = size(parameters);
    [~,ia,ic] = unique(parameters.string); % Unique Elements in string column
    v = accumarray(ic, 1); % Tally Occurrences Of Rows in string column
    parameters = parameters(ia(v==1),:); % Keep Rows for string that Only Appear Once
    parameters=sortrows(parameters,'Item'); %Re-sort parameters by Item #
    homographs_filtered = [num2str(list_size-size(parameters,1)),' words filtered (homographs)'];
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


    %% set variables
    GPProbToPhono = zeros(numOfIterations,1);
    GPProbToLetter = zeros(numOfIterations,1);
    GPProbToSem = zeros(numOfIterations,1);
    GPProbToGraph = zeros(numOfIterations,1);
    GPProbToWordFreq = zeros(numOfIterations,1);
    GPProbToOrthoNeigh = zeros(numOfIterations,1);
    GPProbToPhonoNeigh = zeros(numOfIterations,1);
    GPProbToPhonoGraphNeigh =zeros(numOfIterations,1);
    GPProbToSemNeighDen = zeros(numOfIterations,1);
    GPProbToBigramFreq =zeros(numOfIterations,1);
    GPProbToNLetters = zeros(numOfIterations,1);
    GPProbToNSyll = zeros(numOfIterations,1);
    GPProbToBiphonProb = zeros(numOfIterations,1);
    GPProbToOR = zeros(numOfIterations,1);
    GPProbToGP = zeros(numOfIterations,1);
    GPProbToORProb =zeros(numOfIterations,1);
    ORProbToPhono =zeros(numOfIterations,1);
    ORProbToLetter =zeros(numOfIterations,1);
    ORProbToSem =zeros(numOfIterations,1);
    ORProbToGraph =zeros(numOfIterations,1);
    ORProbToWordFreq =zeros(numOfIterations,1);
    ORProbToOrthoNeigh =zeros(numOfIterations,1);
    ORProbToPhonoNeigh =zeros(numOfIterations,1);
    ORProbToPhonoGraphNeigh =zeros(numOfIterations,1);
    ORProbToSemNeighDen =zeros(numOfIterations,1);
    ORprobToBigramFreq =zeros(numOfIterations,1);
    ORProbToNLetters =zeros(numOfIterations,1);
    ORProbToNSyll =zeros(numOfIterations,1);
    ORProbToBiphonProb =zeros(numOfIterations,1);
    ORProbToOR =zeros(numOfIterations,1);
    ORProbToGP =zeros(numOfIterations,1);
    GPToPhono =zeros(numOfIterations,1);
    GPToLetter =zeros(numOfIterations,1);
    GPToSem = zeros(numOfIterations,1);
    GPToGraph =zeros(numOfIterations,1);
    GPToWordFreq =zeros(numOfIterations,1);
    GPToOrthoNeigh =zeros(numOfIterations,1);
    GPToPhonoNeigh =zeros(numOfIterations,1);
    GPToPhonoGraphNeigh =zeros(numOfIterations,1);
    GPToSemNeighDen =zeros(numOfIterations,1);
    GPToBigramFreq =zeros(numOfIterations,1);
    GPToNLetters = zeros(numOfIterations,1);
    GPToNSyll =zeros(numOfIterations,1);
    GPToBiphonProb =zeros(numOfIterations,1);
    GPToOR =zeros(numOfIterations,1);
    ORToPhono =zeros(numOfIterations,1);
    ORToLetter =zeros(numOfIterations,1);
    ORToSem =zeros(numOfIterations,1);
    ORToGraph =zeros(numOfIterations,1);
    ORToWordFreq =zeros(numOfIterations,1);
    ORToOrthoNeigh =zeros(numOfIterations,1);
    ORToPhonoNeigh =zeros(numOfIterations,1);
    ORToPhonoGrpahNeigh = zeros(numOfIterations,1);
    ORToSemNeighDen =zeros(numOfIterations,1);
    ORToBigramFreq =zeros(numOfIterations,1);
    ORToNLetters =zeros(numOfIterations,1);
    ORToNSyll = zeros(numOfIterations,1);
    ORToBiphonProb =zeros(numOfIterations,1);
    WordFreqToPhono =zeros(numOfIterations,1);
    WordFreqToLetter =zeros(numOfIterations,1);
    WordFreqToSem =zeros(numOfIterations,1);
    WordFreqToGraph =zeros(numOfIterations,1);
    WordFreqToOrthoNeigh =zeros(numOfIterations,1);
    WordFreqToPhonoNeigh =zeros(numOfIterations,1);
    WordFreqToPhonoGraphNegih =zeros(numOfIterations,1);
    WordFreqToSemNeighDen =zeros(numOfIterations,1);
    WordFreqToBigramFreq =zeros(numOfIterations,1);
    WordFreqToNLetters =zeros(numOfIterations,1);
    WordFreqToNSyll =zeros(numOfIterations,1);
    WordFreqToBiphonProb =zeros(numOfIterations,1);
    OrthoNeighToPhonoNeigh =zeros(numOfIterations,1);
    OrthoNeighToPhonoGraphNeigh =zeros(numOfIterations,1);
    OrthoNeighToSemNeighDen =zeros(numOfIterations,1);
    OrthoNeighToBigramFreq = zeros(numOfIterations,1);
    OrthoNeighToNLetters =zeros(numOfIterations,1);
    OrthoNeighToNSyll = zeros(numOfIterations,1);
    OrthoNeighToBiphonProb =zeros(numOfIterations,1);
    OrthoNeighToPhono=zeros(numOfIterations,1);
    OrthoNeighToLetter=zeros(numOfIterations,1);
    OrthoNeighToGraph=zeros(numOfIterations,1);
    OrthoNeighToSem=zeros(numOfIterations,1);
    PhonoNeighToPhonoGraphNeigh =zeros(numOfIterations,1);
    PhonoNeighToSemNeighDen =zeros(numOfIterations,1);
    PhonoNeighToBigramFreq =zeros(numOfIterations,1);
    PhonoNeighToNLetters =zeros(numOfIterations,1);
    PhonoNeighToNSyll = zeros(numOfIterations,1);
    PhonoNeighToBiphonProb =zeros(numOfIterations,1);
    PhonoNeighToSem= zeros(numOfIterations,1);
    PhonoNeighToPhono=zeros(numOfIterations,1);
    PhonoNeighToLetter=zeros(numOfIterations,1);
    PhonoNeighToGraph=zeros(numOfIterations,1);
    PhonoGraphNeighToSemNeighDen =zeros(numOfIterations,1);
    PhonoGraphNeighToBigramFreq =zeros(numOfIterations,1);
    PhonoGraphNeighToNLetters =zeros(numOfIterations,1);
    PhonoGraphNeighToNSyll =zeros(numOfIterations,1);
    PhonoGraphNeighToBiphonProb =zeros(numOfIterations,1);
    PhonoGraphNeighToSem=zeros(numOfIterations,1);
    PhonoGraphNeighToPhono=zeros(numOfIterations,1);
    PhonoGraphNeighToLetter=zeros(numOfIterations,1);
    PhonoGraphNeighToGraph=zeros(numOfIterations,1);
    SemNeighDenToBigramFreq =zeros(numOfIterations,1);
    SemNeighDenToNLetters =zeros(numOfIterations,1);
    SemNeighDenToNSyll =zeros(numOfIterations,1);
    SemNeighDenToBiphonProb =zeros(numOfIterations,1);
    SemNeighDenToSem=zeros(numOfIterations,1);
    SemNeighDenToPhono=zeros(numOfIterations,1);
    SemNeighDenToLetter=zeros(numOfIterations,1);
    SemNeighDenToGraph=zeros(numOfIterations,1);
    BigramFreqToNLetters =zeros(numOfIterations,1);
    BigramFreqToNSyll =zeros(numOfIterations,1);
    BigramFreqToBiphonProb =zeros(numOfIterations,1);
    BigramFreqToSem=zeros(numOfIterations,1);
    BigramFreqToPhono=zeros(numOfIterations,1);
    BigramFreqToLetter=zeros(numOfIterations,1);
    BigramFreqToGraph=zeros(numOfIterations,1);
    GraphToLetter =zeros(numOfIterations,1);
    GraphToPhono =zeros(numOfIterations,1);
    GraphToNLetters=zeros(numOfIterations,1);
    GraphToNSyll=zeros(numOfIterations,1);
    GraphToBiphonProb=zeros(numOfIterations,1);
    GraphToSem=zeros(numOfIterations,1);
    PhonoToLetter =zeros(numOfIterations,1);
    PhonoToNLetters=zeros(numOfIterations,1);
    PhonoToNSyll=zeros(numOfIterations,1);
    PhonoToBiphonProb=zeros(numOfIterations,1);
    PhonoToSem=zeros(numOfIterations,1);
    LetterToNLetters=zeros(numOfIterations,1);
    LetterToNSyll=zeros(numOfIterations,1);
    LetterToBiphonProb=zeros(numOfIterations,1);
    LetterToSem=zeros(numOfIterations,1);
    NLettersToBiphonProb=zeros(numOfIterations,1);
    NLettersToSem=zeros(numOfIterations,1);
    NSyllToBiphonProb =zeros(numOfIterations,1);
    NSyllToSem=zeros(numOfIterations,1);
    BiphonProbToSem =zeros(numOfIterations,1);
    NSyllToNLetters =zeros(numOfIterations,1);

    % standard deviations
    std_phonemeTri=zeros(numOfIterations,1);
    std_gpTri=zeros(numOfIterations,1);
    std_letterTri=zeros(numOfIterations,1);
    std_semanticTri=zeros(numOfIterations,1);
    std_graphemeTri=zeros(numOfIterations,1);
    std_ORTri=zeros(numOfIterations,1);
    std_Freq_SUBTLEXUS=zeros(numOfIterations,1);
    std_OLD20sum=zeros(numOfIterations,1);
    std_PLD20sum=zeros(numOfIterations,1);
    std_Phonographic_N=zeros(numOfIterations,1);
    std_Sem_N_D=zeros(numOfIterations,1);
    std_BigramF_Avg_C_Log=zeros(numOfIterations,1);
    std_Length=zeros(numOfIterations,1);
    std_NSyll =zeros(numOfIterations,1);
    std_Biphon_Prob=zeros(numOfIterations,1);


   % graphemeLength=zeros(1,2);
   % phonemeLength=zeros(1,2);
   % ORLength=zeros(1,2);
   % wLength=zeros(1,2);
  %  gpLength=zeros(1,2);

    wordindex=strings(numOfIterations,[numOfWords+1]);

    string2OR="+";

    iteration=zeros(numOfIterations,1);
    %% Start iterations
    it=0;
    tally=0;
    for a=1:numOfIterations
        iteration(a,1)=a;
        %fprintf("iteration %d\n", a); % printing a start for visual cue; fprintf writes data to text file; %d start\n,a prints row for each n of iteration followed by "start"

        % Perhaps consolidate later
        rng('shuffle');
        r = randperm(size(parameters,1))
        %randomly generating numbers within how many words there are
        newparameters = parameters(r(1:numOfWords),:); %takes N (number of words) numbers from r and grabs respective row from parameters
        randWords = newparameters.string; %grabs the selected words
    
        wordkey = parameters.string; %grabs all words

        wordlist= strings(numOfWords, 1); %creating the list of 200 words used in the triangles
        for i=1:numOfWords %loop respective rows from word key for first N (number of words) numbers from r
            wordlist(i)=wordkey(r(i)); %%wordlist and randWords are the same
        end

       % wordindex(a,1) = sprintf('%d',a); %iteration
       % temp=transpose(wordlist(1:numOfWords));
     %   wordindex(a,2:end) = transpose(wordlist(1:numOfWords)); %ERROR

        wordindex(a,:) =[sprintf('%d',a), transpose(wordlist(1:numOfWords))]
       % wordindex(a,2:[numOfWords+1]) = transpose(wordlist(1:numOfWords)); %ERROR
        
        %     % creates triangles for pairwise distances
        %     w = pdist(randWords); % determines pairwise distance
        %     E=squareform(w); % formats distance matrix
        %     PTriangles=triu(E,-1); % grabs upper triangular part of matrix from one diagonal below the center
        tri_mask=triu(ones(numOfWords)>0,1);

        % Setting up the results matrices
        LetterTri = zeros(numOfWords);
        PhonemeTri = zeros(numOfWords);
        GraphemeMat = zeros(numOfWords);
        SemanticMat = zeros(numOfWords);
        GPMat = zeros(numOfWords);
        % bigpMat = array2table(zeros(numOfWords));
        ORMat = zeros(numOfWords);
        WordFreqMat = zeros(numOfWords);
        OrthNeighborsMat = zeros(numOfWords);
        PhonNeighborsMat = zeros(numOfWords);
        PhonographicNumMat = zeros(numOfWords);
        SemNeighborsMat = zeros(numOfWords);
        BigramF_Avg_C_LogMat = zeros(numOfWords);
        TrigramF_Avg_C_LogMat = zeros(numOfWords);
        Biphon_ProbMat = zeros(numOfWords);
        LengthMat = zeros(numOfWords);
        NSyllMat = zeros(numOfWords);
        GPprobMat = zeros(numOfWords);
        ORprobMat = zeros(numOfWords);

        % get uniqe list
        wp_tally=nchoosek(1:numOfWords,2);

        %finding the word pairings and their distances
        for t=1:size(wp_tally,1)
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
           % wLength(1) = strlength(word1);
           % wLength(2) = strlength(word2);
            letterTri(wp_tally(t,1), wp_tally(t,2)) = editDistance(word1, word2,'SwapCost',1)/max([strlength(word1),strlength(word2)] );

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

               % gpLength(1) = size(tokenDetails(word1GP),1);
                %gpLength(2) = size(tokenDetails(word2GP),1);

                % wordLength(1) = count(string1gp,' ') + 1;
                % wordLength(2) = count(string2gp,' ') + 1;
                gpEditDistance = gpDistance(1,1) / max([size(tokenDetails(word1GP),1), size(tokenDetails(word2GP),1)]);

                gpTri(wp_tally(t,1),  wp_tally(t,2)) = gpEditDistance;
            end

            % Onset rimes EDIT DISTANCE
            index1 = find(strcmp(word1, string(parameters.string)) == 1); %finding where word1 is in the unique word list
            index2 = find(strcmp(word2, string(parameters.string)) == 1); %find word2 in unique word list
            if isempty(index1) == 0 && isempty(index2) == 0 %check if string was found
                %string1OR = strjoin(ORTable(index1,2:end)); % condenses all 29 pg columns from gpTable into one string for word1

                string1ORfull=strings(1,numcolsOR);
                string2ORfull=strings(1,numcolsOR);
                string1OR=strings(1);
                string2OR=strings(1);
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
               % ORLength(1) = size(tokenDetails(string1ORtok),1);
               % ORLength(2) = size(tokenDetails(string2ORtok),1);
                OREditDistance = ORDistance(1,1) / max([size(tokenDetails(string1ORtok),1), size(tokenDetails(string2ORtok),1)]);

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
                %phonemeLength(1) = size(tokenDetails(phonemeWord1),1);
               % phonemeLength(2) = size(tokenDetails(phonemeWord2),1);

               %graphemeLength(1) = size(tokenDetails(graphemeWord1),1);
              %  graphemeLength(2) = size(tokenDetails(graphemeWord2),1);

                pEditDistance = editDistance(phonemeWord1,phonemeWord2,'SwapCost',1);
                gEditDistance = editDistance(graphemeWord1,graphemeWord2,'SwapCost',1);
                phonemeEditDistance = pEditDistance(1,1) / max([size(tokenDetails(phonemeWord1),1),size(tokenDetails(phonemeWord2),1)]);
                graphemeEditDistance = gEditDistance(1,1) / max([size(tokenDetails(graphemeWord1),1), size(tokenDetails(graphemeWord2),1)]);

                phonemeTri(wp_tally(t,1), wp_tally(t,2)) = phonemeEditDistance; % records phonemes and graphemes edit distances
                graphemeTri(wp_tally(t,1),wp_tally(t,2)) = graphemeEditDistance;
            end

            %% SUM DISTANCES
            index1 = find(strcmp(word1, string(scope.Word)) == 1); %finding where word1 is in the unique word list
            index2 = find(strcmp(word2, string(scope.Word)) == 1);

            Freq_SUBTLEXUS(wp_tally(t,1),wp_tally(t,2)) = scope.Freq_SUBTLEXUS(index1) + scope.Freq_SUBTLEXUS(index2);
            OLD20sum(wp_tally(t,1),wp_tally(t,2)) = scope.OLD20(index1) + scope.OLD20(index2);
            PLD20sum(wp_tally(t,1),wp_tally(t,2)) = scope.PLD20(index1) + scope.PLD20(index2);
            Phonographic_N(wp_tally(t,1),wp_tally(t,2)) = scope.Phonographic_N(index1) + scope.Phonographic_N(index2);
            Sem_N_D(wp_tally(t,1),wp_tally(t,2)) = scope.Sem_N_D(index1) + scope.Sem_N_D(index2);
            BigramF_Avg_C_Log(wp_tally(t,1),wp_tally(t,2)) = scope.BigramF_Avg_C_Log(index1) + scope.BigramF_Avg_C_Log(index2);
            TrigramF_Avg_C_Log(wp_tally(t,1),wp_tally(t,2)) = scope.TrigramF_Avg_C_Log(index1) + scope.TrigramF_Avg_C_Log(index2);
            Biphon_Prob(wp_tally(t,1),wp_tally(t,2)) = biphone.BiphonProb(index1) + biphone.BiphonProb(index2);

            GPprob(wp_tally(t,1),wp_tally(t,2)) = parameters.GP_mean(index1) + parameters.GP_mean(index2);
            ORprob(wp_tally(t,1),wp_tally(t,2)) = parameters.ORGP_mean(index1) + parameters.ORGP_mean(index2);

            %% SUBTRACTION DISTANCES
            index1 = find(strcmp(word1, string(ELP.Word)) == 1); %finding where word1 is in the unique word list
            index2 = find(strcmp(word2, string(ELP.Word)) == 1);

            Length(wp_tally(t,1),wp_tally(t,2)) = ELP.Length(index1) - ELP.Length(index2);
            NSyll(wp_tally(t,1),wp_tally(t,2)) = ELP.NSyll(index1) - ELP.NSyll(index2);
        end
        %fprintf("%d after triangle\n", a);

        % Spearman correlation comparing GP edit distance to phonological and orthographic
        GPProbToPhono(a,1) = corr(phonemeTri(tri_mask), GPprob(tri_mask),'Type','Spearman');
        GPProbToLetter(a,1) = corr(letterTri(tri_mask), GPprob(tri_mask),'Type','Spearman' );
        GPProbToSem(a,1) = corr(semanticTri(tri_mask), GPprob(tri_mask),'Type','Spearman' );
        GPProbToGraph(a,1) = corr(graphemeTri(tri_mask), GPprob(tri_mask),'Type','Spearman' );
        GPProbToWordFreq(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToOrthoNeigh(a,1) = corr( OLD20sum(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToPhonoNeigh(a,1) = corr( PLD20sum(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToPhonoGraphNeigh(a,1) = corr( Phonographic_N(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToSemNeighDen(a,1) = corr( Sem_N_D(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToBigramFreq(a,1) = corr( BigramF_Avg_C_Log(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToNLetters(a,1) = corr( Length(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPprobToNSyll(a,1) = corr( NSyll(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToBiphonProb(a,1) = corr( Biphon_Prob(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToOR(a,1) = corr( ORTri(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToGP(a,1) = corr( gpTri(tri_mask),  GPprob(tri_mask),'Type','Spearman' );
        GPProbToORProb(a,1) = corr( ORprob(tri_mask),  GPprob(tri_mask),'Type','Spearman' );

        ORProbToPhono(a,1) = corr(phonemeTri(tri_mask), ORprob(tri_mask),'Type','Spearman');
        ORProbToLetter(a,1) = corr(letterTri(tri_mask), ORprob(tri_mask),'Type','Spearman' );
        ORProbToSem(a,1) = corr(semanticTri(tri_mask), ORprob(tri_mask),'Type','Spearman' );
        ORProbToGraph(a,1) = corr(graphemeTri(tri_mask), ORprob(tri_mask),'Type','Spearman' );
        ORProbToWordFreq(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToOrthoNeigh(a,1) = corr( OLD20sum(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToPhonoNeigh(a,1) = corr( PLD20sum(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToPhonoGraphNeigh(a,1) = corr( Phonographic_N(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToSemNeighDen(a,1) = corr( Sem_N_D(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORprobToBigramFreq(a,1) = corr( BigramF_Avg_C_Log(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToNLetters(a,1) = corr( Length(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToNSyll(a,1) = corr( NSyll(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToBiphonProb(a,1) = corr( Biphon_Prob(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToOR(a,1) = corr( ORTri(tri_mask),  ORprob(tri_mask),'Type','Spearman' );
        ORProbToGP(a,1) = corr( gpTri(tri_mask),  ORprob(tri_mask),'Type','Spearman' );

        GPToPhono(a,1) = corr(phonemeTri(tri_mask), gpTri(tri_mask),'Type','Spearman');
        GPToLetter(a,1) = corr(letterTri(tri_mask), gpTri(tri_mask),'Type','Spearman' );
        GPToSem(a,1) = corr(semanticTri(tri_mask), gpTri(tri_mask),'Type','Spearman' );
        GPToGraph(a,1) = corr(graphemeTri(tri_mask), gpTri(tri_mask),'Type','Spearman' );
        GPToWordFreq(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToOrthoNeigh(a,1) = corr( OLD20sum(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToPhonoNeigh(a,1) = corr( PLD20sum(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToPhonoGraphNeigh(a,1) = corr( Phonographic_N(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToSemNeighDen(a,1) = corr( Sem_N_D(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToBigramFreq(a,1) = corr( BigramF_Avg_C_Log(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToNLetters(a,1) = corr( Length(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToNSyll(a,1) = corr( NSyll(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToBiphonProb(a,1) = corr( Biphon_Prob(tri_mask),  gpTri(tri_mask),'Type','Spearman' );
        GPToOR(a,1) = corr( ORTri(tri_mask),  gpTri(tri_mask),'Type','Spearman' );

        ORToPhono(a,1) = corr(phonemeTri(tri_mask), ORTri(tri_mask),'Type','Spearman' );
        ORToLetter(a,1) = corr(letterTri(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToSem(a,1) = corr(semanticTri(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToGraph(a,1) = corr(graphemeTri(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToWordFreq(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToOrthoNeigh(a,1) = corr( OLD20sum(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToPhonoNeigh(a,1) = corr( PLD20sum(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToPhonoGrpahNeigh(a,1) = corr( Phonographic_N(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToSemNeighDen(a,1) = corr( Sem_N_D(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToBigramFreq(a,1) = corr( BigramF_Avg_C_Log(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToNLetters(a,1) = corr( Length(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToNSyll(a,1) = corr( NSyll(tri_mask),  ORTri(tri_mask),'Type','Spearman' );
        ORToBiphonProb(a,1) = corr( Biphon_Prob(tri_mask),  ORTri(tri_mask),'Type','Spearman' );

        WordFreqToPhono(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  phonemeTri(tri_mask),'Type','Spearman' );
        WordFreqToLetter(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  letterTri(tri_mask),'Type','Spearman' );
        WordFreqToSem(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );
        WordFreqToGraph(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  graphemeTri(tri_mask),'Type','Spearman' );
        WordFreqToOrthoNeigh(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  OLD20sum(tri_mask),'Type','Spearman' );
        WordFreqToPhonoNeigh(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  PLD20sum(tri_mask),'Type','Spearman' );
        WordFreqToPhonoGraphNegih(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  Phonographic_N(tri_mask),'Type','Spearman' );
        WordFreqToSemNeighDen(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  Sem_N_D(tri_mask),'Type','Spearman' );
        WordFreqToBigramFreq(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  BigramF_Avg_C_Log(tri_mask),'Type','Spearman' );
        WordFreqToNLetters(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  Length(tri_mask),'Type','Spearman' );
        WordFreqToNSyll(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  NSyll(tri_mask),'Type','Spearman' );
        WordFreqToBiphonProb(a,1) = corr( Freq_SUBTLEXUS(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );

        OrthoNeighToPhonoNeigh(a,1) = corr( OLD20sum(tri_mask),  PLD20sum(tri_mask), 'Type','Spearman' );
        OrthoNeighToPhonoGraphNeigh(a,1) = corr( OLD20sum(tri_mask),  Phonographic_N(tri_mask), 'Type','Spearman' );
        OrthoNeighToSemNeighDen(a,1) = corr( OLD20sum(tri_mask),  Sem_N_D(tri_mask), 'Type','Spearman' );
        OrthoNeighToBigramFreq(a,1) = corr( OLD20sum(tri_mask),  BigramF_Avg_C_Log(tri_mask), 'Type','Spearman' );
        OrthoNeighToNLetters(a,1) = corr( OLD20sum(tri_mask),  Length(tri_mask), 'Type','Spearman' );
        OrthoNeighToNSyll(a,1) = corr( OLD20sum(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
        OrthoNeighToBiphonProb(a,1) = corr( OLD20sum(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
        OrthoNeighToPhono(a,1)= corr( OLD20sum(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
        OrthoNeighToLetter(a,1)= corr( OLD20sum(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
        OrthoNeighToGraph(a,1)= corr( OLD20sum(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );
        OrthoNeighToSem(a,1)= corr( OLD20sum(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );

        PhonoNeighToPhonoGraphNeigh(a,1) = corr( PLD20sum(tri_mask),  Phonographic_N(tri_mask), 'Type','Spearman' );
        PhonoNeighToSemNeighDen(a,1) = corr( PLD20sum(tri_mask),  Sem_N_D(tri_mask), 'Type','Spearman' );
        PhonoNeighToBigramFreq(a,1) = corr( PLD20sum(tri_mask),  BigramF_Avg_C_Log(tri_mask), 'Type','Spearman' );
        PhonoNeighToNLetters(a,1) = corr( PLD20sum(tri_mask),  Length(tri_mask), 'Type','Spearman' );
        PhonoNeighToNSyll(a,1) = corr( PLD20sum(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
        PhonoNeighToBiphonProb(a,1) = corr( PLD20sum(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
        PhonoNeighToSem(a,1)= corr( PLD20sum(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );
        PhonoNeighToPhono(a,1)= corr( PLD20sum(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
        PhonoNeighToLetter(a,1)= corr( PLD20sum(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
        PhonoNeighToGraph(a,1)= corr( PLD20sum(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );

        PhonoGraphNeighToSemNeighDen(a,1) = corr( Phonographic_N(tri_mask),  Sem_N_D(tri_mask), 'Type','Spearman' );
        PhonoGraphNeighToBigramFreq(a,1) = corr( Phonographic_N(tri_mask),  BigramF_Avg_C_Log(tri_mask), 'Type','Spearman' );
        PhonoGraphNeighToNLetters(a,1) = corr( Phonographic_N(tri_mask),  Length(tri_mask), 'Type','Spearman' );
        PhonoGraphNeighToNSyll(a,1) = corr( Phonographic_N(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
        PhonoGraphNeighToBiphonProb(a,1) = corr( Phonographic_N(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
        PhonoGraphNeighToSem(a,1)= corr( Phonographic_N(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );
        PhonoGraphNeighToPhono(a,1)= corr( Phonographic_N(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
        PhonoGraphNeighToLetter(a,1)= corr( Phonographic_N(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
        PhonoGraphNeighToGraph(a,1)= corr( Phonographic_N(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );

        SemNeighDenToBigramFreq(a,1) = corr( Sem_N_D(tri_mask),  BigramF_Avg_C_Log(tri_mask), 'Type','Spearman' );
        SemNeighDenToNLetters(a,1) = corr( Sem_N_D(tri_mask),  Length(tri_mask), 'Type','Spearman' );
        SemNeighDenToNSyll(a,1) = corr( Sem_N_D(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
        SemNeighDenToBiphonProb(a,1) = corr( Sem_N_D(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
        SemNeighDenToSem(a,1)= corr( Sem_N_D(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );
        SemNeighDenToPhono(a,1)= corr( Sem_N_D(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
        SemNeighDenToLetter(a,1)= corr( Sem_N_D(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
        SemNeighDenToGraph(a,1)= corr( Sem_N_D(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );

        BigramFreqToNLetters(a,1) = corr( BigramF_Avg_C_Log(tri_mask),  Length(tri_mask), 'Type','Spearman' );
        BigramFreqToNSyll(a,1) = corr( BigramF_Avg_C_Log(tri_mask),  NSyll(tri_mask), 'Type','Spearman' );
        BigramFreqToBiphonProb(a,1) = corr( BigramF_Avg_C_Log(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
        BigramFreqToSem(a,1)= corr( BigramF_Avg_C_Log(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );
        BigramFreqToPhono(a,1)= corr( BigramF_Avg_C_Log(tri_mask),  phonemeTri(tri_mask), 'Type','Spearman' );
        BigramFreqToLetter(a,1)= corr( BigramF_Avg_C_Log(tri_mask),  letterTri(tri_mask), 'Type','Spearman' );
        BigramFreqToGraph(a,1)= corr( BigramF_Avg_C_Log(tri_mask),  graphemeTri(tri_mask), 'Type','Spearman' );

        GraphToLetter(a,1) = corr( graphemeTri(tri_mask),  letterTri(tri_mask),'Type','Spearman' );
        GraphToPhono(a,1) = corr( graphemeTri(tri_mask),  phonemeTri(tri_mask),'Type','Spearman' );
        GraphToNLetters(a,1)= corr( graphemeTri(tri_mask),  Length(tri_mask),'Type','Spearman' );
        GraphToNSyll(a,1)= corr( graphemeTri(tri_mask),  NSyll(tri_mask),'Type','Spearman' );
        GraphToBiphonProb(a,1)= corr( graphemeTri(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );
        GraphToSem(a,1)= corr( graphemeTri(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );

        PhonoToLetter(a,1) = corr( phonemeTri(tri_mask),  letterTri(tri_mask),'Type','Spearman' );
        PhonoToNLetters(a,1)= corr( phonemeTri(tri_mask),  Length(tri_mask),'Type','Spearman' );
        PhonoToNSyll(a,1)= corr( phonemeTri(tri_mask),  NSyll(tri_mask),'Type','Spearman' );
        PhonoToBiphonProb(a,1)= corr( phonemeTri(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );
        PhonoToSem(a,1)= corr( phonemeTri(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );

        LetterToNLetters(a,1)= corr( letterTri(tri_mask),  Length(tri_mask),'Type','Spearman' );
        LetterToNSyll(a,1)= corr( letterTri(tri_mask),  NSyll(tri_mask),'Type','Spearman' );
        LetterToBiphonProb(a,1)= corr( letterTri(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );
        LetterToSem(a,1)= corr( letterTri(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );

        NLettersToBiphonProb(a,1)= corr( Length(tri_mask),  Biphon_Prob(tri_mask),'Type','Spearman' );
        NLettersToSem(a,1)= corr( Length(tri_mask),  semanticTri(tri_mask),'Type','Spearman' );

        NSyllToBiphonProb(a,1) = corr( NSyll(tri_mask),  Biphon_Prob(tri_mask), 'Type','Spearman' );
        NSyllToSem(a,1)= corr( NSyll(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );

        BiphonProbToSem(a,1) = corr( Biphon_Prob(tri_mask),  semanticTri(tri_mask), 'Type','Spearman' );

        NSyllToNLetters(a,1) = corr( NSyll(tri_mask),  Length(tri_mask), 'Type','Spearman' );

        %Example of getting 1 stdev for the whole correlation matrix
        w=0; %(w=0 is n-1, w=1 is n for stdev)
        std_phonemeTri(a,1) =std(phonemeTri(tri_mask),w,"all");
        std_gpTri(a,1) =std(gpTri(tri_mask),w,"all");
        std_letterTri(a,1) =std(letterTri(tri_mask),w,"all");
        std_semanticTri(a,1) =std(semanticTri(tri_mask),w,"all");
        std_graphemeTri(a,1) =std(graphemeTri(tri_mask),w,"all");
        std_ORTri(a,1) =std(ORTri(tri_mask),w,"all");
        std_Freq_SUBTLEXUS(a,1) =std(Freq_SUBTLEXUS(tri_mask),w,"all");
        std_OLD20sum(a,1) =std(OLD20sum(tri_mask),w,"all");
        std_PLD20sum(a,1) =std(PLD20sum(tri_mask),w,"all");
        std_Phonographic_N(a,1) =std(Phonographic_N(tri_mask),w,"all");
        std_Sem_N_D(a,1) =std(Sem_N_D(tri_mask),w,"all");
        std_BigramF_Avg_C_Log(a,1) =std(BigramF_Avg_C_Log(tri_mask),w,"all");
        std_Length(a,1) =std(Length(tri_mask),w,"all");
        std_NSyll(a,1) =std(NSyll(tri_mask),w,"all");
        std_Biphon_Prob(a,1) =std(Biphon_Prob(tri_mask),w,"all");
    end
    % print out results every x # of iterations
    % it=it+1;
    %  tally=tally+1;
    % if tally==increment
    %% writes everything in excel spreadsheets
    % iteration(1:it,1)=(1:it)';
    %wordindex(:,1) = iteration; %Number rows 1 through number of iterations for word index
    % set up table based on variables

    finalResults=table(iteration,NSyllToNLetters,GPToPhono	,GPToLetter	,GPToSem,GPToGraph,GPToWordFreq,GPToOrthoNeigh,GPToPhonoNeigh,GPToPhonoGraphNeigh,GPToSemNeighDen,  ...
        GPToBigramFreq,GPToNLetters,GPToNSyll,GPToBiphonProb,GPToOR,ORToPhono,ORToLetter,ORToSem,ORToGraph,ORToWordFreq,  ...
        ORToOrthoNeigh,ORToPhonoNeigh,ORToPhonoGrpahNeigh,ORToSemNeighDen,ORToBigramFreq,ORToNLetters,ORToNSyll,ORToBiphonProb,WordFreqToPhono,  ...
        WordFreqToLetter,WordFreqToSem,WordFreqToGraph,WordFreqToOrthoNeigh,WordFreqToPhonoNeigh,WordFreqToPhonoGraphNegih,  ...
        WordFreqToSemNeighDen,WordFreqToBigramFreq,WordFreqToNLetters,WordFreqToNSyll,  ...
        WordFreqToBiphonProb,OrthoNeighToPhonoNeigh,OrthoNeighToPhonoGraphNeigh,OrthoNeighToSemNeighDen,OrthoNeighToBigramFreq,OrthoNeighToNLetters,OrthoNeighToNSyll,  ...
        OrthoNeighToBiphonProb,OrthoNeighToPhono,OrthoNeighToLetter,OrthoNeighToGraph,OrthoNeighToSem,PhonoNeighToPhonoGraphNeigh,PhonoNeighToSemNeighDen,  ...
        PhonoNeighToBigramFreq,PhonoNeighToNLetters,PhonoNeighToNSyll,PhonoNeighToBiphonProb,PhonoNeighToSem,PhonoNeighToPhono,PhonoNeighToLetter,PhonoNeighToGraph,  ...
        PhonoGraphNeighToSemNeighDen,PhonoGraphNeighToBigramFreq,PhonoGraphNeighToNLetters,PhonoGraphNeighToNSyll,PhonoGraphNeighToBiphonProb,  ...
        PhonoGraphNeighToSem,PhonoGraphNeighToPhono,PhonoGraphNeighToLetter,PhonoGraphNeighToGraph,SemNeighDenToBigramFreq,SemNeighDenToNLetters,SemNeighDenToNSyll,  ...
        SemNeighDenToBiphonProb,SemNeighDenToSem,SemNeighDenToPhono,SemNeighDenToLetter,SemNeighDenToGraph,BigramFreqToNLetters,  ...
        BigramFreqToNSyll,BigramFreqToBiphonProb,BigramFreqToSem,  ...
        BigramFreqToPhono,BigramFreqToLetter,BigramFreqToGraph,GraphToLetter,GraphToPhono,GraphToNLetters,GraphToNSyll,  ...
        GraphToBiphonProb,GraphToSem,PhonoToLetter,PhonoToNLetters,PhonoToNSyll,PhonoToBiphonProb,PhonoToSem,LetterToNLetters,  ...
        LetterToNSyll,LetterToBiphonProb,LetterToSem,NLettersToBiphonProb,NLettersToSem,NSyllToBiphonProb,NSyllToSem,BiphonProbToSem, ...
        GPProbToPhono,GPProbToLetter,GPProbToSem,GPProbToGraph,GPProbToWordFreq	,GPProbToOrthoNeigh,GPProbToPhonoNeigh, ...
        GPProbToPhonoGraphNeigh	,GPProbToSemNeighDen	,GPProbToBigramFreq	,GPProbToNLetters	,GPprobToNSyll	,GPProbToBiphonProb, ...
        GPProbToOR	,GPProbToGP	,GPProbToORProb	,ORProbToPhono	,ORProbToLetter	,ORProbToSem	,ORProbToGraph	,ORProbToWordFreq	, ...
        ORProbToOrthoNeigh	,ORProbToPhonoNeigh	,ORProbToPhonoGraphNeigh	,ORProbToSemNeighDen	,ORprobToBigramFreq	, ...
        ORProbToNLetters	,ORProbToNSyll	,ORProbToBiphonProb	,ORProbToOR	,ORProbToGP, ...
        std_phonemeTri,std_gpTri,std_letterTri, std_semanticTri,std_graphemeTri,std_ORTri,std_Freq_SUBTLEXUS, std_OLD20sum, ...
        std_PLD20sum,std_Phonographic_N,std_Sem_N_D, std_BigramF_Avg_C_Log,std_Length,std_NSyll, std_Biphon_Prob);

    %% output for SCRIPT uncensor if running as SCRIPT
    % writematrix(wordindex, fullfile('output',strcat(output_code,sprintf('_Orthogonal_WordIndex_iterations_%d',it),'.csv')), 'WriteMode', 'Append'); %write word index
    % writetable(finalResults, fullfile('output',strcat(output_code,sprintf('_Orthogonal_Corrs_iterations_%d',it),'.csv')),'WriteMode', 'Append'); %write orthog correlations
    %% output for FUNCTION (censor if running as SCRIPT)
    writematrix(wordindex, fullfile(output_path,'output',strcat(output_code,sprintf('_Orthogonal_WordIndex_iterations%d_group%d',numOfIterations,pa),'.csv')), 'WriteMode', 'Append'); %write word index
    writetable(finalResults, fullfile(output_path,'output',strcat(output_code,sprintf('_Orthogonal_Corrs_iterations%d_group%d',numOfIterations,pa),'.csv')),'WriteMode', 'Append'); %write orthog correlations
    
    %   tally=0;
    % end
    %end

end
disp(datetime)
