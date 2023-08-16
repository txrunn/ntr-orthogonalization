function distanceMatrix = generateDistanceMatrices(wordlist)
    % This function generates distance matrices based on the provided wordlist.

    numOfWords = length(wordlist);

    % Initialize the word index matrix
    wordindex = zeros(numOfWords, numOfWords + 1);

    for a = 1:numOfWords
        wordindex(a, 2:[numOfWords + 1]) = transpose(wordlist(1:numOfWords));
        wordindex(a, 1) = sprintf('%d', a); % iteration
    end

    % Create triangles for pairwise distances
    w = pdist(wordlist); % Determines pairwise distance
    E = squareform(w); % Formats distance matrix
    distanceMatrix = triu(E, -1); % Grabs upper triangular part of matrix from one diagonal below the center
end
