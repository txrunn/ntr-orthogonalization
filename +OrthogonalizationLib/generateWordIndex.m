function output = generateWordIndex(wordlist)
    % Function to produce the word index of the final list.

    % Initialize the word index
    wordindex = cell(size(wordlist, 1), size(wordlist, 2) + 1);

    % Populate the word index with words from a list
    wordindex(:, 2:end) = transpose(wordlist);

    % Assign iteration numbers to the first column of the word index
    wordindex(:, 1) = arrayfun(@(x) sprintf('%d', x), 1:size(wordlist, 1), 'UniformOutput', false);

    output = wordindex; % Return the populated word index
end
