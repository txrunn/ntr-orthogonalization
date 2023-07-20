% Utilizing the OrthogonilizationLib
% ... any setup needed ...
wordListInputFile = 'wordinput_1.csv'; % or other source of word list

% Create NTR_Orthogonalization object
ortho = NTR_Orthogonalization(wordListInputFile);

% perform filtering
ortho.filterData();

% perform orthogonalization
[allMatrices, allWordParameters] = OrthogonalizationLib.NTR_Orthogonalization.performOrthogonalization(ortho.parameters);

% generate correlation data
correlationData = OrthogonalizationLib.generateCorrelationData.performGeneration(allMatrices);

% generate word index
wordIndex = OrthogonalizationLib.generateWordIndex.performGeneration(ortho.parameters);

% visualize and save distance matrices
OrthogonalizationLib.visualizeDistanceMatrix.visualizeAndSave(allMatrices, ortho.parameters.string);

% ... any cleanup ...
