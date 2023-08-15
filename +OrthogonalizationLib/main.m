
% Utilizing the OrthogonalizationLib
% TODO: Add any setup or initializations needed

% Specify the source of word list
wordListInputFile = 'wordinput_1.csv';

% Create NTR_Orthogonalization object
ortho = NTR_Orthogonalization(wordListInputFile);

% Perform data filtering
ortho.filterData();

% TODO: Implement and call the orthogonalization method for the NTR_Orthogonalization object

% Generate correlation data
% TODO: Update the function name and arguments if needed
correlationData = OrthogonalizationLib.generateCorrelationData.performGeneration();

% Generate word index
% TODO: Update the function name and arguments if needed
wordIndex = OrthogonalizationLib.generateWordIndex.performGeneration();

% Visualize and save distance matrices
% TODO: Update the function name and arguments if needed
OrthogonalizationLib.visualizeDistanceMatrix.visualizeAndSave();

% TODO: Add any cleanup or final steps needed
