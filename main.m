% Utilizing the OrthogonalizationLib
% TODO: Add any setup or initializations needed
addpath('./+OrthogonalizationLib/');
% Specify the source of word list
wordListInputFile = '/data/wordinput_1.csv';

% Create NTR_Orthogonalization object
ortho = NTR_Orthogonalization(wordListInputFile);

% Perform data filtering
ortho.filterData();

% Perform orthogonalization logic
ortho.runIterations();
ortho.performOrthogonalization();
ortho.postProcessResults();
% TODO: Update any additional methods or steps needed for orthogonalization

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
