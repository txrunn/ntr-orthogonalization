
% NTR_Orthogonalization - Main function for the Orthogonalization Library
%
% This script serves as the primary entry point for the Orthogonalization 
% process, helping researchers in fMRI block-based studies to create word lists 
% for their experiments. The goal is to ensure the words are minimally correlated.
%
% The library provides an easy way for studies to be executed by other scientists 
% and can be customized for related projects.
%
% The entire process is modularized into distinct functions, each addressing a 
% specific step in the process. The user can run this script after setting up 
% the necessary input files and parameters.
%
% INPUTS:
%   - wordinput: A CSV file containing the words to be processed.
%   - wordoutput: Name of the output MAT file where results will be saved.
%
% OUTPUTS:
%   - MAT file: Contains various variables and matrices resulting from the 
%               orthogonalization process.
%
% DEPENDENCIES:
%   This script relies on multiple helper functions and libraries present in the 
%   '+OrthogonalizationLib' folder.
%
% USAGE:
%   Modify the 'wordinput' and 'wordoutput' variables as per your study requirements.
%   Ensure all dependencies are in the MATLAB path.
%   Run the script.

% --- INPUT SETUP ---

function NTR_Orthogonalization
% Path to the input CSV file containing words for the study
wordinput = 'data/wordinput_1.csv';

% Name of the output MAT file where results will be saved
wordoutput = 'results/wordoutput_1.mat';

% Ensure the library is in the path
addpath('+OrthogonalizationLib');

% --- MAIN ORTHOGONALIZATION PROCESS ---

% Load words from the input CSV
words = OrthogonalizationLib.loadWordsFromCSV(wordinput);

% Generate word index
wordIndex = OrthogonalizationLib.generateWordIndex(words);

% Generate correlation data
correlationData = OrthogonalizationLib.generateCorrelationData(words, wordIndex);

% Process the words and generate distance matrices
distanceMatrices = OrthogonalizationLib.generateDistanceMatrices(correlationData);

% Visualize the generated distance matrices (optional)
% OrthogonalizationLib.visualizeDistanceMatrix(distanceMatrices, words);

% Save the results to the output MAT file
save(wordoutput, 'words', 'wordIndex', 'correlationData', 'distanceMatrices');