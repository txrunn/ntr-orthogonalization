
# Orthogonalization Library

This library contains the core functionalities for the NTR Orthogonalization process. The methods within this library perform a series of operations on wordlists to generate distance matrices, visualize them, and run orthogonalization procedures.

## Classes and Their Descriptions

- `NTR_Orthogonalization`: This is the primary class that orchestrates the orthogonalization process. It calls various methods to handle different segments of the process.
- `generateWordIndex`: Handles the extraction and generation of word indices from the provided datasets.
- `generateDistanceMatrices`: Manages the creation of distance matrices based on word correlations.
- `generateCorrelationData`: Extracts correlation data from the datasets.
- `visualizeDistanceMatrix`: Provides visualization utilities to generate and save heatmaps of distance matrices.
