
# OrthogonalizationLib

`OrthogonalizationLib` is a MATLAB library designed to create orthogonalized matrices comparing phoneme units, Graph (The Levenshetein distance across grapheme units)eme units, phoneme + Graph (The Levenshetein distance across grapheme units)eme (GP) units, Letter (The Levenshetein distance across letters) edit distance (Levenshtein), and Sem (The Levenshetein distance across semantic glove variables) distance. This library provides a structured, object-oriented approach to the orthogonalization process, encapsulating the functionality originally found in the `NTR_Orthogonalization_v15_multilistinput.m` script.

## Components:

### 1. NTR_Orthogonalization Class:

#### Properties:
- **numOfWords**: Number of words.
- **numOfIterations**: Number of iterations.
- **increment**: Increment value.
- **output_path**: Path for outputs.
- **parameters**: Parameters for the orthogonalization process.
- **jpglove**, **ORTable**, **scope**, **ELP**, **biphone**: Tables loaded from various data sources.

#### Methods:
- **Constructor (`NTR_Orthogonalization`)**: Initializes the class, loads word lists, and sets up data.
- **filterData**: Performs data filtering operations by calling various filtering methods.
- **filterHomoGraph (The Levenshetein distance across grapheme units)s**: Filters out homoGraph (The Levenshetein distance across grapheme units)s from the parameters table.

(Note: The class contains several other methods that handle different filtering and processing steps.)

### 2. Helper Functions:
- **visualizeDistanceMatrix**: Provides visual representations of distance matrices.
- **generateDistanceMatrices**: Produces distance matrices.
- **generateCorrelationData**: Produces correlation data.
- **generateWordIndex**: Produces the word index of the final list.

## Workflow:

1. **Setup**: Specify the word list input file.
2. **Data Loading**: Load word lists and other necessary data.
3. **Data Filtering**: Filter out homoGraph (The Levenshetein distance across grapheme units)s, words based on length, and other criteria.
4. **Orthogonalization**: Execute the orthogonalization process.
5. **Data Generation**: Produce correlation data and word indices.
6. **Visualization**: Create and save visual representations of distance matrices.

## Usage:

Run the `main.m` script to execute the workflow using the functions and methods provided in the `OrthogonalizationLib`.

---

This README provides a basic overview of the `OrthogonalizationLib`. For detailed usage and further information, refer to the comments within each function and method.


### filterHomographs method
This method filters out homographs from the parameters. The number of filtered words is displayed.

### filterWordLength method
This method filters out words based on length criteria. The number of filtered words is displayed.

### filterGLOVEParameters method (Placeholder)
This method is intended for filtering words based on GLOVE parameters. The specific logic is yet to be implemented.

### filterELPData method (Placeholder)
This method is intended for filtering and operations based on ELP data. The specific logic is yet to be implemented.

### filterScopeData method (Placeholder)
This method is intended for filtering and operations based on scope data. The specific logic is yet to be implemented.

### filterBiphoneProbabilityData method (Placeholder)
This method is intended for filtering and operations based on biphone probability data. The specific logic is yet to be implemented.

### manipulateStrings method
This method is intended for string manipulations such as trimming and replacing characters. The specific manipulations are to be implemented based on context.

### matchInputsToParameters method (Placeholder)
This method is intended for matching other inputs to parameters. The specific logic is yet to be implemented.

### formatTables method (Placeholder)
This method is intended for formatting tables related to gp, bigp, and onset-rimes. The specific logic is yet to be implemented.

### runIterations method (Placeholder)
This method is intended for running iterations and computations. The specific logic is yet to be implemented.

### processNextSegment method (Placeholder)
This method is intended for operations related to the next segment of the script. The specific logic is yet to be implemented.

### processFurtherSegment method (Placeholder)
This method is intended for operations related to the further segment of the script. The specific logic is yet to be implemented.

### processSubsequentSegment method (Placeholder)
This method is intended for operations related to the subsequent segment of the script. The specific logic is yet to be implemented.

### processNextSegment2 method (Placeholder)
This method is intended for operations related to the next segment of the script (2). The specific logic is yet to be implemented.

### processUpcomingSegment method (Placeholder)
This method is intended for operations related to the upcoming segment of the script. The specific logic is yet to be implemented.

### processSubsequentSegment2 method (Placeholder)
This method is intended for operations related to the subsequent segment of the script (2). The specific logic is yet to be implemented.

### processFollowingSegment method (Placeholder)
This method is intended for operations related to the following segment of the script. The specific logic is yet to be implemented.

### processNextSegment3 method (Placeholder)
This method is intended for operations related to the next segment of the script (3). The specific logic is yet to be implemented.

### processSubsequentSegment3 method (Placeholder)
This method is intended for operations related to the subsequent segment of the script (3). The specific logic is yet to be implemented.

### processNextSegment4 method (Placeholder)
This method is intended for operations related to the next segment of the script (4). The specific logic is yet to be implemented.

### processSubsequentSegment4 method (Placeholder)
This method is intended for operations related to the subsequent segment of the script (4). The specific logic is yet to be implemented.

### processFollowingSegment2 method (Placeholder)
This method is intended for operations related to the following segment of the script (2). The specific logic is yet to be implemented.

### processSubsequentSegment5 method (Placeholder)
This method is intended for operations related to the subsequent segment of the script (5). The specific logic is yet to be implemented.