
# TODO for NTR Orthogonalization Library Integration

## Completed Sections

- **Initialization & Parameters**: Extracted parameters, read word lists, loaded paths, tables, and parameters.
- **Data Filtering**: Integrated methods for filtering out homographs, words by length, and other criteria. Filters applied to various tables like parameters, GLOVE, ELP, etc.
- **Triangular Matrix Initialization**: Identified and integrated sections from the original script related to this.
- **Distance Matrix Logic**: Identified and integrated logic related to the creation of triangular matrices for pairwise distances.

## Pending Sections

- **Generalization**: Make the program adaptable to use any word list.
- **Parallelization Optimization**: Ensure the integrated logic is optimized for parallel processing.
- **Heatmap Visualization**: Logic for visualizing distance matrices as heatmaps. Identification and integration pending.
- **Output & Saving**: Final steps of the algorithm where results are saved and outputs are generated. Yet to be integrated.
- **Review `main.m`**: Ensure it functions similar to the original v15 script.

## Integration Challenges

- **Data Structures**: The original script uses certain MATLAB-specific data structures. Need to ensure compatibility and seamless integration.
- **Logic Distribution**: The original script is a sequential set of instructions. The class-based approach in the library might require redistribution and organization of logic.

## Additional Notes

- Continuous review and validation are essential. After integrating each section, we need to validate the output against the original script to ensure correctness.
