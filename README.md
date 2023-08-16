# NTR Orthogonalization Library

## Overview

The NTR Orthogonalization Library is a specialized MATLAB toolkit designed to aid fMRI block-based studies in constructing word lists for their experiments. By ensuring minimal correlation between words, the library helps scientists optimize the construction of word lists, thereby enhancing the precision of fMRI studies.

## Features

- **Minimized Word Correlation**: Construct word lists with minimal inter-word correlations.
- **Adjustable Word List Length**: Tailor the length of word lists while maintaining minimal correlation.
- **User-Friendly Interface**: Designed for scientists with varying levels of coding experience.
- **Customizable**: Accommodates unique datasets and parameters to fit specific project needs.
- **Visualization Tools**: Visualize distance matrices and correlations to aid study design.
- **Extensibility**: Designed for ease of extension to cater to related projects.

## Algorithm

The core algorithm behind the NTR Orthogonalization Library focuses on minimizing the correlation between words in block-based fMRI studies. The algorithm achieves this by:

1. Calculating pairwise correlations between words based on their semantic vectors.
2. Creating a triangular distance matrix that represents the semantic distances between words.
3. Iteratively adjusting the word list to minimize the overall correlation while maintaining the desired length.

This approach ensures that the word lists generated are both short and minimally correlated, optimizing the precision of fMRI studies.

## Getting Started

### Prerequisites

- MATLAB (R2019a or later recommended).
- Data files (e.g., `wordinput_1.csv`) should be placed in the `data` directory.

### Setup

1. Clone or download the repository to your local machine.
2. Open MATLAB and set the working directory to the root folder of this project.
3. Run the `main.m` script to access the primary functionality.
4. Explore the `+OrthogonalizationLib` directory for detailed implementations of the library's functions.

## Usage

The primary script to run the library's functionality is `main.m`. Executing this script will guide you through the process of creating an optimized word list for your fMRI block-based study.

For users looking to customize or extend the library's functionalities, the `+OrthogonalizationLib` directory contains the core MATLAB classes and methods. Each file within this directory is well-documented to provide clarity on its purpose, usage, inputs, and outputs.

## Customization and Extensions

The library is designed with extensibility in mind. If you have specific requirements or datasets, the library's modular structure allows for easy adaptations and extensions. For guidance on customizing the library, refer to -.

## Contribution and Support

Contributions to enhance the library are welcome. If you encounter any issues or have questions regarding its use, please raise an issue on the GitHub repository.

## License

This project is open-source. Please refer to the `LICENSE` file for more details.

## Acknowledgments

- The initial foundation of this library is based on the work done in the original fMRI block-based study script by Audrey Lyu.
