# NTR Orthogonalization

This repository includes code and algorithm originally developed by Audrey Lyu for orthogonalization of word lists. It generates the distance matrices, correlational data, and parameters of each word from the final list(s) of worths that are minimally correlated. The current version allows input from a single or multiple lists. Our team has adapted and enhanced this original code for better performance and scalability.

 Library

This library provides functionalities related to the NTR Orthogonalization process.

## Directory Structure

- `README.md`: This file.
- `distanceMatFigures.m`: A script related to visualizing distance matrices (needs further integration into the library).
- `+OrthogonalizationLib`: The main library directory containing core functions and classes.
- `data`: Directory containing necessary data files.
- `results`: Directory containing result or output files.
- `archive`: Directory containing older versions or variations of scripts.

## Usage

To use the library, add the `+OrthogonalizationLib` directory to your MATLAB path and call the necessary functions or classes.

## Further Development

- Integrate the functionality of `distanceMatFigures.m` into the main library or utilize it as a utility script.
- Continuously update the library based on new requirements or functionalities.


### MATLAB
To generate figures and edit on MATLAB individually (change letterTri to matrix of interest):

```MATLAB
letterTriTrans = transpose(letterTri);
heatmap(letterTriTrans, 'Colormap', parula)
set(gca,'XData',wordlist, 'YData', wordlist)
```

### Figure Editing
To edit the figure:
1. Click “Property Inspector” on top of Figure menu in MATLAB ![image](https://github.com/txrunn/ntr-orthogonalization/assets/31973391/53ba9a67-9d1b-42f0-b568-b31bd46978b5)

2. Edit font, title, colors, etc. here
3. On this figure menu, you can also zoom into specific regions and export it as a separate file

## Goal
The overall goal is to combine these two algorithms and make it so that NTR_Orthogonalization_v15_multilistinput.m can output a file of all the visualized distance matrix color figures, and ensure that algorithms run with any list of words.

## Credits
The original code and algorithm were developed by Audrey Lyu. This enhanced version was adapted by our team.



### Updates on formatTables(obj) Method
The `formatTables(obj)` method in `NTR_Orthogonalization.m` has been updated to include logic for reading various datasets into tables:
1. Reading a masterlist into a table named `parameters` from `ntr_masterlist_gp.xlsx`.
2. Reading data into a table named `jpglove` from `jpglove.csv`.
3. Reading onset-rime input into a table named `ORTable` from `ntr_masterlist_onset_rimes.xlsx`.

This method ensures the datasets are read into tables and are available for further processing within the library.
