# NeighborhoodCoordination
Code accompanying our manuscript "Coordinated cellular neighborhoods orchestrate antitumoral immunity at the colorectal cancer invasive front"

## Introduction
This repo consists of two parts: 1) Code used to generate the results for our publication “Coordinated cellular neighborhoods orchestrate antitumoral immunity at the colorectal cancer invasive front” and 2) Useful functions from our paper that we have generalized to make more user-friendly for the broader scientific community.  These functions can be found in the “Neighborhoods” directory.  

The following describes how to navigate the code we provide in the “Neighborhoods” directory. The code is implemented in python3 using jupyter notebooks.  Some experience with installing python packages will be required to get up and running.  Only a familiarity with python is required to follow and adapt these functions for your own use.

## Installation
We recommend using the free package manager Anaconda (https://www.anaconda.com) to help install the packages required to use these scripts.
•	Jupyter is most easily installed through Anaconda.  Otherwise, one can follow the instructions here  https://jupyter.org/
•	The used scripts require very few dependencies.  Packages required:
o	statsmodels
o	numpy
o	pandas
o	jupyter
o	seaborn
o	scikitlearn
o	tensorly
o	shapely


## Functions
1.	The single-cell file (with patient, TMA core, neighborhood10 and ClusterName annotations) is in the Mendeley data link provided in the manuscript. Any time you see read in of 'cells2_salil', this can be replaced with read in of the Mendeley data file.

2.	Neighborhood Identification:  This notebook walks a user through identifying Cellular Neighborhoods in high parameter imaging data as described in our publication.

3.	Voronoi Generation:  This notebook helps visualize the Cellular Neighborhoods on the tissue and allows the user to overlay a collection of cells over these neighborhoods to explore their spatial distribution.

4.	tensor_decomposition_cleaned_up:  This notebook describes how to perform tensor decomposition after each single cell has been allocated to a Cellular Neighborhood and Cell Type. This has been updated to start from the Mendeley upload.

5.	Neighborhood Mixing:  This notebook allows the user to describe the spatial contacts between two Cellular Neighborhoods of interest.  

6.	Cell Type Differential Enrichment:  This notebook allows a user to identify cell types that are differentially enriched within specific cellular neighborhoods across a set of conditions or clinical groups.  

7.  app_CRC_contacts.R:  A Shiny application for computing contacts between different cell types.

