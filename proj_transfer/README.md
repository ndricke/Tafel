# Fitting unusual electrokinetic data
The goal of this project is to find a model that fits this electrochemical data, which no current models do in a satisfactory manner. In particular, the relationships in these datasets suggest there are fractional order kinetics occuring for each of these systems in a manner that is not reproducable by the standard steady-state approaches to electrochemical kinetics that are normally used.

## Folder: data\_csv
This folder contains the experimental datasets for several different electrochemical catalysts in projects that I worked on for both Yogi and Mercea. Each dataset provides information about the relationship between any of current, voltage, pH, or O2 concentration. Note: all data also has a corresponding image in the "data\_images" directory.

### Gold HER data
Files:

Jackson\_HER\_data.csv --> data stripped from figure. The pH the data was taken at isn't well ordered

Jackson\_HERfull\_data.csv --> data with pH values added manually. Can be loaded as a pandas dataframe

Notes:

This dataset came from Megan Jackson in Yogi Surendranath's lab. I don't think they've published this yet. I've included both the figure and the group meeting presentation that she had given with that figure. The main jist of it is that they see a lower pH dependence than would be expected for a first-order proton transfer when at high pH.

### N-GCC data
Files:

NGCC\_O2vsRate\_data.csv --> how current depends on the O2 concentration

NGCC\_Tafel\_data.csv --> standard current-voltage relationship in a Tafel plot

NGCC\_pHvsOnset.csv --> how onset potential depends on pH

Notes:

This is data collected by Alex Murray while he was in Yogi Surendranath's lab, and is in the paper that I published with him. There is a slightly lower than 1st order O2 dependence, and a weak pH dependence. 

### Nibis data
Files:

Nibis\_Tafel\_pH8\_data.csv --> Tafel data collected at pH 8

Nibis\_Tafel\_pH13\_data.csv --> Tafel data collected at pH 13

Nibis\_pHvsOnset\_data.csv --> onset potential as a function of pH

Notes:

This is data collected by Elise Miner while she was in Mercea Dinca's lab, and is in the paper that I published with her. Like the GCC's, there is a slightly lower than 1st order O2 dependence, and a weak pH dependence.


## Folder: data\_images
Contains the images that correspond to the files in data\_csv. 

## Python Scripts:

convertJackson.py: convert data on gold HER electrokinetics that was scraped from paper

DisorderMech.py: defines classes for modeling electrokinetics with disorder present

ElecMech.py: defines general classes for modeling electrokinetics

fit\_goldHER.py: python scripts for fitting Megan Jackson's HER data. This could likely be improved with more systematic parameter sweeps, and more sophisticated modeling of disorder

fit\_ngccTafel.py: script for fitting NGCC tafel data only

fit\_ngccTafel\_pcet2.py: Fits NGCC tafel data with pcet2 mechanism without considering pH or O2 concentration

fit\_ngccTafel\_pcetet.py: Fits NGCC tafel data with pcetet mechanism without considering pH or O2 concentration. Results should be exactly the same as pcet2 mechanism when not including pH

fit\_ngccTafelpH\_pcet2.py: Simultaneously fit pH vs onset potential data and tafel data with same mechanism

ngcc\_TafelpH.py: script for fitting NGCC tafel, pH, and O2 data simultaneously. This is not complete.

SimTafel.py: a class for simulating and fitting electrokinetic data

A series of test cases for showcasing the effects of disorder:

test\_disorder\_JvsV.py

test\_disorder\_onset.py

test\_disorder\_pH.py


## Folder: Related\_Documents
RickeAcsCat2017\_NGCC.pdf: a pdf of my paper with Yogi on NGCCs

MinerAcsCat2017\_Nibis.pdf: a pdf of Elise's paper that I collaborated on for the MOF Nibis

Jackson\_HalfOrderSlides.pptx: powerpoint presentation Megan gave during group meeting
