ToothFinder
===========

Analyses images of leaves and finds marginal teeth.

As part of the Morphidas project (http://www.computing.surrey.ac.uk/morphidas), we developed software that finds marginal teeth given the image of a single leaf. It traces the outline of the leaf and locates local maxima and minima to find the teeth. It then calculates the area of each tooth, the tip angle and several related features. A paper describing the work was published in PLoS One (http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0042112)

This software requires Matlab and the Mathworks Image Processing Toolbox; it has been developed and tested using Matlab version 7.10.0 (R2010a). 

An example leaf image, in Matlab's .mat format, can be downloaded here as a demonstration. To run the code using this example leaf, enter the following command:

    [totalArea, angle, frequency, toothCount, trilength, triratio, innerLength, toothList] = toothFinder('edge_1164_1.mat',1); 

 
