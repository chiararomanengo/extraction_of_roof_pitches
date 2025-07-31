# extraction_of_roof_piches

This repository contains a Matlab library for extracting roof pitches form a point cloud that represent a single building and characterizing them via a planes fitting algorithm based on the Hough transform. 

## Content of the repository

The ```code``` directory contains our method that is described in Sections 5.2-5.3 of the paper "Geometry-aware estimation of photovoltaic energy from aerial LiDAR point clouds".

The ```vlp``` directory contains the pont clouds that represent buildings of an area of Genoa (3820) georeferenced with EPSG 7791. 

Finally, the ```data``` directory contains building footprints corresponding to the area 3820 in the form of standard ESRI Shapefiles.

All the data are publicly available in the public geoportal (https://mappe.comune.genova.it/MapStore2/#/viewer/1000003072), hosted by the Genova Municipality.


## How to use it
To use this method, you can simply run the ```main.m``` file in Matlab. 
The input are .las files containing the point cloud representing a building. 
The point cloud is then processed and the method find the best fitting primitive type and its geometric descriptors.  
The output is, for each building, a folder that contains:
- a set of .off files that represent polygonal meshes encoding the geometry of each roof pitch;
- a .txt file containing the tilt angles of each roof pitch (the tilt with respect to the horizontal and azimuthal planes).
