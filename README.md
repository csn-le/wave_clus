Wave_clus
=========

Introduction
------------

Wave_clus is a fast and unsupervised algorithm for spike detection and sorting using wavelets and super-paramagnetic clustering. Although it gives a first unsupervised solution, this can be further modified according to
the experimenters’ preference (semi-automatic sorting). Wave_clus is free (and therefore without any warranty) for any non commercial applications. For any commercial application please contact Rodrigo Quian Quiroga.

#### Requirements
Wave_clus runs under Windows, Linux and Mac. It requires MATLAB 7.6 (R2008a) or higher. It uses the functions clusterXX.exe, provided by Eytan Domani, which is an executable that does the superparamagnetic clustering of the data. The wavelet and the signal processing toolboxes are not necessary.

How-to
------

#### Installation
In MATLAB go to the menu File/Set Path and add the directory wave\_clus with subfolders to the MATLAB path. The functions from [FilterM](https://www.mathworks.com/matlabcentral/fileexchange/32261-filterm "FilterM") can be added to the path for using them instead of the functions of the Signal Processing Toolbox.


#### Basic Gui Instructions
1. Open the GUI of Wave_clus, typing `wave_clus` at the MATLAB command prompt. 
2. (Optional) Edit the parameters with the *Set\_parameter* button.
3. Press *Load* and select the raw data. Read the *Input Files* section in this document.
4. Explore the clustering output at different temperatures (you will select a minimun cluster size with the same click). Remember to fix clear classes before move to another temperature.
* (Optional) Select spikes manually to create new classes.
* (Optional) Merge split classes (select classes with *fix* button and press *Merge*).
5. (Optional) Reject spurious classes.
6. Force the rest of the spikes.
* Remove non-spike events from classes, selecting them manually and rejecting the new class.
7. Press *Save* to save the current result.

#### Batch_files
The automatic detection and clustering can be performed by the functions Get\_spikes and Do\_clustering respectively (see the .m files in batch_files/ for more instructions).


#### Output Files


The output of Wave_clus (obtained either using the *Save* clusters button in the GUI or the one given automatically by the batch files) is times_[filename].mat, which is a MATLAB file containing the following variables: **par** (parameters used for clustering), **spikes** (a matrix with the spike shapes), **inspk** (a matrix with the features of the spike shapes) and **cluster_class** (a matrix with the clustering results). The variable **cluster_class** has 2 columns and nspk rows (nspk is the number of spikes). The first column is the cluster class, with integers denoting the clusters membership and a value of 0 for those spikes not assigned to any cluster. The second column is the spike times in ms.

#### Input Files

Wave_clus can read MATLAB files (extension .mat) with continuous data or spikes for clustering spike shapes that have already been detected (e.g. detected on-line by the acquisition system). It should have either a vector named **data** (the continuous signal) or a matrix named **spikes** (nr. of spikes x length of the spike shape) plus a vector **index** with the spike times. If the variable **sr** is inside the file, it will set the sampling rate. Otherwise **par.sr** inside the file `set_parameters` will be use.

All the supported formats (.mat, .int, .NSx, .pl2, .tdt and .ncs) use the codes in the folder `Raw_data_readers` to get the data from the files. Some of them require to run the codes in the folder `tools` before.


Important links
---------------

Questions, problems and suggestions: [issues section](https://github.com/csn-le/wave_clus/issues "Issues").

More instructions, FAQ and developer information: [wiki](https://github.com/csn-le/wave_clus/wiki "Wiki").

Documentation and sample data: [official website](http://www2.le.ac.uk/centres/csn/research-2/spike-sorting "Official Website").


References
----------

#### How to cite

__Unsupervised spike detection and sorting with wavelets and superparamagnetic clustering.__
R. Quian Quiroga, Z. Nadasdy and Y. Ben-Shaul
Neural Computation 16, 1661-1687; 2004.


######For a non technical reference about spike sorting see:

Quick guide: [Spike Sorting](http://www2.le.ac.uk/departments/engineering/research/bioengineering/neuroengineering-lab/Publications/spike%20sorting%20quick%20guide.pdf "quick guide")<br/>
Quian Quiroga, R.<br/>
Current Biology, Vol 22. R45–R46, 2012.

[Spike Sorting](https://www.scholarpedia.org/article/Spike_sorting "spike sorting in Scholarpedia")<br/>
R. Quian Quiroga<br/>
Scholarpedia 2 (12): 3583. 2007


