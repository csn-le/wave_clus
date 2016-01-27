Wave_clus
=========

Introduction
------------

Wave_clus is a fast and unsupervised algorithm for spike detection and sorting using wavelets and super-paramagnetic clustering. Although it gives a first unsupervised solution, this can be further modified according to
the experimenters’ preference (semi-automatic sorting). Wave_clus is free (and therefore without any warranty) for any non commercial applications. For any commercial application please contact Rodrigo Quian Quiroga.

#### Requirements
Wave_clus runs under Windows, Linux and Mac. It requires Matlab 7.6 (R2008a) or higher. It uses the functions clusterXX.exe, provided by Eytan Domani, which is an executable that does the superparamagnetic clustering of the data. The wavelet and the signal processing toolboxes are not necessary.

How-to
------

#### Installation
In matlab go to the menu File/Set Path and add the directory wave\_clus with subfolders to the matlab path. The functions from [FilterM](http://www.mathworks.com/matlabcentral/fileexchange/32261-filterm "FilterM") can be added to the path for using them instead of the functions of the Signal Processing Toolbox.


#### Usage

###### Basic Gui Instructions
1. Open the GUI of Wave_clus, typing `wave\_clus` at the MATLAB command prompt. 
2. (Optional) Edit the parameters with the *Set\_parameter* button.
3. Load the raw data. Simplest case: a .mat file with the single channel recording in a variable *data* and the sample rate in a variable *sr*.
4. Explore the temperature diagram fixing clear classes.
⋅⋅* (Optional) Select spikes manually to create new classes.
⋅⋅* (Optional) Merge splited classes (select classes with *fix* button and push *Merge*).
5. (Optional) Reject spurious classes.
6. Force the rest of the spikes.
⋅⋅* Remove non-spike events from classes, selecting them manually and rejecting the new class.
7. Save the results.

The automatic detection and clustering can be performed by the functions Get\_spikes and Do\_clustering respectively(see the .m files in batch_files/ for more instructions).

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


