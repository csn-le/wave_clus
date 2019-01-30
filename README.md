Wave_clus 3
===========

Wave_clus is a fast and unsupervised algorithm for spike detection and sorting that
runs under Windows, Mac or Linux operating systems.

To install, download this repository into a folder. In MATLAB (R2009b or higher) go to *Set Path* and add the directory wave_clus with subfolders to the MATLAB path.

How to use
------
To open the GUI type `wave_clus` in the MATLAB command prompt. To use the batch functions type `Get_spikes(‘filename’)` to do the spike detection
and `Do_clustering(‘filename)` to do the sorting.
Wave_clus generates a file *times_filename.mat*, with a variable **cluster_class** of two columns: the first column
is the class of the spike and the second one is the spike time in ms.

Wave_clus can read raw data or already detected spikes generated with electrophysiology data acquisition systems (Blackrock, Neuralynx, Plexon, TDT,
Intan, etc) or saved as matlab files (.mat). It can also deal with tetrodes and high-
density probes.

Wave_clus is free (and therefore without any warranty) for any non-commercial applications. For any commercial application please contact Prof. Rodrigo Quian Quiroga.

Important links
---------------

Questions, problems and suggestions: [issues section](https://github.com/csn-le/wave_clus/issues "Issues").

More instructions, details, FAQ and developer information: [wiki](https://github.com/csn-le/wave_clus/wiki "Wiki").



References
----------

#### How to cite
__A novel and fully automatic spike sorting implementation with variable number of features.__
F. J. Chaure, H. G. Rey and R. Quian Quiroga. Journal of Neurophysiology; 2018. 
[https://doi.org/10.1152/jn.00339.2018](https://www.physiology.org/doi/full/10.1152/jn.00339.2018)

###### For a non technical reference about spike sorting see:

Quick guide: [Spike Sorting](http://www2.le.ac.uk/departments/engineering/research/bioengineering/neuroengineering-lab/Publications/spike%20sorting%20quick%20guide.pdf "quick guide")<br/>
Quian Quiroga, R.<br/>
Current Biology, Vol 22. R45–R46, 2012.

[Spike Sorting](https://www.scholarpedia.org/article/Spike_sorting "spike sorting in Scholarpedia")<br/>
R. Quian Quiroga<br/>
Scholarpedia 2 (12): 3583. 2007

