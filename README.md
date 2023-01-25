# Content
ImageJ and Matlab code to plot quantitative average localization graphs (localization profiles, polar localization index, symmetry index, whole-cell fluorescence) of fluorescent proteins from (rod-shaped) bacteria segmented by BacStalk. 

# Requirements

  * [**ImageJ/Fiji**](https://fiji.sc/) - version 1.53 or higher
   * plugins contained in Fiji
   * [StackReg](http://bigwww.epfl.ch/thevenaz/stackreg/) plugin
   * [TurboReg](http://bigwww.epfl.ch/thevenaz/turboreg/) plugin
   * [MultiStackReg](https://biii.eu/multistackreg) plugin
 * [**MATLAB**](https://ch.mathworks.com/products/matlab.html) - version R2019b or higher
   * [Image Processing Toolbox](https://ch.mathworks.com/products/image.html)
   * [Parallel Computing Toolbox](https://ch.mathworks.com/products/parallel-computing.html?s_tid=srchtitle_Parallel%20Processing%20Toolbox_1)
 * [**BacStalk**](https://drescherlab.org/data/bacstalk/docs/index.html) - version 1.8

# Preparation

* Prepare folders for data storage and results as specified in the various scripts
  * One folder for each of your proteins of interest containing:
    * Input_Profile_Analysis.xlsx (what you enter here must match exactly how you name the BacStalk file, see below for details)
      * each line represents one set of data
      * **column 1** - strain number
      * **column 2** - name (whatever is in between "Dataset_" and "_condition"
      * **column 3** - date
      * **column 4** - condition
    * folder "data"
    * folder "mat"
    * folder "svg_files"
    * folder "fig_files"
  * jpg files of the plots are just saved in the protein of interest directory
  * all of the final folders must be contained in a parent folder with the name of your protein of interest, either all together or separated in different protein of interest folders (for example to separate the final graphs from the data)
      
* Make sure all directories in the ImageJ and Matlab scripts are matching your personal system's directories
* You can search for "directory" in the scripts and modify them, but make sure to only modify strings and not change variable names!


# How to run this code

* Run one of the SplitImage macros with [ImageJ/Fiji](https://fiji.sc/) to prepare for BacStalk analysis
  * **SplitImage_MJK.ijm** - for single images with first channel phase contrast + second channel fluorescence
  * **SplitImage_MJK_mNG+mScI.ijm** - for single images with first channel phase contrast + second channel fluorescence + third channel fluorescence (the Matlab code will only use one of the fluorescence channels at a time but you can run BacStalk with both channels together)
  * **..._firstPic_MJK.ijm** - same but for image sequences, only saves the first image

* Run [BacStalk](https://drescherlab.org/data/bacstalk/docs/index.html) to segment the cells
* Save BacStalk mat file
	* mat files must be saved in dir_data directory (see localization_profile_static.m)
  * the file name must follow the following format: **date_"DataSet"_number_condition** (e.g. 20220106_DataSet_923_mNG_PilG_sol0h)
  * file name elements must match the fields in Input_Profile_Analysis.xlsx

* Enter in the Input_Profile_Analysis.xlsx file the data you want to plot together
  * It is recommended to plot maximum 3 strains and 3 conditions at a time
  * If plotting more, the plots get a bit crowded

* Run localization_profile_static.m with [MATLAB](https://ch.mathworks.com/products/matlab.html) located in graph_plotting
  * change lines according to the protein (folder name) you want to plot and the name you gave the fluorescence channel
  * modify the lines containing "directory" (search for it) to match your system
