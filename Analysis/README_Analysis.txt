GENERAL INFORMATION: Analysis
-----------------------------

1. Title of dataset:
	Information integration for decision-making in desert locusts

2. Author information:
	2.1.	Yannick Günzel
			University of Konstanz
			yannick.guenzel@uni-konstanz.de
			ORCID: 0000-0001-7553-4742
		
	2.2.	Felix B. Oberhauser
			University of Konstanz
			felix.oberhauser@outlook.com
			ORCID: 0000-0002-9278-2497
	
	2.3. 	Einat Couzin-Fuchs
			University of Konstanz
			einat.couzin@uni-konstanz.de
			ORCID: 0000-0001-5269-345X

3. Date of data collection:
	2019-2020

4. Geographic location of data collection:
	Konstanz, Germany

5. Funding sources that supported the collection of the data:
	This work was completed with the support of the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy – EXC 2117 – 422037984.

6. Data and code availability
	Analysis code and processed trajectories are available at https://doi.org/10.5281/zenodo.7541780. Custom-written GUI for manually supervised tracking is available at https://github.com/YannickGuenzel/BlobMaster3000




DATA & FILE OVERVIEW
--------------------

1. Description of the dataset:
	We generated these data to investigate how desert locusts (Schistocerca gregaria) - either alone or in groups of varying sizes - accumulate evidence and integrate different classes of information for effective foraging.
	Primary video data have not been included here as this would exceed space limitations. Instead, we provide normalized 2-D trajectory data for each trial (for details, see data-specific information below).

2. Description of data analysis:
	Data were analyzed and plotted using custom-written MATLAB (R2022a) scripts (a detailed description of each script's behavior can be found below). Files that contain several sub-functions for the scripts listed above and which must be within the same folder are SubFcn.m and SubFcn_LocustDecisionSystem.m.
	To reproduce our results, run the scripts in the following order:
	1. mainFcn_LocustFeeding.m
	2. mainFcn_LocustDecisionSystem.m
	3. mainFcn_LocustDecisionSystem_plot.m		
	4. mainFcn_LocustFeeding_statistics.m		
	5. mainFcn_ExampleVideos.m		
	



METHODOLOGICAL INFORMATION & SCOPE 
----------------------------------
Locust swarms can extend over several hundred kilometers, and starvation compels this ancient pest to devour everything in its path. Theory suggests that gregarious behavior benefits foraging efficiency, yet the role of social cohesion in locust foraging decisions remains elusive. To this end, we collected high-resolution tracking data of individual and grouped gregarious desert locusts in a 2-choice behavioral assay with animals deciding between patches of either similar or different quality. Carefully maintaining the animals' identities allowed us to monitor what each individual has experienced and to estimate the leaky accumulation process of personally acquired and, when available, socially derived evidence. We fitted these data to a model based on Bayesian estimation to gain insight into the locust social decision-making system for patch selection. By disentangling the relative contribution of each information class, our study suggests that locusts balance incongruent evidence but reinforce congruent ones. We provide insight into the collective foraging decisions of social (but non-eusocial) insects and present locusts as a powerful empirical system to study individual choices and their consequent collective dynamics.




FILE-SPECIFIC INFORMATION:
--------------------------
	1. ExampleVideo (folder)
	   This folder contains the raw data for one trial to create an example animation that shows how animals were tracked while maintaining their identities over time.
		
		1.1. 15UE20191206.mp4
		     Primary data from December 6th 2019. A video at 25fps with 15 animals foraging under the unequal patch conditions for 30 minutes. The high-quality food is placed in the top left corner of the arena. The low-quality food is placed in the bottom right corner of the arena.
		
		1.2. 15UE20191206_annotation.mat
			1.2.1. Number of variables: 7
			1.2.2. Variable list:
				Arena			  	: information on the arena
									- radius_cm (radius of the arena [cm]
				Comment				: information on which food patch had which quality
				PatchA				: information on the location and radius of food patch A (see "Comment" for respective quality)
									- loc (location [x,y])
									- radius
				PatchB				: information on the location and radius of food patch B (see "Comment" for respective quality)
									- loc (location [x,y])
									- radius				
				Tracker				: tracking software used (either TRex or BlobMaster3000). 
				Transformation		: information on the transformation applied to the data to normalize them (centring, scaling, rotating)
									- rotation_deg (how was the data rotated to align food patches horizontally [degree])
									- px_per_cm (conversion factor between pixels and centimeters)
									- translation_px (translation of data to center them [x,y])
				Valid				: indicate whether everything is ok (Boolean)
		
		1.3. 15UE20191206_tracked.csv
		     NOTE: tracking data has not been normalized and is in video coordinates.
			1.3.1. Number of variables	: 5
			1.3.2. Number of rows		: N*45000 frames (N=number of animals=15)
			1.3.3. Variable list:
				cnt					: row counter (0 - N*45000-1)
				frame				: corresponding video frame (blocks of N)
				pos_x				: animal's normalized x-position in video frame (blocks of N)
				pos_y				: animal's normalized y-position in video frame (blocks of N)
				id					: current animal ID (N unique entries)
	
	2. mainFcn_ExampleVideos.m
	   Script for producing a video exemplifying how animals were tracked while maintaining their identities over time
		
	3. mainFcn_LocustDecisionSystem.m
	   This script infers the locust decision-making rule used in patch choice based on a Bayesian formalism. It assumes that mainFcn_LocustFeeding.m has finished successfully and saved pooled data as PooledData.mat
	
	4. mainFcn_LocustDecisionSystem_plot.m
	   This script plots the results of inferring the locust decision-making rule and calculates statistics accordingly. It assumes that mainFcn_LocustDecisionSystem.m finished successfully. 
	
	5. mainFcn_LocustFeeding.m
	   This script pools trials corresponding to group size and experimental condition, performs all analysis steps, and exports raw figures as pdf vector graphics. We saved figures using a toolbox for exporting publication quality figures https://github.com/altmany/export_fig
	
	6. mainFcn_LocustFeeding_statistics.m
	   This script calculates statistics according to the figures created by mainFcn_LocustFeeding.m (assuming it finished successfully).
	
	7. SubFcn.m
	   This class contains several function that are being called during the execution of above-listed scripts.
		
	8. SubFcn_LocustDecisionSystem.m
	   This class contains several function that are being called during the execution of mainFcn_LocustDecisionSystem.m