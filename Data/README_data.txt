GENERAL INFORMATION
-------------------

1. Title of dataset:
	Information integration for nutritional decision-making in desert locusts

2. Author information:
	2.1.	Yannick Günzel
			University of Konstanz
			yannick.guenzel@uni-konstanz.de
		
	2.2.	Felix B. Oberhauser
			University of Konstanz
			felix.oberhauser@outlook.com
	
	2.3. 	Einat Couzin-Fuchs
			University of Konstanz
			einat.couzin@uni-konstanz.de

3. Date of data collection:
	2019-2020

4. Geographic location of data collection:
	Konstanz, Germany

5. Funding sources that supported the collection of the data:
	Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy – EXC 2117 – 42203798




DATA & FILE OVERVIEW
--------------------

1. Description of the dataset:
	We generated these data to investigate how desert locusts (Schistocerca gregaria) - either alone or in groups of varying sizes - accumulate evidence and integrate different classes of information for effective foraging.
	Raw video data have not been included as this would exceed space limitations. Instead, we provide normalized 2D trajectory data for each trial (for details, see data-specific information below).




METHODOLOGICAL INFORMATION & SCOPE 
----------------------------------
Swarms of the migratory desert locust can extend over several hundred square kilometres, and starvation compels this ancient pest to devour everything in its path. Theory suggests that gregarious behaviour benefits foraging efficiency over a wide range of spatial food distributions. However, despite the importance of identifying the processes by which swarms locate and select feeding sites to predict their progression, the role of social cohesion during foraging remains elusive. We investigated the evidence accumulation and information integration processes that underlie locusts' nutritional decision-making by employing a Bayesian formalism on high-resolution tracking data from foraging locusts. We tested individual gregarious animals and groups of different sizes in a 2-choice behavioural assay in which food patch qualities were either different or similar. We then predicted the decisions of individual locusts based on personally acquired and socially derived evidence by disentangling the relative contributions of each information class. Our study suggests that locusts balance incongruent evidence but reinforce congruent ones, resulting in more confident assessments when evidence aligns. We provide new insights into the interplay between personal experience and social context in locust foraging decisions which constitute a powerful empirical system to study local individual decisions and their consequent collective dynamics.




DATA-SPECIFIC INFORMATION:
--------------------------
	1. Tracking (folder)
	   Each trial's data with all tracked animals were saved as a CSV table and a corresponding MAT file containing additional information. File names give insight into the tested group size and experimental condition.
	   File names constructed as follows:  <group size><condition><date(yymmdd)>_* 
	   Thus, the pair 01EQ20191202_tracked.csv and 01EQ20191202_annotation.mat belonged to a trial with one animal (01) under the equal food condition (EQ) and were collected on February 2nd 2019 (20191202).
		
		1.1. _tracked.csv
			1.1.1. Number of variables	: 5
			1.1.2. Number of rows		: N*45000 frames (N=number of animals)
			1.1.3. Variable list:
				cnt			: row counter (0 - N*45000-1)
				frame			: corresponding video frame (blocks of N)
				pos_x			: animal's normalized x-position in video frame (blocks of N)
				pos_y			: animal's normalized y-position in video frame (blocks of N)
				id			: current animal ID (N unique entries)
			
		1.2. _annotation.mat
			1.2.1. Number of variables: 7
			1.2.2. Variable list:
				Arena			  : information on the arena
								- radius_cm (radius of the arena [cm]
				
				Comment			: information on which food patch had which quality
				
				PatchA			: information on the location and radius of food patch A (see "Comment" for respective quality)
								- loc (location [x,y])
								- radius
				
				PatchB			: information on the location and radius of food patch B (see "Comment" for respective quality)
								- loc (location [x,y])
								- radius				
				
				Tracker			: tracking software used. 
				
				Transformation		: information on the transformation applied to the data to normalize them (centring, scaling, rotating)
								- rotation_deg (how was the data rotated to align food patches horizontallz [degree])
								- px_per_cm (conversion factor between pixels and centimeters)
								- translation_px (translation of data to center them [x,y])
				
				Valid			: indicate whether everything is ok (Boolean)
		
		
	2. PatchWeights.csv (table)
	   Table with before and after weights of food patches.
		
		2.1. Number of variables		: 10
		2.2. Number of rows	  		: 109
		2.3. Variable list:
			id		  		: trial (<group size><condition><date>)
			date  				: date (yymmdd)
			A_loc	  			: position ID of patch A
			A_cond				: condition of patch A (HQ:1, LQ:-1)
			A_before			: weight of patch A [g] before beginning of trial
			A_after				: weight of patch A [g] after beginning of trial 
			B_loc	  			: position ID of patch A
			B_cond				: condition of patch B (HQ:1, LQ:-1)
			B_before			: weight of patch B [g] before beginning of trial
			B_after				: weight of patch B [g] after beginning of trial
		2.4. Missing data codes : NaN 
