Dataset for Günzel et al. "Information integration for nutritional decision-making in desert locusts"

Abstract
--------
Swarms of the migratory desert locust can extend over several hundred square kilometres, and starvation compels this ancient pest to devour everything in its path. Theory suggests that gregarious behaviour benefits foraging efficiency over a wide range of spatial food distributions. However, despite the importance of identifying the processes by which swarms locate and select feeding sites to predict their progression, the role of social cohesion during foraging remains elusive. We investigated the evidence accumulation and information integration processes that underlie locusts' nutritional decision-making by employing a Bayesian formalism on high-resolution tracking data from foraging locusts. We tested individual gregarious animals and groups of different sizes in a 2-choice behavioural assay in which food patch qualities were either different or similar. We then predicted the decisions of individual locusts based on personally acquired and socially derived evidence by disentangling the relative contributions of each information class. Our study suggests that locusts balance incongruent evidence but reinforce congruent ones, resulting in more confident assessments when evidence aligns. We provide new insights into the interplay between personal experience and social context in locust foraging decisions which constitute a powerful empirical system to study local individual decisions and their consequent collective dynamics.


tracking_data
-------------
Each trial's data with all tracked animals is saved as <group size><condition><date(yymmdd)>_tracked.csv
with a corresponding annotation *annotation.mat stating the patchs' locations and arena's dimensions.
Note that data is already preprocessed by centring, scaling, and rotating the trajectories.


PatchWeights
------------
Table with before and after weights of food patches.