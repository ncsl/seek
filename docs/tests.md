# Tests
1. test failed data
	- Macauley/Kristin to insert fake data
2. test patient
	- LA02
	- EFRI03
3. test:
	- test plotting brain regions
	- test plotting seeg channels
	- test getting brain_region -> seeg channels


1. Input
- take CT Images
- take MRI Image, brain envelope 
- 

Methods:
- map brain envelope to the CT space and mask the skull out
- applying simple threshold 
	- can get clusters of bright spots that appear for varying different thresholds
	- 
- computing distance between each dot with all other clusters
	- distance matrix
	- use another threshold to decide which are connected and not
	- loop through and check if they are in a line or not
	- now each grouping is a vector of channels that represent an electrode
(note that it does not handle the case when the channels are not put into a line)


2. Output
- output list of channels grouped by electrodes
- excel sheet that can group these electrodes
	- missing electrode when algorithm fails
	-> can loop back to read into the function again with an entry point
	