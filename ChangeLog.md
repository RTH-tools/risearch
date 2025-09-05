##v 2.1 - 2018-09-07
- Updated mismatched seeds option, NOTE: -m option changed behavior in comparison to previous version
- Added --noGUseed option to disable wobble pairs in the seed
- Added -p3 output option
##v 2.0 - 2017-01-13
- First official release

013-04-22  Anne Wenzel  <wenzel@rth.dk>
	* include energy threshold
	* added alternative output formats

020-03-29 Giulia Corsi <giulia@rth.dk>
	* added scoring matrices for RNA-DNA and DNA-DNA interactions
	* new option -f to force start and end of interactions
	* new option -w to weight stacking base pair contributions

021-06-15 Giulia Corsi <giulia@rth.dk>
	* added error if any of the options -f -s -n -l -e -p are used
	together with -f and -w
