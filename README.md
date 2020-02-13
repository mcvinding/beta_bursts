# Beta burst analysis
General scripts and functions for time-domain analysis of beta bursts/events in MEG/(EEG).

**STILL WORK IN PROGRESS: NO WARRANTY OR GUARANTEE WHATSOEVER**

Supposed to do the following:
1. Find beta events in continuous data (e.g. resting state)
2. Find beta events in event-related data.

The functions will take a single channel band-pass filtered time-series to find bursts. Can be for either a continuous time-series or epochs.

All functions and example scripts are written in MATLAB and are/should be compatible with [FieldTrip](http://www.fieldtriptoolbox.org/).

## Warnings
Paths etc. might be completely messed up. Use with caution!

Example scripts are copy/paste from other projects. Filenames and examples are not generalizable. No compatibility guaranteed. Use with caution!

## Status of development:
* Finding beta events in continuous time-series (1) implemented for resting state (i.e. continuous data) in [PD_beta_burst](https://github.com/mcvinding/PD_beta_bursts). The scripts in this repo are still somewhat specific to this project.
* Finding beta events for event related data is implemented but not properly tested on real data.
* Supports methods for defining thresholds based on Shin et al. and Little et al. - add more if appropriate.
* Missing proper documentation.

## Possible future developments
* Gluing option in definfing beta bursts
* Functionality for other bands that beta; e.g. alpha fluctuations.

## Publications
Most of this was developed for the following publication:

Vinding, M. C., Tsitsi, P., Waldthaler, J., Oostenveld, R., Ingvar, M., Svenningsson, P., & Lundqvist, D. (2019). Reduction of spontaneous cortical beta bursts in Parkinson's disease is linked to symptom severity. *biorXiv.org*. https://doi.org/10.1101/775494

## Contact
For questions or suggestions

