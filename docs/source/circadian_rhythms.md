# Rest-Activity Rhythms

## Overview
Use `oak.circadian_rhythms` to compute circadian rest-activity rhythm (RAR) summary
statistics from Beiwe accelerometer data. This complements the gait/step
output of `oak` (see [oak](oak.md)): the same raw accelerometer stream,
summarised as 24-hour activity patterning rather than walking.

## Usage:
```
from forest.oak.circadian_rhythms import run
from forest.constants import Frequency


# Determine study folder and output folder
study_folder = "project/data"
output_folder = "project/results"

# Determine study timezone and time frames for data analysis
tz_str = "America/New_York"
time_start = "2018-01-01 00_00_00"
time_end = "2022-01-01 00_00_00"

# Determine output resolution. Frequency.DAILY or Frequency.HOURLY_AND_DAILY
# additionally writes a per-participant daily file; any other value writes
# only the recording-level summary. See forest.constants.Frequency here:
# https://github.com/onnela-lab/forest/blob/develop/forest/constants.py
frequency = Frequency.HOURLY_AND_DAILY
users = None

# Call the main function
run(study_folder, output_folder, tz_str, frequency,
    time_start, time_end, users)
```
### Default tuning parameters:
```
# interval length in seconds; must divide a 24-hour day evenly - interval_seconds
# (this is the length of one interval, not the count per day: the number of
# intervals in a 24-hour day is 86_400 / interval_seconds, e.g. 1440 at 60 s)
interval_seconds = 60

# numeric value of 1 g in the data's units (1.0 if the accelerometer is
# already in g, consistent with oak; ~9.81 if the data are in m/s^2) - gravity
gravity = 1.0

# participants with fewer valid days than this are skipped - min_valid_days
min_valid_days = 3
```

## Activity metric

Per sample, movement intensity is summarised as ENMO (the Euclidean norm of
the acceleration vector minus one gravitational unit, clipped at zero):
`ENMO = max(sqrt(x^2 + y^2 + z^2) / gravity - 1, 0)`. ENMO is averaged into
fixed intervals on a timezone-local grid, with missing intervals held as `NaN`.

## Output

### Recording-level summary

The recording-level output is written to `rar_summary.csv` with one row per
participant. The following variables are created:

|     Variable            	|     Type     	|     Description of Variable                                                            	|
|-------------------------	|--------------	|----------------------------------------------------------------------------------------	|
|     backend_id          	|       str     	|     Beiwe backend identifier of the participant                                        	|
|     IS                  	|      float    	|     Interdaily stability, in (0, 1]; higher means a more reproducible daily pattern    	|
|     IV                  	|      float    	|     Intradaily variability, typically in [0, ~2]; higher means more fragmentation      	|
|     RA                  	|      float    	|     Relative amplitude, (M10 - L5) / (M10 + L5)                                         	|
|     L5                  	|      float    	|     Mean activity over the least-active 5 hours of the average day                     	|
|     M10                 	|      float    	|     Mean activity over the most-active 10 hours of the average day                     	|
|     L5_onset_h          	|      float    	|     Clock hour at which the L5 window starts                                           	|
|     M10_onset_h         	|      float    	|     Clock hour at which the M10 window starts                                          	|
|     cosinor_mesor       	|      float    	|     MESOR (rhythm-adjusted mean) of the 24 h cosinor fit                               	|
|     cosinor_amplitude   	|      float    	|     Amplitude of the 24 h cosinor fit                                                  	|
|     cosinor_acrophase_h 	|      float    	|     Clock hour of the fitted peak, in [0, 24)                                          	|
|     cosinor_r2          	|      float    	|     Coefficient of determination of the cosinor fit                                    	|
|     n_epochs            	|       int     	|     Number of intervals on the day-aligned grid                                           	|
|     n_valid_epochs      	|       int     	|     Number of non-missing intervals                                                    	|
|     n_days              	|      float    	|     Number of days spanned by the recording                                            	|

### Daily summary

When `frequency` is `Frequency.DAILY` or `Frequency.HOURLY_AND_DAILY`, a
per-participant `<backend_id>_rar_daily.csv` is also written with the
day-decomposable metrics (IS is recording-level only and is omitted):

|     Variable            	|     Type     	|     Description of Variable                                                	|
|-------------------------	|--------------	|----------------------------------------------------------------------------	|
|     date                	|       str     	|     Day of observation (yyyy-mm-dd)                                        	|
|     IV                  	|      float    	|     Intradaily variability for the day                                     	|
|     RA                  	|      float    	|     Relative amplitude for the day                                         	|
|     L5                  	|      float    	|     Least-active 5-hour mean for the day                                   	|
|     M10                 	|      float    	|     Most-active 10-hour mean for the day                                   	|
|     L5_onset_h          	|      float    	|     Clock hour at which the L5 window starts                               	|
|     M10_onset_h         	|      float    	|     Clock hour at which the M10 window starts                              	|
|     cosinor_mesor       	|      float    	|     MESOR of the day's cosinor fit                                         	|
|     cosinor_amplitude   	|      float    	|     Amplitude of the day's cosinor fit                                     	|
|     cosinor_acrophase_h 	|      float    	|     Clock hour of the day's fitted peak                                    	|
|     cosinor_r2          	|      float    	|     Coefficient of determination of the day's cosinor fit                  	|
|     n_valid_epochs      	|       int     	|     Number of non-missing intervals in the day                               	|

## Reading the output

Interdaily stability (IS) says how consistently the activity pattern repeats at
the same clock times across days. It sits near 1 when the same daily pattern
repeats almost exactly and near `1/n_days` when activity is independent random
noise, with no repeatable daily pattern.

Intradaily variability (IV) says how fragmented activity is within a day,
measured by changes between consecutive samples relative to overall
variability. White noise sits at about 2. A smooth 24-hour sinusoidal pattern
falls towards 0, and how far depends on the sampling interval: about 0.07 at
hourly intervals, almost nothing at one-minute intervals.

Lower IS and higher IV together mean a less consistent daily rhythm with more
abrupt fluctuations within each day.

These metrics have been studied in relation to mortality in prospective
cohorts. Xu et al. found that only one of the two held up in NHANES (n=7,252,
US adults, 674 deaths): IV associated with higher all-cause mortality, HR 1.13
(1.01-1.26) per quartile, fully adjusted, while IS was not significant. Wendt
et al. found the opposite split in a Brazilian cohort of older adults (n=965,
80 deaths): IS associated with lower mortality, HR 0.74 (0.55-0.98) per
standard deviation, while IV was not, and neither association survived
adjustment for physical activity. So the direction is consistent, more
disrupted rhythms track with worse outcomes, but which specific metric reaches
significance depends on the population and the covariates, and physical
activity can absorb the effect entirely. Treat IS and IV as descriptive rhythm
measures, not validated standalone mortality predictors.

## Assumptions and limitations

- **24-hour local days.** Days are treated as fixed 24-hour windows with a
  constant number of intervals per day (`86_400 / interval_seconds`). Days that are
  not 24 hours long -- daylight-saving transitions, which produce 23- or
  25-hour local days -- are not specially handled. On those specific days this
  can bias valid-day coverage counts, day alignment, and the timing of L5,
  M10, and the cosinor acrophase. This mirrors a wider assumption in the
  codebase and is a known limitation.
- **Local-midnight alignment.** Recording-level folding and the cosinor phase
  assume the series starts at local midnight; `run()` guarantees this via the
  day-grid alignment, but callers using the metric functions directly are
  responsible for that alignment (see the `cosinor` docstring).

## References

- Van Someren EJW, et al. Bright light therapy: improved sensitivity to its
  effects on rest-activity rhythms in Alzheimer patients by application of
  nonparametric methods. *Chronobiology International* 16(4):505-518 (1999).
- Gonçalves BSB, et al. Nonparametric methods in actigraphy: an update.
  *Sleep Science* 7(3):158-164 (2014).
- Cornelissen G. Cosinor-based rhythmometry. *Theoretical Biology and
  Medical Modelling* 11:16 (2014).
- Xu Y, Su S, Li X, Mansuri A, McCall WV, Wang X. Blunted rest-activity
  circadian rhythm increases the risk of all-cause, cardiovascular disease and
  cancer mortality in US adults. *Scientific Reports* 12:20665 (2022).
- Wendt A, Bielemann RM, Wehrmeister FC, et al. Is rest-activity rhythm
  prospectively associated with all-cause mortality in older people regardless
  of sleep and physical activity level? The 'Como Vai?' Cohort study.
  *PLoS ONE* 19(2):e0298031 (2024).
