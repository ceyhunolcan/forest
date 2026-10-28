# Rest-Activity Rhythms

## Overview
Use `oak.circadian_rhythms` to compute circadian rest-activity rhythm (RAR) summary
statistics from Beiwe accelerometer data. This complements the gait/step
output of `oak` (see [oak](oak.md)): the same raw accelerometer stream,
summarised as 24-hour activity patterning rather than walking.

## Usage
```python
from forest.oak.circadian_rhythms import run
from forest.constants import Frequency


# Determine study folder and output folder
study_folder = "project/data"
output_folder = "project/results"

# Determine study timezone and time frames for data analysis
tz_str = "America/New_York"
time_start = "2018-01-01 00_00_00"
time_end = "2022-01-01 00_00_00"

# Write the recording-level summary and a daily file for each participant.
frequency = Frequency.DAILY
users = None

# Call the main function
run(study_folder, output_folder, tz_str, frequency,
    time_start, time_end, users)
```
`Frequency.DAILY` is the default. `Frequency.HOURLY_AND_DAILY` is accepted
as an alias for these same two outputs; it does not produce hourly metrics.
Other `Frequency` values write only the recording-level summary. This module
uses the shared Forest enum, whose names describe other analyses too.

### Default tuning parameters
```python
# interval length in seconds; must divide a 24-hour day evenly - interval_seconds
# (this is the length of one interval, not the count per day: the number of
# intervals in a 24-hour day is 86_400 / interval_seconds, e.g. 1440 at 60 s)
interval_seconds = 60

# numeric value of 1 g in the data's units (1.0 if the accelerometer is
# already in g, consistent with oak; ~9.81 if the data are in m/s^2) - gravity
gravity = 1.0

# participants with fewer valid days than this are skipped - min_valid_days
min_valid_days = 3

# each counted day must have at least this many observed hours
min_valid_hours_per_day = 16.0
```

A day counts once only when its own observations reach 16 hours. Six days
with 12 hours each therefore do not pass the default three-day requirement.
The 16-hour default follows the reviewer's suggested inclusion convention;
it is adjustable, rather than a validated threshold for smartphone sensing.
For a near-complete-day protocol, use 24 hours; a protocol expecting shorter
observation may need a lower threshold and a sensitivity analysis.

Daily output has a separate, more permissive rule: a date needs at least 12
observed hours (`MIN_DAY_COVERAGE = 0.5`) to receive daily metrics. Dates below
that threshold retain their date and count, with metric cells left blank.
This 12-hour rule is a reporting choice, not an evidence-based inclusion
criterion. Passing the participant gate allows all its observed intervals
on ordinary 24-hour days to enter the recording-level metrics, including
intervals from days below the 16-hour inclusion threshold.

## Activity metric

Per sample, movement intensity is summarised as ENMO (the Euclidean norm of
the acceleration vector minus one gravitational unit, clipped at zero):
`ENMO = max(sqrt(x^2 + y^2 + z^2) / gravity - 1, 0)`. ENMO is averaged into
fixed intervals on a timezone-local grid, with missing intervals held as `NaN`.
The result is expressed in gravitational units after subtracting gravity:
`0.02` means an average excess magnitude of 0.02 g. It is a movement proxy;
it does not identify sleep, device wear, or physical activity categories.

An interval is observed if it has at least one valid sample. The coverage
count therefore measures availability of intervals, not continuous sampling
or verified wear time. A sparsely sampled minute counts the same as a densely
sampled minute. Inspect raw sample counts and gaps before using these summaries.

### Platform checks

The module does not infer operating system or acceleration units from the CSV.
Before pooling iOS and Android recordings, check units and resting magnitude,
sample rates, within-interval sample counts, missing intervals by clock time,
and ENMO distributions separately for each OS, app version, and OS version.
Use `gravity=1.0` for data in g and approximately `9.81` for data in m/s²;
process groups with different units separately. Synthetic unit-scaling tests
verify the arithmetic, but real iOS/Android equivalence has not been validated
by this PR. A phone left on a desk can produce an apparently regular rest
period that says little about the participant's activity.

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

### Worked output example

This example is synthetic, with three identical, completely observed days at
one-minute intervals. Activity follows `0.05 + 0.04*cos(2*pi*(hour-14)/24)`.
The following excerpt from `rar_summary.csv` is rounded for display:

```csv
backend_id,IS,IV,RA,L5,M10,L5_onset_h,M10_onset_h,cosinor_mesor,cosinor_amplitude,cosinor_acrophase_h,cosinor_r2,n_epochs,n_valid_epochs,n_days
example,1.000000,0.000019,0.722781,0.012795,0.079516,23.5,9.0,0.05,0.04,14.0,1.0,4320,4320,3.0
```

IS is 1 because every day repeats exactly. The fitted curve peaks at 14:00,
with a mean of 0.05 g and an amplitude of 0.04 g. L5 begins at 23:30 and
wraps across midnight; it is a five-hour average, so it exceeds the single
lowest point of the curve. M10 begins at 09:00. The low IV follows from the
smooth curve sampled every minute, not from a clinical cutoff.

For each of these three dates the daily file contains the same metric values
except IS and recording-level counts, which are omitted. Its count is 1440
rather than 4320. For example, selected daily columns read:

```csv
date,IV,RA,L5,M10,L5_onset_h,M10_onset_h,cosinor_acrophase_h,n_valid_epochs
2026-03-01,0.000019,0.722781,0.012795,0.079516,23.5,9.0,14.0,1440
```

The CSV column names `n_epochs` and `n_valid_epochs` are retained pending the
maintainer's answer about renaming output fields. They count intervals, with
no difference in meaning from `interval_seconds`. Likewise, `rar` in the
filenames abbreviates rest-activity rhythm; it does not indicate an archive.
`n_days` counts the padded calendar-day span, including missing dates and
dates with daylight-saving changes, rather than the number meeting the inclusion threshold.

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

Relative amplitude (RA) compares the busiest ten-hour mean with the quietest
five-hour mean. Larger values indicate a stronger contrast in the observed
data, but neither RA nor the L5/M10 onsets identify sleep. The cosinor MESOR
is the fitted mean, amplitude is half the fitted peak-to-trough range, and
acrophase is the fitted peak's local clock hour. A low cosinor R² means that
a single smooth 24-hour wave explains little of the recorded variability;
its peak time may then be uninformative.

Check coverage before interpreting an unusual value. Constant activity over
24 hours warrants inspecting the raw stream: a flat signal has no activity
variance, so IS and IV are undefined, while cosinor amplitude is near zero
and its phase has no useful interpretation. All-zero activity also makes RA
undefined. Blank CSV cells represent undefined metrics, not zero activity.
A missing nighttime region can prevent estimating the quietest window, or
force L5 onto a different fully observed window. Compare sampling intervals,
coverage patterns, and device units before comparing participants.

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

### Data-problem signals

Check coverage before interpreting unusual values. These patterns warrant inspecting the raw accelerometer data:

* Flat activity throughout the day. A constant signal may reflect a stationary phone or a recording problem. With no activity variance, IS and IV are undefined. Cosinor amplitude is approximately zero, and its reported peak time has no useful interpretation. All-zero activity also makes RA undefined; constant non-zero activity gives RA = 0.
* Near-zero activity with irregular fluctuations. A stationary phone can still record sensor noise. If that noise varies independently between intervals, IV can be near 2 and IS near `1/n_days`, about 0.14 over seven days, with cosinor amplitude and R² near zero. A synthetic noise example produces this pattern; it is not a validated way to identify a stationary phone. Check the raw stream before interpreting it as a participant’s irregular activity.
* A participant missing from `rar_summary.csv`. Check the processing logs. Participants with fewer than `min_valid_days` days meeting `min_valid_hours_per_day` are skipped with a warning. The defaults require at least three days with 16 observed hours each.
* Blank metric cells. These represent undefined or unavailable estimates, not zero activity. Daily metrics are blank for dates with fewer than 12 observed hours (`MIN_DAY_COVERAGE = 0.5`). L5, M10, their onset hours, and RA are also blank when either window lacks a complete candidate.
* Gaps at recurring clock times. Missing nighttime observations can exclude the true quietest window and shift L5 to another observed period. Recording-level windows use an average-day profile: each clock position needs an observation on at least one day, but there is no minimum number of contributing days per position.
* High coverage despite sparse sampling. An interval counts as observed with just one valid sample. Check raw sample counts and gaps; coverage does not establish continuous recording or that the participant carried the phone.
* Unexpected differences between recordings, including iOS and Android. Check acceleration units, the `gravity` setting, sampling intervals, and missingness patterns before interpreting differences as participant behaviour.

A low cosinor R² alone does not establish a data problem: the activity pattern may simply fit a single 24-hour wave poorly. Likewise, a strong apparent rhythm can reflect when someone carries or sets down their phone. These summaries describe recorded phone movement; they do not independently establish sleep or device wear.


## Assumptions and limitations

- **Daylight-saving changes.** The module uses local clock time
  and keeps a fixed 24-hour grid for each day. When clocks move
  forward, the missing hour remains `NaN`. When clocks move
  backward, observations from the repeated hour are averaged.
  Forest logs a warning when it detects these changes.

  DST dates are not automatically excluded. They follow the
  same coverage rules as other days, although missing or repeated
  hours can affect the resulting metrics. These dates should be
  checked when interpreting the results.

- **Missing observations.** IV uses only consecutive observed
  pairs and calculates variance from the available values.
  It does not bridge gaps. Each clock position in the average-day
  profile needs at least one observation across the recording.
  Positions without observations remain `NaN`.

  L5 and M10 are calculated only from windows where every clock
  position has data. If no complete window is available, the
  affected metrics are undefined. Missing observations can also
  shift the selected windows, especially when missingness is
  related to activity.

- **Local-midnight alignment.** Recording-level metrics assume
  that the series starts at local midnight. The `run()` function
  handles this alignment automatically. Anyone calling the
  metric functions directly must ensure the data is aligned
  correctly (see the `cosinor` docstring).

## Validation against existing packages

We compared Forest's calculations with ActCR 0.4.0 and CosinorPy 3.1
using six synthetic activity patterns. These included regular and noisy
cosine rhythms, a square wave, a square wave with naps, and white noise.
Each recording covered seven days and was tested at one-minute and
hourly resolution.

We used ActCR to check IS, IV, L5, M10, and RA, and CosinorPy to check
the four cosinor metrics.

Of the 108 comparisons, 100 agreed within our chosen tolerances.
Six differences came from how the packages handle windows crossing
midnight: Forest allows them, while ActCR searches within a single
day. We considered two white-noise acrophase results separately because
peak timing is difficult to interpret when there is little rhythmic
signal.

IS and IV agreed almost exactly. Cosinor mesor, amplitude, and R-squared
also matched closely, with a maximum acrophase difference of about
1.2 minutes.

These results support the calculations on complete synthetic data.
They do not establish how well the metrics perform with missing
observations, daylight-saving transitions, or real-world differences
in smartphone use.


## Implementation and provenance

These formulas build on the author's
[actrhythm library](https://github.com/ceyhunolcan/actrhythm).
The related NHANES analysis code is available in
[od-activity-rhythm](https://github.com/ceyhunolcan/od-activity-rhythm).

Applying these methods to NHANES does not mean that wearable MIMS
activity values and smartphone ENMO are directly comparable. We
also have not established whether iOS and Android recordings produce
equivalent measurements.

Direct implementation keeps the few required formulas readable and avoids
adding another activity data model or device reader. CosinorPy and pyActigraphy
provide broader analyses; wrapping one would be reasonable if Forest later
needs multiple cosinor components, diagnostic plots, or their wider feature
sets. This module provides a single 24-hour component and the listed summaries.

The module remains self-contained within Oak. Its input preparation, metric
functions, and runner can be separated along Oak's preprocess/analysis/runner
boundaries later; ENMO preparation is independent of the gait pipeline's
resampling and walking detection.

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
  https://doi.org/10.1038/s41598-022-24894-z
- Wendt A, Bielemann RM, Wehrmeister FC, et al. Is rest-activity rhythm
  prospectively associated with all-cause mortality in older people regardless
  of sleep and physical activity level? The 'Como Vai?' Cohort study.
  *PLoS ONE* 19(2):e0298031 (2024).
  https://doi.org/10.1371/journal.pone.0298031
