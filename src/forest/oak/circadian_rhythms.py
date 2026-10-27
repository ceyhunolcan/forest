"""Rest-activity rhythm (RAR) metrics for the accelerometer (Oak) tree.

Derives circadian rest-activity summaries from Beiwe accelerometer data: interdaily stability (IS),
intradaily variability (IV), relative amplitude (RA) with L5/M10, and a single-component 24 h
cosinor fit (MESOR, amplitude, acrophase, R^2).

These complement oak's gait/step outputs: the same raw accelerometer stream, a different (24 h
patterning) question. Activity is summarised per fixed interval as ENMO -- the Euclidean norm of the
acceleration vector minus one gravitational unit, clipped at zero -- which is the standard
actigraphy proxy for movement intensity.

The public entry point ``run`` mirrors ``forest.oak.runners.run``. It writes a recording-level
``rar_summary.csv`` (one row per participant) and, when a daily frequency is requested,
per-participant ``<backend_id>_rar_daily.csv`` files holding the day-decomposable metrics.

References:
    Van Someren EJW, et al. Chronobiol Int 16(4):505-518 (1999).
    Goncalves BSB, et al. Sleep Sci 7(3):158-164 (2014).
    Cornelissen G. Theor Biol Med Model 11:16 (2014).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from os.path import isdir, join as pathjoin

import numpy as np
import pandas as pd

from forest.constants import Frequency, SECONDS_IN_DAY
from forest.utils import get_ids


logger = logging.getLogger(__name__)

ACCELEROMETER_DIR = "accelerometer"
# Coverage a single day needs before the daily file reports metrics for it. This is
# deliberately looser than run()'s min_valid_hours_per_day, which decides whether a
# participant's recording-level summary is reported at all.
MIN_DAY_COVERAGE = 0.5
L5_WINDOW_HOURS = 5     # the "5" in L5: width of the least-active window
M10_WINDOW_HOURS = 10   # the "10" in M10: width of the most-active window


def as_float_array(values: np.ndarray) -> np.ndarray:
    """Turns the input into a flat array of floats."""
    arr = np.asarray(values, dtype=float).ravel()

    if arr.size == 0:
        raise ValueError("Array is empty")

    return arr


def interdaily_stability(values: np.ndarray, intervals_per_day: int) -> float:
    """Measures how consistent activity patterns are across days.

    Returns NaN if it can't be calculated.
    """
    arr = as_float_array(values)
    period = int(intervals_per_day)

    count = np.sum(~np.isnan(arr))
    if count == 0:
        return float("nan")

    mean = np.nanmean(arr)
    positions = np.arange(arr.size) % period

    position_means: list[float] = []
    for h in range(period):
        same_time = arr[positions == h]

        if np.any(~np.isnan(same_time)):
            position_means.append(float(np.nanmean(same_time)))
        else:
            position_means.append(np.nan)

    average_day = np.asarray(position_means, dtype=float)

    numerator = count * np.nansum((average_day - mean) ** 2)
    denominator = period * np.nansum((arr - mean) ** 2)

    if denominator == 0:
        return float("nan")

    return float(numerator / denominator)


def intradaily_variability(activity: np.ndarray) -> float:
    """Intradaily variability (IV), typically in [0, ~2].

    IV is the mean squared difference between successive intervals divided by
    the variance of the series, so it asks how much the signal jumps from one
    interval to the next relative to how much it varies overall. Higher values
    mean a more fragmented rhythm: IV approaches ~2 for white noise and ~0 for
    a smooth sinusoid sampled densely. With no missing intervals this is the
    Van Someren form, ``n * sum (x_i - x_{i-1})^2 / ((n - 1) * sum (x_i - m)^2)``;
    the implementation below normalises over observed terms instead so that
    gaps do not bias the result.

    Args:
        activity: Per-interval activity vector; NaNs are omitted.

    Returns:
        The IV value, or NaN if it is undefined.
    """
    arr = as_float_array(activity)
    if arr.size < 2:
        return float("nan")

    valid = ~np.isnan(arr)

    # A successive difference is only defined when both of its endpoints are
    # present; np.diff across a gap yields NaN, which np.nansum then drops. To
    # keep the normalisation consistent we divide each sum by the count of
    # terms that actually contributed rather than by the raw interval count, so
    # missing intervals neither inflate nor deflate IV. With no gaps this reduces
    # exactly to the Van Someren n / (n - 1) form.
    valid_pairs = valid[1:] & valid[:-1]
    n_pairs = int(np.count_nonzero(valid_pairs))
    n_valid = int(np.count_nonzero(valid))
    if n_pairs < 1 or n_valid < 2:
        return float("nan")

    grand_mean = np.nanmean(arr)
    diffs = np.diff(arr)
    mean_sq_diff = np.nansum(diffs ** 2) / n_pairs
    variance = np.nansum((arr - grand_mean) ** 2) / n_valid
    return float(mean_sq_diff / variance) if variance > 0 else float("nan")


def compute_average_day_profile(activity: np.ndarray, intervals_per_day: int) -> np.ndarray:
    """Average one day's shape by folding every day onto the same clock grid.

    Each position in the returned profile is the mean activity at that point in
    the day across all the days in the recording. A position with no data on
    any day stays NaN rather than being filled in, because substituting an
    average there would make an unobserved part of the day look as active as
    the observed parts and would bias L5, M10 and the relative amplitude that
    are computed from this profile.

    Args:
        activity: Per-interval activity vector.
        intervals_per_day: Number of intervals in a 24 h day.

    Returns:
        Array of length intervals_per_day, NaN wherever a clock position was
        never observed.
    """
    arr = as_float_array(activity)
    period = int(intervals_per_day)
    # Fold every day onto one clock grid: pos gives each interval its position
    # within a day, so arr[pos == h] gathers clock position h across all days.
    pos = np.arange(arr.size) % period

    return np.array([
        np.nanmean(arr[pos == h]) if np.any(~np.isnan(arr[pos == h]))
        else np.nan
        for h in range(period)
    ])


def compute_circular_window_mean(
    profile: np.ndarray, window: int
) -> np.ndarray:
    """Mean over a length-``window`` window at every circular start."""
    width = max(1, int(window))
    extended = np.concatenate([profile, profile[: width - 1]]) if width > 1 else profile
    kernel = np.ones(width) / width
    return np.convolve(extended, kernel, mode="valid")


@dataclass
class ActivityWindows:
    """Least- and most-active window averages from an average-day profile.

    Args:
        l5: mean activity over the least-active five hours of the average day
        m10: mean activity over the most-active ten hours of the average day
        l5_onset_hours: clock hour at which the least-active window starts
        m10_onset_hours: clock hour at which the most-active window starts
    """
    l5: float
    m10: float
    l5_onset_hours: float
    m10_onset_hours: float


def compute_l5_m10(activity: np.ndarray, intervals_per_day: int) -> ActivityWindows:
    """Compute L5 and M10 from the average-day profile, with their onset hours.

    L5 and M10 are the conventional names for the mean activity over the
    least-active five hours and the most-active ten hours of the average day.
    A low L5 means the quietest part of the day is genuinely quiet, and a high
    M10 means the busiest part is genuinely busy, so the two together describe
    how far the day swings between rest and activity. The onset hours say when
    each of those windows begins.

    Args:
        activity: Per-interval activity vector.
        intervals_per_day: Number of intervals in a 24 h day.

    Returns:
        An ActivityWindows holding both window means and their onset hours.
    """
    period = int(intervals_per_day)
    intervals_per_hour = period / 24.0
    profile = compute_average_day_profile(activity, period)

    # Each window is slid around the clock, so L5 is the minimum of every
    # 5-hour window mean and M10 the maximum of every 10-hour window mean.
    means_5 = compute_circular_window_mean(
        profile, round(L5_WINDOW_HOURS * intervals_per_hour)
    )
    means_10 = compute_circular_window_mean(
        profile, round(M10_WINDOW_HOURS * intervals_per_hour)
    )
    # A window that overlaps an unobserved clock position comes out NaN from the
    # convolution above, so only fully observed windows are candidates here. If
    # the recording is gappy enough that no window is complete, every metric is
    # left undefined rather than reported from partial coverage.
    if np.all(np.isnan(means_5)) or np.all(np.isnan(means_10)):
        return ActivityWindows(
            l5=float("nan"),
            m10=float("nan"),
            l5_onset_hours=float("nan"),
            m10_onset_hours=float("nan"),
        )

    # nanargmin/nanargmax rather than argmin/argmax: the plain versions would
    # return the position of a NaN window instead of skipping it.
    idx_5 = int(np.nanargmin(means_5))
    idx_10 = int(np.nanargmax(means_10))

    return ActivityWindows(
        l5=float(means_5[idx_5]),
        m10=float(means_10[idx_10]),
        l5_onset_hours=float(idx_5 / intervals_per_hour),
        m10_onset_hours=float(idx_10 / intervals_per_hour),
    )


def relative_amplitude(l5_value: float, m10_value: float) -> float:
    """Relative amplitude ``(M10 - L5) / (M10 + L5)``, in [0, 1]."""
    denom = m10_value + l5_value
    if denom == 0:
        return float("nan")
    return float((m10_value - l5_value) / denom)


@dataclass
class CosinorFit:
    """Result of a single-component cosinor fit to an activity series.

    Args:
        mesor: rhythm-adjusted mean of the fitted curve
        amplitude: half the peak-to-trough difference of the fitted curve
        acrophase_hours: time of the fitted peak, in hours measured from the
            start of the series; a local clock hour only when the series
            starts at local midnight, as run() ensures
        r_squared: fraction of the variance in the activity series explained
            by the fit
    """
    mesor: float
    amplitude: float
    acrophase_hours: float
    r_squared: float


def cosinor(
    activity: np.ndarray, interval_seconds: float, period_hours: float = 24.0
) -> CosinorFit:
    """Single-component cosinor fit ``M + A cos(w t - phi)`` via OLS.

    The phase is measured from the first element of ``activity``: ``t = 0`` is
    index 0, so ``cosinor_acrophase_h`` is the peak time relative to the start
    of the series. It equals a local *clock* hour only when the series begins
    at local midnight, which ``run()`` guarantees by aligning to a
    midnight-anchored day grid before calling this. A caller passing an
    unaligned series gets a start-relative phase, not a clock hour. Days that
    are not 24 h (DST transitions) are not yet handled; see the module notes.

    Args:
        activity: Per-interval activity vector; NaNs are omitted.
        interval_seconds: Interval length in seconds.
        period_hours: Rhythm period; defaults to 24 h.

    Returns:
        A CosinorFit. Every field is NaN when fewer than three intervals are
        present, since the three-parameter fit is undefined below that.
    """
    nan_result = CosinorFit(float("nan"), float("nan"), float("nan"), float("nan"))
    arr = as_float_array(activity)
    hours = np.arange(arr.size) * (float(interval_seconds) / 3600.0)
    valid = ~np.isnan(arr)
    if int(valid.sum()) < 3:
        return nan_result
    times = hours[valid]
    values = arr[valid]
    omega = 2.0 * np.pi / float(period_hours)
    design = np.column_stack([
        np.ones_like(times), np.cos(omega * times), np.sin(omega * times)
    ])
    coef, _res, _rank, _sv = np.linalg.lstsq(design, values, rcond=None)
    mesor, beta, gamma = coef
    amplitude = float(np.hypot(beta, gamma))
    acrophase = float((np.arctan2(gamma, beta) / omega) % period_hours)
    resid = values - design @ coef
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((values - values.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return CosinorFit(
        mesor=float(mesor),
        amplitude=amplitude,
        acrophase_hours=acrophase,
        r_squared=float(r_squared),
    )


def compute_rhythm_metrics(
    rows: list[dict[str, float | str]],
    backend_id: str,
    activity: np.ndarray,
    interval_seconds: float,
) -> None:
    """Full recording-level RAR metric set.

    Args:
        activity: Per-interval activity vector spanning whole days.
        interval_seconds: Interval length in seconds; ``86400 / interval_seconds``
            must be a whole number.

    Returns:
        Dict of IS, IV, RA, L5/M10 (+ onsets), cosinor terms and counts.
    """
    intervals_per_day = SECONDS_IN_DAY / float(interval_seconds)
    if not float(intervals_per_day).is_integer():
        raise ValueError(f"interval_seconds={interval_seconds} does not divide a 24 h day")
    intervals_per_day = int(intervals_per_day)
    arr = as_float_array(activity)
    windows = compute_l5_m10(arr, intervals_per_day)
    fit = cosinor(arr, interval_seconds)

    rows.append({
        "backend_id": backend_id,
        "IS": interdaily_stability(arr, intervals_per_day),
        "IV": intradaily_variability(arr),
        "RA": relative_amplitude(windows.l5, windows.m10),
        "L5": windows.l5,
        "M10": windows.m10,
        "L5_onset_h": windows.l5_onset_hours,
        "M10_onset_h": windows.m10_onset_hours,
        "cosinor_mesor": fit.mesor,
        "cosinor_amplitude": fit.amplitude,
        "cosinor_acrophase_h": fit.acrophase_hours,
        "cosinor_r2": fit.r_squared,
        "n_epochs": int(arr.size),
        "n_valid_epochs": int(np.sum(~np.isnan(arr))),
        "n_days": round(arr.size / intervals_per_day, 3),
    })


def find_axis_columns(data_frame: pd.DataFrame) -> dict[str, str]:
    """Look up a frame's columns by a normalised version of their names.

    Maps the lower-cased, stripped name to the column's real name, so callers
    can ask for "x" or "utc time" without depending on how the header was
    capitalised or spaced in the file they were handed, and then use the value
    to index the frame.
    """
    return {c.strip().lower(): c for c in data_frame.columns}


def build_utc_index(data_frame: pd.DataFrame) -> pd.DatetimeIndex:
    """Build a UTC DatetimeIndex from a Beiwe accelerometer frame.

    Beiwe accelerometer files use the same header on both platforms and carry
    a numeric "timestamp" in epoch milliseconds as well as a human-readable
    "UTC time" string. The numeric column is preferred because it needs no
    date parsing; the string column is a fallback for a file that is missing
    it. A file with neither cannot be placed on a timeline, so that raises
    rather than returning something empty.
    """
    columns = find_axis_columns(data_frame)

    if "timestamp" in columns:
        series = pd.to_datetime(data_frame[columns["timestamp"]], unit="ms", utc=True)
    elif "utc time" in columns:
        series = pd.to_datetime(data_frame[columns["utc time"]], utc=True)
    else:
        raise ValueError("no 'timestamp' or 'utc time' column")

    return pd.DatetimeIndex(series)


def read_accelerometer_file(filepath: str) -> pd.DataFrame:
    """Read one raw accelerometer CSV, validating required columns."""
    data_frame = pd.read_csv(filepath)
    columns = find_axis_columns(data_frame)
    if not {"x", "y", "z"}.issubset(columns):
        raise ValueError(f"{filepath}: missing x/y/z columns")
    if "timestamp" not in columns and "utc time" not in columns:
        raise ValueError(f"{filepath}: missing timestamp column")
    data_frame.columns = [c.strip().lower() for c in data_frame.columns]
    return data_frame


def compute_enmo_intervals(
    data_frame: pd.DataFrame, interval_seconds: int, tz_str: str, gravity: float = 1.0
) -> pd.Series:
    """Mean ENMO per fixed interval on a gap-filled, tz-local grid.

    ENMO is ``max(|a| / gravity - 1, 0)`` per sample. Set ``gravity`` to the
    numeric value of 1 g in the data's units (1.0 for g, ~9.81 for m/s^2).

    Args:
        data_frame: Raw Beiwe accelerometer rows (``timestamp`` in ms UTC or
            a ``UTC time`` column, plus ``x``/``y``/``z``; case-insensitive).
        interval_seconds: Interval length in seconds.
        tz_str: IANA timezone for localisation, e.g. ``America/New_York``.
        gravity: Value of 1 g in the data's units.

    Returns:
        Activity Series indexed by tz-aware interval-start timestamps, spanning
        the first to last observed interval with missing intervals set to NaN.
    """
    if data_frame.empty:
        return pd.Series(dtype=float)
    columns = find_axis_columns(data_frame)
    for axis in ("x", "y", "z"):
        if axis not in columns:
            raise ValueError(f"missing accelerometer axis '{axis}'")
    x = data_frame[columns["x"]].to_numpy(dtype=float)
    y = data_frame[columns["y"]].to_numpy(dtype=float)
    z = data_frame[columns["z"]].to_numpy(dtype=float)
    magnitude = np.sqrt(x ** 2 + y ** 2 + z ** 2) / float(gravity)
    enmo = np.clip(magnitude - 1.0, 0.0, None)
    series = pd.Series(enmo, index=build_utc_index(data_frame)).sort_index()
    series = series.tz_convert(tz_str)
    rule = f"{int(interval_seconds)}s"
    interval_means = series.resample(rule).mean()
    grid = pd.date_range(interval_means.index.min(), interval_means.index.max(), freq=rule, tz=tz_str)
    return interval_means.reindex(grid)


def parse_time_bound(value: str | None, tz_str: str) -> pd.Timestamp | None:
    """Parse a Beiwe ``%Y-%m-%d %H_%M_%S`` bound into a tz-aware stamp."""
    if value is None:
        return None
    naive = pd.to_datetime(value, format="%Y-%m-%d %H_%M_%S")
    return naive.tz_localize(tz_str)


def align_to_day_grid(activity: pd.Series, interval_seconds: int) -> tuple[np.ndarray, pd.Timestamp]:
    """Pad to local-midnight-aligned whole days for interval-of-day folding."""
    # Assumes 24-hour local days: intervals_per_day is fixed regardless of DST, so
    # 23-/25-hour transition days are not specially handled (known limitation,
    # see circadian_rhythms.md).
    intervals_per_day = SECONDS_IN_DAY // int(interval_seconds)
    index = pd.DatetimeIndex(activity.index)
    start = index[0].normalize()
    end = index[-1].normalize() + pd.Timedelta(days=1)
    rule = f"{int(interval_seconds)}s"
    grid = pd.date_range(start, end, freq=rule, tz=index.tz, inclusive="left")
    values = activity.reindex(grid).to_numpy(dtype=float)
    usable = (values.size // intervals_per_day) * intervals_per_day
    return values[:usable], start


def load_participant_activity(
    study_folder: str,
    backend_id: str,
    interval_seconds: int,
    tz_str: str,
    gravity: float,
    time_start: str | None,
    time_end: str | None,
) -> pd.Series:
    """Read, concatenate and interval-average one participant's accelerometer data."""
    acc_dir = pathjoin(study_folder, backend_id, ACCELEROMETER_DIR)
    if not isdir(acc_dir):
        return pd.Series(dtype=float)
    files = sorted(f for f in os.listdir(acc_dir) if f.endswith(".csv"))
    frames: list[pd.DataFrame] = []
    for filename in files:
        try:
            frames.append(read_accelerometer_file(pathjoin(acc_dir, filename)))
        except (ValueError, pd.errors.EmptyDataError, pd.errors.ParserError):
            logger.warning("skipping unreadable file %s", filename)
    if not frames:
        return pd.Series(dtype=float)
    raw = pd.concat(frames, ignore_index=True)
    activity = compute_enmo_intervals(raw, interval_seconds, tz_str, gravity)
    lower = parse_time_bound(time_start, tz_str)
    upper = parse_time_bound(time_end, tz_str)
    if lower is not None:
        activity = activity[activity.index >= lower]
    if upper is not None:
        activity = activity[activity.index <= upper]
    return activity


def compute_daily_metrics(
    values: np.ndarray, start: pd.Timestamp, interval_seconds: int
) -> pd.DataFrame:
    """Per-day day-decomposable metrics (IS is recording-level only)."""
    intervals_per_day = SECONDS_IN_DAY // int(interval_seconds)
    n_days = values.size // intervals_per_day
    rows: list[dict[str, float | str]] = []
    for day in range(n_days):
        chunk = values[day * intervals_per_day:(day + 1) * intervals_per_day]
        date = (start + pd.Timedelta(days=day)).date().isoformat()
        valid_fraction = np.sum(~np.isnan(chunk)) / intervals_per_day
        if valid_fraction < MIN_DAY_COVERAGE:
            rows.append({"date": date, "n_valid_epochs": int(np.sum(~np.isnan(chunk)))})
            continue
        windows = compute_l5_m10(chunk, intervals_per_day)
        fit = cosinor(chunk, interval_seconds)

        row: dict[str, float | str] = {
            "date": date,
            "IV": intradaily_variability(chunk),
            "RA": relative_amplitude(windows.l5, windows.m10),
            "L5": windows.l5,
            "M10": windows.m10,
            "L5_onset_h": windows.l5_onset_hours,
            "M10_onset_h": windows.m10_onset_hours,
            "cosinor_mesor": fit.mesor,
            "cosinor_amplitude": fit.amplitude,
            "cosinor_acrophase_h": fit.acrophase_hours,
            "cosinor_r2": fit.r_squared,
            "n_valid_epochs": int(np.sum(~np.isnan(chunk))),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def get_ids_with_accelerometer_data(study_folder: str) -> list[str]:
    """List the participant IDs in a study folder that have accelerometer data.

    Builds on forest.utils.get_ids, so the registry folder and hidden folders
    are excluded, then keeps only the participants who actually have an
    accelerometer subfolder to read.
    """
    return sorted(
        backend_id for backend_id in get_ids(study_folder)
        if isdir(pathjoin(study_folder, backend_id, ACCELEROMETER_DIR))
    )


def count_valid_days(
    values: np.ndarray, intervals_per_day: int, min_valid_hours_per_day: float
) -> int:
    """Count the days that individually have enough data to be usable.

    A day counts only if that day on its own reaches the required number of
    observed hours. Totalling the observed intervals across the whole
    recording and dividing by a day's length would let several half-observed
    days stand in for one whole one, which is not the same thing: six days at
    twelve hours each would pass a three-day requirement without a single day
    being well covered.

    Args:
        values: Per-interval activity on a whole number of day-aligned days.
        intervals_per_day: Number of intervals in a 24 h day.
        min_valid_hours_per_day: Hours a day needs before it counts.

    Returns:
        The number of days meeting the requirement.
    """
    n_days = values.size // intervals_per_day
    if n_days == 0:
        return 0

    days = values[: n_days * intervals_per_day].reshape(n_days, intervals_per_day)
    coverage_per_day = np.sum(~np.isnan(days), axis=1) / intervals_per_day
    required_coverage = float(min_valid_hours_per_day) / 24.0

    return int(np.sum(coverage_per_day >= required_coverage))


def run(
    study_folder: str,
    output_folder: str,
    tz_str: str,
    frequency: Frequency = Frequency.DAILY,
    time_start: str | None = None,
    time_end: str | None = None,
    users: list[str] | None = None,
    interval_seconds: int = 60,
    gravity: float = 1.0,
    min_valid_days: float = 3.0,
    min_valid_hours_per_day: float = 16.0,
) -> pd.DataFrame:
    """Compute rest-activity rhythm summaries over a Beiwe study folder.

    Args:
        study_folder: Path with one subfolder per backend ID, each holding
            an ``accelerometer`` directory of raw CSVs.
        output_folder: Directory for outputs (created if absent).
        tz_str: IANA timezone of the study site.
        frequency: Output resolution; ``Frequency.DAILY`` or
            ``Frequency.HOURLY_AND_DAILY`` additionally writes per-day files.
        time_start: Optional lower bound, ``%Y-%m-%d %H_%M_%S``.
        time_end: Optional upper bound, ``%Y-%m-%d %H_%M_%S``.
        users: Optional explicit backend IDs; otherwise all are used.
        interval_seconds: Interval length in seconds; must divide 86_400 evenly.
        gravity: Value of 1 g in the data's units (1.0 = data already in g).
        min_valid_days: Participants with fewer valid days are skipped.
        min_valid_hours_per_day: Hours of data a single day needs before it
            counts as one of those valid days. The default of 16 h follows the
            usual wear-time convention in accelerometry studies, which ask for
            at least 16 hours a day over at least three or four days. Lower it
            for protocols with shorter expected wear, raise it if you want only
            near-complete days to count.

    Returns:
        The recording-level summary DataFrame (also written to disk).
    """
    os.makedirs(output_folder, exist_ok=True)
    if users is None:
        users = get_ids_with_accelerometer_data(study_folder)
    intervals_per_day = SECONDS_IN_DAY // int(interval_seconds)
    write_daily = frequency in (Frequency.DAILY, Frequency.HOURLY_AND_DAILY)
    rows: list[dict[str, float | str]] = []
    for backend_id in users:
        activity = load_participant_activity(
            study_folder, backend_id, interval_seconds, tz_str, gravity, time_start, time_end
        )
        if activity.empty:
            continue
        values, start = align_to_day_grid(activity, interval_seconds)
        valid_days = count_valid_days(values, intervals_per_day, min_valid_hours_per_day)
        if valid_days < min_valid_days:
            logger.warning(
                "skipping %s: %.1f valid days (< %.1f)", backend_id, valid_days, min_valid_days
            )
            continue
        compute_rhythm_metrics(rows, backend_id, values, interval_seconds)
        if write_daily:
            daily = compute_daily_metrics(values, start, interval_seconds)
            daily.to_csv(pathjoin(output_folder, f"{backend_id}_rar_daily.csv"), index=False)
    summary = pd.DataFrame(rows)
    if summary.empty:
        # No participant met the coverage threshold; still return a frame with the identifier column
        # so callers can rely on its presence.
        summary = pd.DataFrame(columns=["backend_id"])
    summary.to_csv(pathjoin(output_folder, "rar_summary.csv"), index=False)
    return summary
