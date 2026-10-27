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

    average_day = compute_average_day_profile(arr, period)

    numerator = count * np.nansum((average_day - mean) ** 2)
    denominator = period * np.nansum((arr - mean) ** 2)

    if denominator == 0:
        return float("nan")

    return float(numerator / denominator)


def intradaily_variability(values: np.ndarray) -> float:
    """Measures how choppy activity is between neighboring intervals.

    Compares changes to overall variation and returns NaN if it can't be calculated.
    """
    arr = as_float_array(values)

    if arr.size < 2:
        return float("nan")

    observed = ~np.isnan(arr)
    valid_pairs = observed[1:] & observed[:-1]

    pair_count = np.count_nonzero(valid_pairs)
    observed_count = np.count_nonzero(observed)

    if pair_count == 0 or observed_count < 2:
        return float("nan")

    mean = np.nanmean(arr)
    differences = np.diff(arr)

    mean_squared_difference = np.nansum(differences ** 2) / pair_count
    variance = np.nansum((arr - mean) ** 2) / observed_count

    if variance == 0:
        return float("nan")

    return float(mean_squared_difference / variance)


def compute_average_day_profile(values: np.ndarray, intervals_per_day: int) -> np.ndarray:
    """Calculates the average activity at each position across days.

    Positions with no observed values stay as NaN.
    """
    arr = as_float_array(values)
    period = int(intervals_per_day)

    positions = np.arange(arr.size) % period
    position_means: list[float] = []

    for h in range(period):
        same_time = arr[positions == h]

        if np.any(~np.isnan(same_time)):
            position_means.append(float(np.nanmean(same_time)))
        else:
            position_means.append(np.nan)

    return np.asarray(position_means, dtype=float)


def compute_circular_window_mean(profile: np.ndarray, window: int) -> np.ndarray:
    """Averages values from each starting position, wrapping around the end of the day.

    Returns NaN if a window contains missing values.
    """
    arr = as_float_array(profile)
    window_size = max(1, int(window))

    offsets = np.arange(window_size)
    window_means: list[float] = []

    for start in range(arr.size):
        positions = (start + offsets) % arr.size
        values = arr[positions]

        window_means.append(float(np.mean(values)))

    return np.asarray(window_means, dtype=float)


@dataclass
class ActivityWindows:
    """Stores the least and most active periods of the average day.

    Attributes:
        l5: Average activity during the least active 5 hours.
        m10: Average activity during the most active 10 hours.
        l5_onset_hours: Clock hour when the least active period starts.
            For example, 23.5 means 23:30.
        m10_onset_hours: Clock hour when the most active period starts.
            For example, 8.5 means 08:30.
    """
    l5: float
    m10: float
    l5_onset_hours: float
    m10_onset_hours: float


def compute_l5_m10(values: np.ndarray, intervals_per_day: int) -> ActivityWindows:
    """Finds the least active 5 hours (L5) and most active 10 hours (M10).

    Uses the average day and returns NaN if the results can't be calculated.
    """
    period = int(intervals_per_day)
    intervals_per_hour = period / 24

    average_day = compute_average_day_profile(values, period)

    l5_window = round(L5_WINDOW_HOURS * intervals_per_hour)
    m10_window = round(M10_WINDOW_HOURS * intervals_per_hour)

    l5_means = compute_circular_window_mean(average_day, l5_window)
    m10_means = compute_circular_window_mean(average_day, m10_window)

    if np.all(np.isnan(l5_means)) or np.all(np.isnan(m10_means)):
        missing = float("nan")
        return ActivityWindows(
            l5=missing,
            m10=missing,
            l5_onset_hours=missing,
            m10_onset_hours=missing,
        )

    l5_start = int(np.nanargmin(l5_means))
    m10_start = int(np.nanargmax(m10_means))

    return ActivityWindows(
        l5=float(l5_means[l5_start]),
        m10=float(m10_means[m10_start]),
        l5_onset_hours=float(l5_start / intervals_per_hour),
        m10_onset_hours=float(m10_start / intervals_per_hour),
    )


def relative_amplitude(l5_value: float, m10_value: float) -> float:
    """Measures the contrast between the most and least active periods.

    The score ranges from 0 to 1 for nonnegative activity and returns
    NaN when the total is zero or a value is missing.
    """
    total = m10_value + l5_value

    if total == 0:
        return float("nan")

    return float((m10_value - l5_value) / total)


@dataclass
class CosinorFit:
    """Stores the results of fitting a wave to daily activity patterns.

    Attributes:
        mesor: Average level of the fitted wave.
        amplitude: How far the wave rises above its average level.
        acrophase_hours: Time in hours from the start of the series when
            the wave peaks. This matches local clock time when the series
            starts at local midnight, as run() ensures.
        r_squared: How much of the activity variation the fitted wave
            explains, from 0 to 1.
    """
    mesor: float
    amplitude: float
    acrophase_hours: float
    r_squared: float


def cosinor(
    values: np.ndarray,
    interval_seconds: float,
    period_hours: float = 24.0,
) -> CosinorFit:
    """Fits a cosine wave to activity data and measures how well it matches.

    Peak time is measured from the first value. It matches local clock time
    when run() aligns the series to local midnight.
    Returns all NaN values if fewer than three observations are available.
    """
    missing = float("nan")
    empty_fit = CosinorFit(
        mesor=missing,
        amplitude=missing,
        acrophase_hours=missing,
        r_squared=missing,
    )

    arr = as_float_array(values)
    hours = np.arange(arr.size) * interval_seconds / 3600

    observed = ~np.isnan(arr)
    observed_values = arr[observed]
    observed_hours = hours[observed]

    if observed_values.size < 3:
        return empty_fit

    omega = 2 * np.pi / period_hours

    matrix = np.column_stack([
        np.ones_like(observed_hours),
        np.cos(omega * observed_hours),
        np.sin(omega * observed_hours),
    ])

    coefficients, _, _, _ = np.linalg.lstsq(
        matrix, observed_values, rcond=None
    )

    mesor, beta, gamma = coefficients

    amplitude = np.hypot(beta, gamma)
    peak_time = (np.arctan2(gamma, beta) / omega) % period_hours

    fitted_values = matrix @ coefficients

    residual_sum = np.sum((observed_values - fitted_values) ** 2)
    total_sum = np.sum((observed_values - np.mean(observed_values)) ** 2)

    r_squared = float("nan") if total_sum == 0 else float(1 - residual_sum / total_sum)
    return CosinorFit(
        mesor=float(mesor),
        amplitude=float(amplitude),
        acrophase_hours=float(peak_time),
        r_squared=r_squared,
    )


def compute_rhythm_metrics(
    rows: list[dict[str, float | str]],
    backend_id: str,
    activity: np.ndarray,
    interval_seconds: float,
) -> None:
    """Calculates rhythm metrics and adds one summary row to the results.

    Args:
        rows: List where the new results will be added.
        backend_id: Identifier for the participant.
        activity: Activity values covering whole days, starting at local midnight.
        interval_seconds: Length of each activity interval in seconds.

    Raises:
        ValueError: If the interval length does not divide a day evenly.
    """
    intervals_per_day = SECONDS_IN_DAY / interval_seconds

    if not intervals_per_day.is_integer():
        raise ValueError(
            f"interval_seconds={interval_seconds} does not divide "
            "a 24-hour day evenly"
        )

    period = int(intervals_per_day)
    arr = as_float_array(activity)

    windows = compute_l5_m10(arr, period)
    fit = cosinor(arr, interval_seconds)

    rows.append({
        "backend_id": backend_id,
        "IS": interdaily_stability(arr, period),
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
        "n_valid_epochs": int(np.count_nonzero(~np.isnan(arr))),
        "n_days": round(arr.size / period, 3),
    })


def find_axis_columns(data_frame: pd.DataFrame) -> dict[str, str]:
    """Maps cleaned-up column names to the original ones.

    This makes it easier to find columns with different capitalization
    or extra spaces.
    """
    return {
        column.strip().lower(): column
        for column in data_frame.columns
    }


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
    values: np.ndarray,
    start: pd.Timestamp,
    interval_seconds: int,
) -> pd.DataFrame:
    """Calculates activity metrics for each day of the recording.

    Days with less than half their data observed keep only their date
    and valid interval count.

    Args:
        values: Activity values starting at local midnight.
        start: Starting date and time of the recording.
        interval_seconds: Length of each interval in seconds.

    Returns:
        A DataFrame containing the daily activity metrics.
    """
    arr = np.asarray(values, dtype=float).ravel()
    period = int(SECONDS_IN_DAY / interval_seconds)
    whole_days = arr.size // period

    rows: list[dict[str, float | str]] = []

    for day in range(whole_days):
        daily_values = arr[day * period:(day + 1) * period]
        date = (start + pd.DateOffset(days=day)).date().isoformat()

        observed_count = int(np.count_nonzero(~np.isnan(daily_values)))
        coverage = observed_count / period

        if coverage < MIN_DAY_COVERAGE:
            rows.append({
                "date": date,
                "n_valid_epochs": observed_count,
            })
            continue

        windows = compute_l5_m10(daily_values, period)
        fit = cosinor(daily_values, float(interval_seconds))

        rows.append({
            "date": date,
            "IV": intradaily_variability(daily_values),
            "RA": relative_amplitude(windows.l5, windows.m10),
            "L5": windows.l5,
            "M10": windows.m10,
            "L5_onset_h": windows.l5_onset_hours,
            "M10_onset_h": windows.m10_onset_hours,
            "cosinor_mesor": fit.mesor,
            "cosinor_amplitude": fit.amplitude,
            "cosinor_acrophase_h": fit.acrophase_hours,
            "cosinor_r2": fit.r_squared,
            "n_valid_epochs": observed_count,
        })

    columns = [
        "date",
        "IV",
        "RA",
        "L5",
        "M10",
        "L5_onset_h",
        "M10_onset_h",
        "cosinor_mesor",
        "cosinor_amplitude",
        "cosinor_acrophase_h",
        "cosinor_r2",
        "n_valid_epochs",
    ]

    return pd.DataFrame(rows, columns=columns)


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
    values: np.ndarray,
    intervals_per_day: int,
    min_valid_hours_per_day: float,
) -> int:
    """Counts how many days have enough observed activity data.

    Each day is checked separately, so six half-covered days
    aren't combined to count as full days.
    """
    arr = np.asarray(values, dtype=float).ravel()
    period = int(intervals_per_day)

    whole_days = arr.size // period

    if whole_days == 0:
        return 0

    daily_values = arr[:whole_days * period].reshape(whole_days, period)

    observed_per_day = np.count_nonzero(~np.isnan(daily_values), axis=1)
    observed_hours = observed_per_day * 24 / period

    valid_days = np.count_nonzero(observed_hours >= min_valid_hours_per_day)

    return int(valid_days)


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
