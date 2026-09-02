"""Ground-truth tests for forest.oak.circadian_rhythms (rest-activity rhythms).

Each test builds a synthetic signal with known RAR properties and checks that the estimators recover
them, plus an end-to-end ``run`` over a small synthetic study folder.

Original Author: Ceyhun Olcan
"""

import numpy as np
import pandas as pd
import pytest

from forest.constants import Frequency
from forest.oak.circadian_rhythms import (compute_average_day_profile, count_valid_days,
    compute_enmo_intervals, compute_l5_m10, compute_rhythm_metrics, cosinor,
    interdaily_stability, intradaily_variability, run)


INTERVAL_SECONDS = 60
INTERVALS_PER_DAY = 86_400 // INTERVAL_SECONDS


@pytest.fixture
def rng():
    return np.random.default_rng(20260709)


def _hours(n_days):
    return np.arange(n_days * INTERVALS_PER_DAY) * (INTERVAL_SECONDS / 3600.0)


def test_strong_rhythm_recovers_cosinor(rng):
    times = _hours(7)
    omega = 2 * np.pi / 24
    activity = 50 + 40 * np.cos(omega * (times - 14))
    activity = np.clip(activity + rng.normal(0, 2, times.size), 0, None)
    rows: list[dict[str, float | str]] = []
    compute_rhythm_metrics(rows, "test", activity, INTERVAL_SECONDS)
    metrics = rows[0]
    assert metrics["IS"] > 0.9
    assert metrics["IV"] < 0.05
    assert metrics["RA"] > 0.6
    assert abs(metrics["cosinor_mesor"] - 50) < 1.5
    assert abs(metrics["cosinor_amplitude"] - 40) < 1.5
    assert abs(metrics["cosinor_acrophase_h"] - 14) < 0.3
    assert 8 <= metrics["M10_onset_h"] <= 12


def test_white_noise_is_at_floor(rng):
    activity = np.clip(rng.normal(50, 15, 7 * INTERVALS_PER_DAY), 0, None)
    rows: list[dict[str, float | str]] = []
    compute_rhythm_metrics(rows, "test", activity, INTERVAL_SECONDS)
    metrics = rows[0]
    assert abs(metrics["IS"] - 1.0 / 7.0) < 0.05
    assert abs(metrics["IV"] - 2.0) < 0.15
    assert metrics["RA"] < 0.1
    assert metrics["cosinor_r2"] < 0.02


def test_identical_days_give_is_one():
    one_day = np.abs(np.sin(np.linspace(0, 2 * np.pi, INTERVALS_PER_DAY, endpoint=False))) * 30
    activity = np.tile(one_day, 5)
    assert abs(interdaily_stability(activity, INTERVALS_PER_DAY) - 1.0) < 1e-6


def test_cosinor_flat_signal_is_degenerate():
    flat = np.full(3 * INTERVALS_PER_DAY, 5.0)
    result = cosinor(flat, INTERVAL_SECONDS)
    assert result.amplitude < 1e-6


def test_intradaily_variability_needs_two_points():
    assert np.isnan(intradaily_variability(np.array([1.0])))


def test_enmo_adapter_separates_rest_and_activity(rng):
    fs, seconds = 10, 3600
    start = pd.Timestamp("2026-03-01 00:00:00", tz="UTC")
    count = fs * seconds
    t_rest = start + pd.to_timedelta(np.arange(count) / fs, unit="s")
    rest = pd.DataFrame({
        "timestamp": t_rest.astype("int64") // 1_000_000,
        "x": rng.normal(0, 0.02, count),
        "y": rng.normal(0, 0.02, count),
        "z": 1 + rng.normal(0, 0.02, count),
    })
    t_active = t_rest + pd.Timedelta(hours=1)
    active = pd.DataFrame({
        "timestamp": t_active.astype("int64") // 1_000_000,
        "x": rng.normal(0, 0.5, count),
        "y": rng.normal(0, 0.5, count),
        "z": 1 + rng.normal(0, 0.5, count),
    })
    raw = pd.concat([rest, active], ignore_index=True)
    interval_means = compute_enmo_intervals(raw, INTERVAL_SECONDS, "UTC", gravity=1.0)
    assert bool(np.all(interval_means.dropna() >= 0))
    assert 118 <= int(interval_means.notna().sum()) <= 121
    rest_mean = float(np.nanmean(interval_means.iloc[:60]))
    active_mean = float(np.nanmean(interval_means.iloc[60:120]))
    assert active_mean > 5 * rest_mean


def test_run_end_to_end(tmp_path, rng):
    fs, days = 2, 4
    start = pd.Timestamp("2026-03-01 00:00:00", tz="America/New_York").tz_convert("UTC")
    count = fs * 86400 * days
    times = start + pd.to_timedelta(np.arange(count) / fs, unit="s")
    local_hour = (
        times.tz_convert("America/New_York").hour + times.tz_convert("America/New_York").minute / 60
    ).to_numpy()
    amp = 0.05 + 0.6 * np.clip(np.cos(2 * np.pi * (local_hour - 14) / 24), 0, None)
    frame = pd.DataFrame({
        "timestamp": times.astype("int64") // 1_000_000,
        "x": rng.normal(0, amp),
        "y": rng.normal(0, amp),
        "z": 1 + rng.normal(0, amp),
    })
    acc_dir = tmp_path / "study" / "participantA" / "accelerometer"
    acc_dir.mkdir(parents=True)
    per_day = fs * 86400
    for day in range(days):
        chunk = frame.iloc[day * per_day:(day + 1) * per_day]
        chunk.to_csv(acc_dir / f"2026-03-0{day + 1} 00_00_00.csv", index=False)
    
    out_dir = tmp_path / "out"
    summary = run(
        str(tmp_path / "study"),
        str(out_dir),
        "America/New_York",
        frequency=Frequency.HOURLY_AND_DAILY,
        interval_seconds=60,
    )
    assert (out_dir / "rar_summary.csv").exists()
    assert (out_dir / "participantA_rar_daily.csv").exists()
    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["backend_id"] == "participantA"
    assert abs(row["cosinor_acrophase_h"] - 14) < 0.5
    assert row["IS"] > 0.8
    daily = pd.read_csv(out_dir / "participantA_rar_daily.csv")
    assert len(daily) == days
    assert "IS" not in daily.columns  # IS is recording-level only


def _write_participant(acc_dir, start_local, n_days, fs, rng):
    """Write n_days of synthetic accelerometer CSVs for one participant."""
    acc_dir.mkdir(parents=True)
    start = pd.Timestamp(start_local, tz="America/New_York").tz_convert("UTC")
    count = fs * 86400 * n_days
    times = start + pd.to_timedelta(np.arange(count) / fs, unit="s")
    local_hour = (
        times.tz_convert("America/New_York").hour + times.tz_convert("America/New_York").minute / 60
    ).to_numpy()
    amp = 0.05 + 0.6 * np.clip(np.cos(2 * np.pi * (local_hour - 14) / 24), 0, None)
    frame = pd.DataFrame({
        "timestamp": times.astype("int64") // 1_000_000,
        "x": rng.normal(0, amp),
        "y": rng.normal(0, amp),
        "z": 1 + rng.normal(0, amp),
    })
    per_day = fs * 86400
    for day in range(n_days):
        chunk = frame.iloc[day * per_day:(day + 1) * per_day]
        fname = f"2026-03-{day + 1:02d} 00_00_00.csv"
        chunk.to_csv(acc_dir / fname, index=False)


def test_run_skips_participants_below_min_valid_days(tmp_path, rng):
    """min_valid_days must exclude participants with too little data."""
    study = tmp_path / "study"
    _write_participant(study / "rich" / "accelerometer", "2026-03-01 00:00:00", 5, 2, rng)
    _write_participant(study / "sparse" / "accelerometer", "2026-03-01 00:00:00", 1, 2, rng)
    summary = run(
        str(study),
        str(tmp_path / "out"),
        "America/New_York",
        frequency=Frequency.DAILY,
        interval_seconds=60,
        min_valid_days=3.0
    )
    ids = set(summary["backend_id"])
    assert "rich" in ids
    assert "sparse" not in ids


def test_run_time_bounds_clip_data(tmp_path, rng):
    """time_start/time_end restrict the window; clipping below min_valid_days then skips the
    participant."""
    study = tmp_path / "study"
    _write_participant(study / "rich" / "accelerometer", "2026-03-01 00:00:00", 5, 2, rng)
    unbounded = run(
        str(study), str(tmp_path / "out_all"), "America/New_York",
        frequency=Frequency.DAILY, interval_seconds=60, min_valid_days=3.0,
    )
    assert "rich" in set(unbounded["backend_id"])
    bounded = run(
        str(study),
        str(tmp_path / "out_clip"),
        "America/New_York",
        frequency=Frequency.DAILY,
        interval_seconds=60,
        min_valid_days=3.0,
        time_start="2026-03-01 00_00_00",
        time_end="2026-03-02 12_00_00",
    )
    bounded_ids = set(bounded["backend_id"]) if "backend_id" in bounded.columns else set()
    assert "rich" not in bounded_ids


def test_intradaily_variability_normalizes_over_valid_pairs():
    """IV with gaps must divide each sum by the count of terms that actually
    contributed, not by the raw interval count.

    An alternating 0, 2, 0, 2, ... series has every squared successive
    difference equal to 4 and variance 1, so its IV is exactly 4.0. Dropping
    one 0 and one 2 (equal counts, so the mean stays 1.0) leaves every
    surviving pair at squared-difference 4 and every present point at deviation
    1, so IV must remain exactly 4.0 under valid-pairs normalization. The
    pre-fix formula divides the gap-reduced numerator by the full interval count
    and so returns a biased value below 4.0.
    """
    base = np.array([0.0, 2.0] * 12)          # 24 intervals; even idx 0, odd idx 2
    assert intradaily_variability(base) == pytest.approx(4.0)

    gapped = base.copy()
    gapped[4] = np.nan                         # drops a 0
    gapped[9] = np.nan                         # drops a 2 -> mean unchanged
    assert intradaily_variability(gapped) == pytest.approx(4.0)


def test_average_day_profile_leaves_unobserved_positions_nan():
    """An unobserved clock position must stay NaN in the profile.

    The profile used to be filled in with the mean of the observed positions,
    which made a part of the day nobody was recorded during look as active as
    the rest and pulled L5 and the relative amplitude with it.
    """
    quiet_then_active = np.concatenate([np.full(8, 1.0), np.full(16, 10.0)])
    gappy = quiet_then_active.copy()
    gappy[3] = np.nan

    profile = compute_average_day_profile(gappy, 24)
    assert np.isnan(profile[3])
    assert int(np.isnan(profile).sum()) == 1


def test_l5_m10_skips_windows_that_overlap_a_gap():
    """L5 and M10 come only from windows where every position was observed.

    With hour 3 missing, every 5-hour window that covers it is undefined, so
    the earliest complete quiet window starts at hour 4 and has to reach into
    the active part of the day: four quiet hours and one active hour give
    (1 + 1 + 1 + 1 + 10) / 5 = 2.8. Reporting 1.0 here would mean reading a
    value off a position that was never recorded.
    """
    quiet_then_active = np.concatenate([np.full(8, 1.0), np.full(16, 10.0)])

    complete = compute_l5_m10(quiet_then_active, 24)
    assert complete.l5 == pytest.approx(1.0)
    assert complete.m10 == pytest.approx(10.0)
    assert complete.l5_onset_hours == pytest.approx(0.0)

    gappy = quiet_then_active.copy()
    gappy[3] = np.nan
    windows = compute_l5_m10(gappy, 24)
    assert windows.l5 == pytest.approx(2.8)
    assert windows.l5_onset_hours == pytest.approx(4.0)
    assert windows.m10 == pytest.approx(10.0)


def test_l5_m10_undefined_when_no_window_is_complete():
    """Every metric is NaN when the gaps leave no fully observed window."""
    quiet_then_active = np.concatenate([np.full(8, 1.0), np.full(16, 10.0)])
    scattered_gaps = quiet_then_active.copy()
    scattered_gaps[::5] = np.nan

    windows = compute_l5_m10(scattered_gaps, 24)
    assert np.isnan(windows.l5)
    assert np.isnan(windows.m10)
    assert np.isnan(windows.l5_onset_hours)
    assert np.isnan(windows.m10_onset_hours)


def test_count_valid_days_requires_each_day_to_qualify_on_its_own():
    """Partial days must not add up to whole ones.

    Six days recorded for twelve hours each used to count as three valid days,
    because the old check totalled every observed interval in the recording
    and divided by the length of a day. None of those six days is adequately
    covered, so none of them should count.
    """
    intervals_per_day = INTERVALS_PER_DAY

    six_half_days = np.full(6 * intervals_per_day, np.nan)
    for day in range(6):
        day_start = day * intervals_per_day
        six_half_days[day_start:day_start + intervals_per_day // 2] = 1.0

    # what the old check would have produced, for contrast
    assert np.sum(~np.isnan(six_half_days)) / intervals_per_day == pytest.approx(3.0)
    assert count_valid_days(six_half_days, intervals_per_day, 16.0) == 0

    three_full_days = np.full(3 * intervals_per_day, 1.0)
    assert count_valid_days(three_full_days, intervals_per_day, 16.0) == 3


def test_count_valid_days_is_inclusive_at_the_threshold():
    """A day with exactly the required hours counts; just under does not."""
    intervals_per_day = INTERVALS_PER_DAY
    minutes_in_sixteen_hours = 16 * 60

    exactly_sixteen = np.full(3 * intervals_per_day, np.nan)
    just_under = np.full(3 * intervals_per_day, np.nan)
    for day in range(3):
        day_start = day * intervals_per_day
        exactly_sixteen[day_start:day_start + minutes_in_sixteen_hours] = 1.0
        just_under[day_start:day_start + minutes_in_sixteen_hours - 1] = 1.0

    assert count_valid_days(exactly_sixteen, intervals_per_day, 16.0) == 3
    assert count_valid_days(just_under, intervals_per_day, 16.0) == 0
