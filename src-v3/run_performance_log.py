"""
Lightweight per-run wall-clock timings and size snapshots for edna2obis.
Writes a plain-text report under the run output_dir; does not alter pipeline results.
"""

from __future__ import annotations

import os
import time
from contextlib import contextmanager, nullcontext
from typing import Any


def _safe_filename_fragment(name: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(name))[:120]


def _file_size(path: str | None) -> int | None:
    if not path or not os.path.isfile(path):
        return None
    try:
        return os.path.getsize(path)
    except OSError:
        return None


def _collect_config_input_bytes(params: dict[str, Any]) -> tuple[int, int]:
    """Returns (total_bytes, files_counted)."""
    total = 0
    n = 0
    if params.get("use_excel"):
        b = _file_size(params.get("excel_file"))
        if b is not None:
            total += b
            n += 1
    else:
        for key in (
            "sampleMetadata_file",
            "experimentRunMetadata_file",
            "projectMetadata_file",
        ):
            b = _file_size(params.get(key))
            if b is not None:
                total += b
                n += 1
    for _run, paths in (params.get("datafiles") or {}).items():
        if not isinstance(paths, dict):
            continue
        for k in ("taxonomy_file", "occurrence_file", "analysisMetadata_file"):
            b = _file_size(paths.get(k))
            if b is not None:
                total += b
                n += 1
    return total, n


def _collect_table_metrics(raw_data_tables: dict[str, Any]) -> dict[str, Any]:
    total_tax_rows = 0
    total_occ_rows = 0
    max_occ_cols = 0
    n_runs = len(raw_data_tables)
    for _run, tables in raw_data_tables.items():
        if not isinstance(tables, dict):
            continue
        tax = tables.get("taxonomy")
        occ = tables.get("occurrence")
        if tax is not None and hasattr(tax, "__len__"):
            total_tax_rows += len(tax)
        if occ is not None and hasattr(occ, "shape"):
            total_occ_rows += len(occ)
            max_occ_cols = max(max_occ_cols, int(occ.shape[1]))
    return {
        "n_analysis_runs_loaded": n_runs,
        "sum_taxonomy_table_rows": total_tax_rows,
        "sum_occurrence_table_rows": total_occ_rows,
        "max_occurrence_table_columns": max_occ_cols,
    }


class RunPerformanceLog:
    def __init__(self, output_dir: str, run_name: str) -> None:
        self.output_dir = output_dir
        self.run_name = run_name
        self._run_start = time.perf_counter()
        self._step_start = self._run_start
        self.steps: list[tuple[str, float]] = []
        self.metrics: dict[str, Any] = {}
        self.written = False

    def _elapsed_since_run_start(self) -> float:
        return time.perf_counter() - self._run_start

    def mark(self, step_name: str) -> None:
        now = time.perf_counter()
        dt = now - self._step_start
        self.steps.append((step_name, dt))
        self._step_start = now

    @contextmanager
    def step(self, step_name: str):
        t0 = time.perf_counter()
        try:
            yield
        finally:
            self.steps.append((step_name, time.perf_counter() - t0))
            self._step_start = time.perf_counter()

    def set_pre_pipeline_metrics(self, params: dict[str, Any]) -> None:
        b, n = _collect_config_input_bytes(params)
        self.metrics["config_referenced_input_bytes"] = b
        self.metrics["config_referenced_input_files_counted"] = n
        self.metrics["n_datafile_runs_in_config"] = len(params.get("datafiles") or {})

    def set_post_prep_table_metrics(
        self, params: dict[str, Any], raw_data_tables: dict[str, Any]
    ) -> None:
        t = _collect_table_metrics(raw_data_tables)
        self.metrics.update(t)
        b, n = _collect_config_input_bytes(params)
        self.metrics["config_referenced_input_bytes"] = b
        self.metrics["config_referenced_input_files_counted"] = n

    def set_output_file_metrics(self, params: dict[str, Any]) -> None:
        out = params.get("output_dir", "processed-v3/")
        api = str(params.get("taxonomic_api_source", "WoRMS")).lower()
        hc_seconds = params.get("worms_higher_classification_seconds")
        if hc_seconds is not None:
            try:
                self.metrics["worms_higher_classification_seconds"] = float(hc_seconds)
            except (TypeError, ValueError):
                pass
        paths = {
            "occurrence_core": os.path.join(out, f"occurrence_core_{api}.csv"),
            "dna_derived_extension": os.path.join(out, "dna_derived_extension.csv"),
            "eMoF": os.path.join(out, "eMoF.csv"),
            "taxa_xlsx": os.path.join(
                out, f"taxa_assignment_INFO_{params.get('taxonomic_api_source', 'WoRMS')}.xlsx"
            ),
        }
        snap: dict[str, Any] = {}
        for label, p in paths.items():
            b = _file_size(p)
            if b is None:
                snap[f"{label}_bytes"] = None
                continue
            snap[f"{label}_bytes"] = b
            if label in ("occurrence_core", "dna_derived_extension", "eMoF") and p.endswith(
                ".csv"
            ):
                try:
                    with open(p, "rb") as f:
                        snap[f"{label}_line_count_est"] = sum(1 for _ in f)
                except OSError:
                    snap[f"{label}_line_count_est"] = None
        self.metrics["output_files"] = snap

    def write(self, success: bool | None = None, error_message: str | None = None) -> str | None:
        if self.written:
            return None
        os.makedirs(self.output_dir, exist_ok=True)
        frag = _safe_filename_fragment(self.run_name)
        path = os.path.join(self.output_dir, f"edna2obis_run_performance_{frag}.txt")
        total_wall = time.perf_counter() - self._run_start
        lines: list[str] = []
        lines.append("edna2obis run performance log")
        lines.append(f"run_name: {self.run_name}")
        lines.append(f"total_wall_seconds_end_to_end: {total_wall:.4f}")
        if success is not None:
            lines.append(f"pipeline_completed: {success}")
        if error_message:
            lines.append(f"error_message: {error_message}")
        lines.append("")
        lines.append("--- step_timings_seconds (exclusive per step) ---")
        for name, sec in self.steps:
            lines.append(f"{name}\t{sec:.4f}")
        if self.steps:
            lines.append(
                f"sum_step_seconds\t{sum(s for _, s in self.steps):.4f}  (excludes gaps between steps and startup before first timed block)"
            )
        hc_seconds = self.metrics.get("worms_higher_classification_seconds")
        if isinstance(hc_seconds, (int, float)) and hc_seconds >= 0:
            lines.append(
                f"worms_higher_classification\t{hc_seconds:.4f}  (subset of assign_taxonomy when enabled)"
            )
        lines.append("")
        lines.append("--- size_metrics ---")
        for k in sorted(self.metrics.keys()):
            v = self.metrics[k]
            if isinstance(v, dict):
                lines.append(f"{k}:")
                for sk, sv in sorted(v.items()):
                    lines.append(f"  {sk}\t{sv}")
            else:
                lines.append(f"{k}\t{v}")
        lines.append("")
        lines.append("--- rough_scaling_notes ---")
        lines.append(
            "These are single-run linear proxies from THIS run only. Network/API steps (especially "
            "taxonomic assignment) do not scale linearly with table row counts."
        )
        m = self.metrics
        tax_rows = m.get("sum_taxonomy_table_rows") or 0
        occ_rows = m.get("sum_occurrence_table_rows") or 0
        in_bytes = m.get("config_referenced_input_bytes") or 0
        step_map = dict(self.steps)

        def _rate(label: str, numer: float, denom: float, unit: str) -> None:
            if denom and denom > 0 and numer and numer > 0:
                lines.append(
                    f"{label}: {numer:.4f}s / {denom:.4g} {unit} => {numer / denom:.6e} s per 1 {unit}"
                )

        _rate(
            "create_occurrence_core_per_1000_taxonomy_rows",
            step_map.get("create_occurrence_core", 0),
            tax_rows / 1000.0,
            "1k_taxonomy_rows",
        )
        _rate(
            "create_occurrence_core_per_1000_occurrence_rows",
            step_map.get("create_occurrence_core", 0),
            occ_rows / 1000.0,
            "1k_occurrence_rows",
        )
        _rate(
            "assign_taxonomy_per_1000_taxonomy_rows",
            step_map.get("assign_taxonomy", 0),
            tax_rows / 1000.0,
            "1k_taxonomy_rows",
        )
        hc_rate_seconds = m.get("worms_higher_classification_seconds", 0)
        if not isinstance(hc_rate_seconds, (int, float)):
            hc_rate_seconds = 0
        _rate(
            "worms_higher_classification_per_1000_taxonomy_rows",
            hc_rate_seconds,
            tax_rows / 1000.0,
            "1k_taxonomy_rows",
        )
        _rate(
            "load_asv_data_per_gib_input_bytes",
            step_map.get("load_asv_data", 0),
            in_bytes / (1024**3) if in_bytes else 0,
            "GiB_config_paths",
        )
        lines.append("")
        lines.append("--- extrapolation_examples_linear_only ---")
        lines.append(
            "If you assume the same per-unit rates hold, multiply measured seconds by the factor "
            "in parentheses (not valid for API-bound steps under rate limits)."
        )
        for factor in (0.5, 2.0, 5.0):
            occ_s = step_map.get("create_occurrence_core", 0) or 0
            tax_s = step_map.get("assign_taxonomy", 0) or 0
            load_s = step_map.get("load_asv_data", 0) or 0
            lines.append(
                f"factor_{factor}_x_taxonomy_rows => create_occurrence_core_est_s={occ_s * factor:.4f} ; "
                f"assign_taxonomy_est_s={tax_s * factor:.4f} (only if work ~proportional to taxonomy rows)"
            )
            lines.append(
                f"factor_{factor}_x_config_input_bytes => load_asv_data_est_s={load_s * factor:.4f} "
                f"(only if dominated by bytes read from disk)"
            )

        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        self.written = True
        return path


class _NoOpRunPerformanceLog:
    def __init__(self, output_dir: str, run_name: str) -> None:
        pass

    def step(self, step_name: str):
        return nullcontext()

    def set_pre_pipeline_metrics(self, params: dict[str, Any]) -> None:
        pass

    def set_post_prep_table_metrics(
        self, params: dict[str, Any], raw_data_tables: dict[str, Any]
    ) -> None:
        pass

    def set_output_file_metrics(self, params: dict[str, Any]) -> None:
        pass

    def write(self, success: bool | None = None, error_message: str | None = None) -> None:
        return None


def performance_log_for_config(enabled: bool, output_dir: str, run_name: str):
    if enabled:
        return RunPerformanceLog(output_dir, run_name)
    return _NoOpRunPerformanceLog(output_dir, run_name)
