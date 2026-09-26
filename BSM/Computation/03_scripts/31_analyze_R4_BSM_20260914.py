#!/usr/bin/env python3
"""Deterministic analysis of the newly completed M1-R4 BSM, no simulation.

Requires only already-installed numpy, pandas and Biopython. All original
inputs are read-only; new output and report must not already exist.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import platform
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from Bio import Phylo

JOB = Path(__file__).resolve().parents[1]
SRC = None
N = None
CTX = None
OUT = JOB / "07_tables/bsm_analysis_20260914_v1"
REPORT = JOB / "09_SUMMARIES/BSM_RESULTS_ANALYSIS_20260914.md"
INTERVAL = "empirical 2.5th and 97.5th percentiles across conditional maps; not a confidence interval for the mean"
STATS = ["mean", "sd", "q025", "q975", "mcse"]
SOURCES: dict[str, str] = {}
QA: list[dict] = []
TABLES: dict[str, pd.DataFrame] = {}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bind(path: Path) -> None:
    SOURCES[str(path.relative_to(JOB))] = sha(path)


def read(name: str, folder: Path = None) -> pd.DataFrame:
    p = (SRC if folder is None else folder) / name
    bind(p)
    return pd.read_csv(p, sep="\t", keep_default_na=True)


def check(name: str, passed: bool, detail: str = "") -> None:
    QA.append({"check": name, "pass": bool(passed), "detail": str(detail)})
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def stats(values) -> dict:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if len(a) == 0:
        return {"n_maps": 0, **{k: np.nan for k in STATS}}
    sd = float(np.std(a, ddof=1)) if len(a) > 1 else np.nan
    return dict(n_maps=len(a), mean=float(np.mean(a)), sd=sd,
                q025=float(np.quantile(a, .025, method="linear")),
                q975=float(np.quantile(a, .975, method="linear")),
                mcse=sd / math.sqrt(len(a)))


def grouped(df: pd.DataFrame, keys: list[str], value="count") -> pd.DataFrame:
    rows = []
    for key, g in df.groupby(keys, sort=True, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        rows.append({**dict(zip(keys, key)), **stats(g[value])})
    return pd.DataFrame(rows)


def compare(label: str, actual: pd.DataFrame, original: pd.DataFrame,
            keys: list[str], columns: list[tuple[str, str]] | None = None) -> None:
    if columns is None:
        columns = [(k, k) for k in STATS]
    a = actual.set_index(keys).sort_index()
    b = original.set_index(keys).sort_index()
    check(label + "_keys", a.index.equals(b.index), f"{len(a)} groups")
    for ka, kb in columns:
        av, bv = a[ka].to_numpy(float), b[kb].to_numpy(float)
        delta = float(np.nanmax(np.abs(av - bv)))
        check(label + "_" + ka, np.allclose(av, bv, rtol=0, atol=1e-9, equal_nan=True),
              f"max absolute difference={delta:.17g}; absolute tolerance=1e-9")


def count_grid(label: str, df: pd.DataFrame, keys: list[str], expected: int) -> None:
    check(label + "_unique_grid", len(df) == expected and not df.duplicated(keys).any(), f"rows={len(df)}")
    check(label + "_integer_counts", bool(np.isfinite(df["count"]).all() and
          (df["count"] >= 0).all() and (df["count"] == np.floor(df["count"])).all()))
    check(label + "_map_ids", sorted(df.map_id.unique().tolist()) == list(range(1, N + 1)))


def equality(label: str, a: pd.Series, b: pd.Series) -> None:
    a, b = a.sort_index(), b.sort_index()
    check(label, a.index.equals(b.index) and np.array_equal(a.to_numpy(), b.to_numpy()))


def precision(s: dict) -> float:
    return s["mcse"] / s["mean"] if s["mean"] >= 1 else s["mcse"]


def change(old: float, new: float) -> float:
    return abs(new - old) / max(abs(old), 1e-12) if max(abs(old), abs(new)) >= 1 else abs(new - old)


def main() -> None:
    if OUT.exists() or REPORT.exists():
        raise FileExistsError("Refusing to overwrite existing analysis output/report; use a new version")
    contract_path = JOB / "02_config/conditional_bridge_bsm_R4_20260914.json"
    bind(contract_path)
    contract = json.loads(contract_path.read_text())
    check("frozen_contract", sha(contract_path) == CTX["contract_sha"])
    check("M1_R4_conditional_scope", contract["model"] == "M1" and contract["j"] == 0 and
          contract["scientific_acceptance"] == "NONE" and contract["KKT1"] is False and
          contract["seed_base"] == 202609140)
    areas = read("area_order.tsv", JOB / "02_config").sort_values("analysis_bit_index")
    code_area = dict(zip(areas.analysis_internal_code, areas.abbrev))
    check("frozen_area_order", areas.abbrev.tolist() == contract["area_order"])
    macro = {a: "Old_World" for a in contract["old_world_areas"]}
    macro.update({a: "New_World" for a in contract["new_world_areas"]})
    check("frozen_macro_definition", set(contract["old_world_areas"]) == {"AF","AUS","EPA","IND","MDG","ORI","WPA"} and set(contract["new_world_areas"]) == {"CAM","ENA","SAM","WNA"})
    areas["macroregion"] = areas.abbrev.map(macro)
    TABLES["area_dictionary.tsv"] = areas

    status = read("summary_status.tsv")
    check("completed_checkpoint_stable", len(status) == 1 and status.loc[0, "n_maps"] == N and
          status.loc[0, "status"] == "STABLE" and bool(status.loc[0, "overall_stable"]))
    pm = read("per_map_counts.tsv").sort_values("map_id")
    pp = read("per_map_period_counts.tsv")
    routes = read("dispersal_routes_by_map.tsv")
    rp = read("dispersal_routes_by_map_period.tsv")
    ext = read("extinction_by_area_by_map.tsv")
    ep = read("extinction_by_area_by_map_period.tsv")
    exposure = read("epoch_lineage_exposure.tsv").sort_values("young_ma").reset_index(drop=True)
    exposure["epoch_index"] = np.arange(1, 6)
    exposure["duration_myr"] = exposure.old_ma - exposure.young_ma
    exposure["exposure_fraction"] = exposure.lineage_myr_exposure / exposure.lineage_myr_exposure.sum()
    exposure["denominator_note"] = "total fixed-tree lineage-Myr, not occupied-area or dispersal-opportunity exposure"
    check("contiguous_maps_attempts", pm.map_id.tolist() == list(range(1, N + 1)) and
          pm.attempt_id.tolist() == list(range(1, N + 1)))
    check("fixed_j_and_cladogenesis", (pm.j == 0).all() and (pm.cladogenetic_events == 416).all())
    check("seed_schedule", (pm.map_seed == contract["seed_base"] + pm.attempt_id).all() and
          (pm.source_seed == contract["seed_base"] + 1000000 + pm.attempt_id).all())
    check("positive_map_totals", (pm.d >= 0).all() and (pm.e >= 0).all())
    for name, frame, keys, n in [
        ("period", pp, ["map_id", "epoch", "event"], N*10),
        ("routes", routes, ["map_id", "from_code", "to_code"], N*110),
        ("routes_period", rp, ["map_id", "epoch", "from_code", "to_code"], N*550),
        ("extinction", ext, ["map_id", "area_code"], N*11),
        ("extinction_period", ep, ["map_id", "epoch", "area_code"], N*55)]:
        count_grid(name, frame, keys, n)
    for name, frame in [("routes", routes), ("routes_period", rp)]:
        check(name + "_code_dictionary", frame.from_code.map(code_area).equals(frame.from_area) and
              frame.to_code.map(code_area).equals(frame.to_area) and (frame.from_code != frame.to_code).all())
    for name, frame in [("extinction", ext), ("extinction_period", ep)]:
        check(name + "_code_dictionary", frame.area_code.map(code_area).equals(frame.area))

    # Independent branch-overlap exposure from the frozen Newick, no inference.
    tree_path = JOB / "01_inputs/frozen/tree_scenarioA417.tre"
    bind(tree_path)
    tree = Phylo.read(str(tree_path), "newick")
    depths = tree.depths()
    root_age = max(depths[t] for t in tree.get_terminals())
    raw_exposure = []
    for row in exposure.itertuples():
        v = sum(max(0., min(root_age - depths[parent], row.old_ma) -
                    max(root_age - depths[child], row.young_ma))
                for parent in tree.find_clades(order="preorder") for child in parent.clades)
        raw_exposure.append(v)
    check("independent_tree_exposure", np.allclose(raw_exposure, exposure.lineage_myr_exposure, rtol=0, atol=1e-8),
          f"max abs diff={max(abs(np.array(raw_exposure)-exposure.lineage_myr_exposure)):.17g}; root age={root_age:.9f} Ma")
    TABLES["epoch_exposure.tsv"] = exposure

    equality("d_route_conservation_each_map", routes.groupby("map_id")["count"].sum(), pm.set_index("map_id").d)
    equality("e_area_conservation_each_map", ext.groupby("map_id")["count"].sum(), pm.set_index("map_id").e)
    equality("d_route_conservation_each_map_epoch", rp.groupby(["map_id", "epoch"])["count"].sum(),
             pp[pp.event == "d"].set_index(["map_id", "epoch"])["count"])
    equality("e_area_conservation_each_map_epoch", ep.groupby(["map_id", "epoch"])["count"].sum(),
             pp[pp.event == "e"].set_index(["map_id", "epoch"])["count"])
    equality("routes_epoch_to_whole_tree", rp.groupby(["map_id", "from_code", "to_code"])["count"].sum(),
             routes.set_index(["map_id", "from_code", "to_code"])["count"])
    equality("extinction_epoch_to_whole_tree", ep.groupby(["map_id", "area_code"])["count"].sum(),
             ext.set_index(["map_id", "area_code"])["count"])
    equality("event_epoch_to_whole_tree", pp.groupby(["map_id", "event"])["count"].sum(),
             pm.melt(id_vars="map_id", value_vars=["d", "e"], var_name="event", value_name="count").set_index(["map_id", "event"])["count"])
    nx = read("native_count_crosscheck.tsv")
    check("native_total_crosscheck", nx.totals_match.all() and all((nx[f"independent_{k}"] == nx[f"native_{k}"]).all() for k in ["d", "e", "j"]))

    # Reaggregate the long event table independently, including empty grids.
    ana = read("anagenetic_events_long.tsv")
    check("long_event_ids", not ana.duplicated(["map_id", "event_id"]).any() and len(ana) == int((pm.d + pm.e).sum()))
    check("long_event_types", set(ana.event_type) == {"d", "e"})
    derived_epoch = np.searchsorted(exposure.old_ma.to_numpy(), ana.age_ma.to_numpy(), side="right")
    check("long_event_age_epoch", np.all((derived_epoch >= 0) & (derived_epoch < 5)) and
          np.array_equal(exposure.epoch.to_numpy()[derived_epoch], ana.epoch.to_numpy()) and (ana.age_ma >= 0).all())
    for event, columns, target in [("d", ["map_id", "epoch", "from_code", "to_code"], rp),
                                   ("e", ["map_id", "epoch", "affected_code"], ep.rename(columns={"area_code": "affected_code"}))]:
        observed = ana[ana.event_type == event].groupby(columns).size()
        expected = target.set_index(columns)["count"].sort_index()
        equality("long_" + event + "_reaggregation", observed.reindex(expected.index, fill_value=0), expected)

    totals = pd.DataFrame([{"event": k, **stats(pm[k])} for k in ["d", "e", "j"]])
    totals["interval_definition"] = INTERVAL
    compare("existing_event_totals", totals, read("event_totals_summary.tsv"), ["event"])
    TABLES["event_totals.tsv"] = totals
    TABLES["per_map_event_totals.tsv"] = pm
    pp = pp.drop(columns=["lineage_myr_exposure", "rate_per_lineage_myr"]).merge(exposure, on=["epoch", "young_ma", "old_ma"], validate="many_to_one")
    pp["per_myr"] = pp["count"] / pp.duration_myr
    pp["per_lineage_myr"] = pp["count"] / pp.lineage_myr_exposure
    temporal = grouped(pp, ["epoch", "event"]).rename(columns={k: k + "_count" for k in STATS})
    for value in ["per_myr", "per_lineage_myr"]:
        temp = grouped(pp, ["epoch", "event"], value).drop(columns="n_maps").rename(columns={k: k + "_" + value for k in STATS})
        temporal = temporal.merge(temp, on=["epoch", "event"], validate="one_to_one")
    temporal = temporal.merge(exposure, on="epoch", validate="many_to_one").sort_values(["event", "epoch_index"])
    temporal["interval_definition"] = INTERVAL
    compare("existing_temporal", temporal, read("temporal_event_rate_summary.tsv"), ["epoch", "event"],
            [(k + "_count", k + "_count") for k in STATS] + [(k + "_per_lineage_myr", k + "_rate") for k in STATS])
    TABLES["temporal_events.tsv"] = temporal
    TABLES["per_map_temporal_rates.tsv"] = pp.sort_values(["map_id", "event", "epoch_index"])

    def route_analysis(frame: pd.DataFrame, by_epoch=False) -> tuple[pd.DataFrame, pd.DataFrame]:
        keys = (["epoch"] if by_epoch else []) + ["from_code", "to_code", "from_area", "to_area"]
        z = frame.copy()
        den = pp[pp.event == "d"][["map_id", "epoch", "count"]].rename(columns={"count": "d_total"}) if by_epoch else pm[["map_id", "d"]].rename(columns={"d": "d_total"})
        z = z.merge(den, on=["map_id", "epoch"] if by_epoch else ["map_id"], validate="many_to_one")
        z["share_of_map_d"] = z["count"] / z.d_total.replace(0, np.nan)
        result = grouped(z, keys)
        extra = z.groupby(keys, sort=True).agg(nonzero_fraction=("count", lambda a: float(np.mean(a > 0))),
              total_route_events=("count", "sum"), total_d_events=("d_total", "sum")).reset_index()
        result = result.merge(extra, on=keys, validate="one_to_one")
        result["pooled_share_of_d"] = result.total_route_events / result.total_d_events
        share = grouped(z, keys, "share_of_map_d").rename(columns={"n_maps": "n_maps_with_positive_d", **{k: k + "_share_of_map_d" for k in STATS}})
        result = result.merge(share, on=keys, validate="one_to_one")
        result["route"] = result.from_area + "->" + result.to_area
        result["from_macroregion"] = result.from_area.map(macro)
        result["to_macroregion"] = result.to_area.map(macro)
        result["cross_world"] = result.from_macroregion != result.to_macroregion
        result["interval_definition"] = INTERVAL
        if by_epoch:
            result = result.merge(exposure, on="epoch", validate="many_to_one")
            for k in STATS:
                result[k + "_per_myr"] = result[k] / result.duration_myr
                result[k + "_per_lineage_myr"] = result[k] / result.lineage_myr_exposure
            result = result.sort_values(["epoch_index", "mean", "from_code", "to_code"], ascending=[True, False, True, True])
        else:
            result = result.sort_values(["mean", "from_code", "to_code"], ascending=[False, True, True]).reset_index(drop=True)
            result["rank_by_mean"] = np.arange(1, len(result) + 1)
            result["eligible_nonzero_fraction_ge_0_1"] = result.nonzero_fraction >= .1
            result["plot_top5"] = False
            result.loc[result[result.eligible_nonzero_fraction_ge_0_1].head(5).index, "plot_top5"] = True
        return result, z

    route_summary, route_shares = route_analysis(routes)
    route_period, route_period_shares = route_analysis(rp, True)
    compare("existing_route_summary", route_summary, read("d_route_summary.tsv"), ["from_code", "to_code"])
    TABLES["routes_all.tsv"] = route_summary
    top5 = route_summary[route_summary.plot_top5].copy()
    top5["display_rank"] = np.arange(1, len(top5) + 1)
    TABLES["top5_routes.tsv"] = top5
    TABLES["routes_by_epoch.tsv"] = route_period
    TABLES["per_map_route_shares.tsv"] = route_shares
    top_pairs = pd.MultiIndex.from_frame(top5[["from_code", "to_code"]])
    selected = route_shares[pd.MultiIndex.from_frame(route_shares[["from_code", "to_code"]]).isin(top_pairs)]
    selected_by_map = selected.groupby("map_id").agg(count=("count", "sum"), share_of_map_d=("share_of_map_d", "sum")).reset_index()
    TABLES["top5_combined_share_by_map.tsv"] = selected_by_map
    TABLES["top5_combined_share_summary.tsv"] = pd.DataFrame([{
        "selection_rule": "descending mean; nonzero_fraction >= 0.1; ties by source then destination internal code; first five",
        **stats(selected_by_map.share_of_map_d), "pooled_share_of_d": selected_by_map["count"].sum() / pm.d.sum(),
        "interval_definition": INTERVAL}])

    exsum = grouped(ext, ["area_code", "area"]).sort_values(["mean", "area_code"], ascending=[False, True])
    exsum = exsum.merge(ext.groupby(["area_code", "area"])["count"].apply(lambda a: float(np.mean(a > 0))).rename("nonzero_fraction").reset_index(), on=["area_code", "area"])
    exsum["pooled_share_of_e"] = exsum["mean"] / pm.e.mean()
    exsum["macroregion"] = exsum.area.map(macro)
    exsum["interval_definition"] = INTERVAL
    compare("existing_extinction_summary", exsum, read("e_area_summary.tsv"), ["area_code"])
    TABLES["extinction_by_area.tsv"] = exsum
    experiod = grouped(ep, ["epoch", "area_code", "area"]).merge(exposure, on="epoch", validate="many_to_one")
    for k in STATS:
        experiod[k + "_per_myr"] = experiod[k] / experiod.duration_myr
        experiod[k + "_per_lineage_myr"] = experiod[k] / experiod.lineage_myr_exposure
    TABLES["extinction_by_area_epoch.tsv"] = experiod.sort_values(["epoch_index", "area_code"])

    def macro_analysis(frame: pd.DataFrame, by_epoch=False) -> tuple[pd.DataFrame, pd.DataFrame]:
        z = frame.copy()
        z["from_macroregion"] = z.from_area.map(macro)
        z["to_macroregion"] = z.to_area.map(macro)
        keys = (["epoch"] if by_epoch else []) + ["from_macroregion", "to_macroregion"]
        bymap = z.groupby(["map_id"] + keys)["count"].sum().reset_index()
        den = pp[pp.event == "d"][["map_id", "epoch", "count"]].rename(columns={"count": "d_total"}) if by_epoch else pm[["map_id", "d"]].rename(columns={"d": "d_total"})
        bymap = bymap.merge(den, on=["map_id", "epoch"] if by_epoch else ["map_id"], validate="many_to_one")
        bymap["share_of_map_d"] = bymap["count"] / bymap.d_total.replace(0, np.nan)
        result = grouped(bymap, keys)
        shares = grouped(bymap, keys, "share_of_map_d").rename(columns={"n_maps": "n_maps_with_positive_d", **{k: k + "_share_of_map_d" for k in STATS}})
        result = result.merge(shares, on=keys, validate="one_to_one")
        nz = bymap.groupby(keys)["count"].apply(lambda a: float(np.mean(a > 0))).rename("nonzero_fraction").reset_index()
        result = result.merge(nz, on=keys, validate="one_to_one")
        result["interpretation"] = "imputed-source anagenetic range-expansion events; not independent colonizations or proof of historical direction"
        if by_epoch:
            result = result.merge(exposure, on="epoch", validate="many_to_one")
        return result, bymap

    TABLES["macro_route_summary.tsv"], TABLES["macro_routes_by_map.tsv"] = macro_analysis(routes)
    TABLES["macro_route_summary_by_epoch.tsv"], TABLES["macro_routes_by_map_epoch.tsv"] = macro_analysis(rp, True)
    equality("macro_route_conservation_each_map", TABLES["macro_routes_by_map.tsv"].groupby("map_id")["count"].sum(), pm.set_index("map_id").d)

    # Recompute original stopping statistics from the same map prefixes.
    native_conv = read("bsm_convergence.tsv").set_index("n_maps")
    previous = {}
    checkpoints = []
    for n in range(100, N + 1, 100):
        p = pm[pm.map_id <= n]
        tstats = {k: stats(p[k]) for k in ["d", "e"]}
        rstats = grouped(routes[routes.map_id <= n], ["from_area", "to_area"])
        rstats["route"] = rstats.from_area + "->" + rstats.to_area
        estats = grouped(ext[ext.map_id <= n], ["area"])
        rtop = rstats[rstats["mean"] > 0].sort_values(["mean", "route"], ascending=[False, True]).head(3)
        etop = estats[estats["mean"] > 0].sort_values(["mean", "area"], ascending=[False, True]).head(3)
        temporal_prefix = grouped(pp[pp.map_id <= n], ["epoch", "event"])
        temporal_prefix["precision_metric"] = temporal_prefix.apply(precision, axis=1)
        limiting = temporal_prefix.sort_values("precision_metric", ascending=False).iloc[0]
        max_total = max(precision(v) for v in tstats.values())
        max_temp = float(temporal_prefix.precision_metric.max())
        max_route = max(precision(row) for _, row in rtop.iterrows())
        max_ext = max(precision(row) for _, row in etop.iterrows())
        route_set, ext_set = sorted(rtop.route), sorted(etop.area)
        prevchange = basechange = np.nan
        stable_r_prev = stable_e_prev = stable_r_base = stable_e_base = False
        if n > 100:
            prevchange = max(change(previous[n-100]["totals"][k]["mean"], tstats[k]["mean"]) for k in tstats)
            basechange = max(change(previous[100]["totals"][k]["mean"], tstats[k]["mean"]) for k in tstats)
            stable_r_prev = route_set == previous[n-100]["routes"]
            stable_e_prev = ext_set == previous[n-100]["ext"]
            stable_r_base = route_set == previous[100]["routes"]
            stable_e_base = ext_set == previous[100]["ext"]
        stable = bool(n >= 200 and max(max_total, max_temp, max_route, max_ext, prevchange, basechange) <= .05 and
                      stable_r_prev and stable_e_prev and stable_r_base and stable_e_base)
        row = dict(n_maps=n, max_total_precision_metric=max_total, max_temporal_precision_metric=max_temp,
                   max_top3_route_precision_metric=max_route, max_top3_extinction_precision_metric=max_ext,
                   max_total_change_metric_vs_previous_100=prevchange, max_total_change_metric_vs_100=basechange,
                   overall_stable=stable, top3_routes=";".join(route_set), top3_extinction_areas=";".join(ext_set),
                   limiting_temporal_epoch=limiting.epoch, limiting_temporal_event=limiting.event,
                   limiting_temporal_mean=float(limiting["mean"]), limiting_temporal_mcse=float(limiting.mcse),
                   limiting_metric_kind="relative_mcse" if limiting["mean"] >= 1 else "absolute_mcse",
                   threshold=.05, status="STABLE" if stable else "EXTEND" if n > 100 else "NEEDS_EXTENSION_FOR_CHANGE_CHECK")
        for key in ["max_total_precision_metric", "max_temporal_precision_metric", "max_top3_route_precision_metric", "max_top3_extinction_precision_metric", "max_total_change_metric_vs_previous_100", "max_total_change_metric_vs_100"]:
            check(f"checkpoint_{n}_{key}", np.isclose(row[key], native_conv.loc[n, key], rtol=0, atol=1e-12, equal_nan=True))
        check(f"checkpoint_{n}_stability", stable == bool(native_conv.loc[n, "overall_stable"]))
        check(f"checkpoint_{n}_top_sets", row["top3_routes"] == native_conv.loc[n, "top3_routes"] and row["top3_extinction_areas"] == native_conv.loc[n, "top3_extinction_areas"])
        folder = SRC.parent / f"n_{n:04d}"
        own = read("bsm_convergence.tsv", folder).iloc[-1]
        check(f"checkpoint_{n}_historical_record", int(own.n_maps) == n and bool(own.overall_stable) == stable)
        for k in ["d", "e"]:
            row[k + "_mean"] = tstats[k]["mean"]
            row[k + "_mcse"] = tstats[k]["mcse"]
        checkpoints.append(row)
        previous[n] = dict(totals=tstats, routes=route_set, ext=ext_set)
    TABLES["stability_checkpoints.tsv"] = pd.DataFrame(checkpoints)
    TABLES["qa_checks.tsv"] = pd.DataFrame(QA)

    oldest = exposure.iloc[-1]
    text = f"""# New R4 conditional BSM analysis — 2026-09-14

{N} newly simulated, accepted M1–R4 histories reached the first prespecified passing checkpoint. All {len(QA)} independent aggregation checks passed. No additional model fitting or histories were generated by this analysis.

These results use the current 417-tip tree and revised geography, with the explicitly authorized Pelopidas mathias legacy-code exception. They do not use historical histories, old geography, or R6 outputs. M1 remains the predesignated main time-stratified DEC scenario; R6 is a separate joint geography/range scenario, not a pure range-ceiling test.

## Counts and uncertainty

Expansion d: mean {pm.d.mean():.6f}, SD {pm.d.std(ddof=1):.6f}, empirical 95% history interval {totals.iloc[0].q025:.6f}–{totals.iloc[0].q975:.6f}, MCSE {totals.iloc[0].mcse:.6f}.

Contraction e: mean {pm.e.mean():.6f}, SD {pm.e.std(ddof=1):.6f}, empirical 95% history interval {totals.iloc[1].q025:.6f}–{totals.iloc[1].q975:.6f}, MCSE {totals.iloc[1].mcse:.6f}.

Empirical intervals are 2.5th–97.5th percentiles among conditional histories, not confidence intervals for a mean. MCSE is sample SD / sqrt(n), expressing simulation precision, not biological sampling or model uncertainty. Histories are not independent biological replicates. j is fixed at zero: zero founder events cannot establish biological absence. Local contraction is area loss from a range, not extinction of a lineage.

## Conditional directions and temporal exposure

All 110 ordered inter-area routes and all 11 contraction areas are retained in source tables, including zero cells. Five displayed directions are selected by descending mean count among routes with nonzero frequency >=0.1, ties by internal source then destination code. A route frequency is not a hypothesis-test significance or probability of a uniquely observed colonization. Weighted source-area imputation for widespread ancestors adds model-dependent directional uncertainty. Repeated expansions are not independent colonizations.

Time rates divide counts by fixed-tree lineage-Myr, not by occupied-area opportunity. Age bins are half-open [young, old). The oldest nominal 33.9–45 Ma bin contains {oldest.lineage_myr_exposure:.6f} lineage-Myr, {100*oldest.exposure_fraction:.4f}% of tree exposure. Root age is {root_age:.9f} Ma; the portion older than the root has no lineage exposure. Counts/Myr use nominal stratum duration, explicitly distinct from lineage-time normalization. No fine-scale dispersal pulse or significance claim is made from these fixed bins.

## Numerical and inferential limitations

Original KKT1=FALSE is retained. scientific_acceptance=NONE; technical QA and Monte Carlo stability do not certify global optimization, model robustness, or integration over trees and parameter uncertainty. State probabilities remain the original native postfit precision. Tiny root/tip timing discrepancies from the frozen rounded tree are disclosed, not repaired by altering the tree.

Producer R code decodes each transaction seal on rebuilding and validates map records. The independent Python audit hashes opaque RDS objects and checks every branch audit and flat event record, but does not independently reconstruct every raw RDS event chain. The endpoint-conditioned bridge preserves native effective-generator behavior with a <=1e-12 conditional Poisson-tail bound per branch; native endpoint preparation and effective bridge Q have documented rounding-level differences.

## Reproducibility

The analysis manifest binds all read sources and derived tables by SHA256. Source tables use full precision, sample SD (ddof=1), linear empirical quantiles (R type 7), and explicit missing denominators. Prefix checkpoints are 100-map blocks through {N}, with first-pass 5% MCSE/change and top-three-set rules independently reproduced. No stopping threshold was added after seeing results. No R5 run is part of the user-approved scope.
"""
    OUT.mkdir(parents=True, exist_ok=False)
    for filename, frame in TABLES.items():
        frame.to_csv(OUT / filename, sep="\t", index=False, float_format="%.17g", na_rep="NA")
    REPORT.write_text(text, encoding="utf-8")
    manifest = dict(schema_version="bsm-deterministic-analysis-1.0", created_utc=datetime.now(timezone.utc).isoformat(),
        status="DERIVED_TABLE_QA_PASS_CONDITIONAL_ONLY", scientific_acceptance="NONE", n_maps=N,
        model="M1", maximum_range_size=4, R6="SEPARATE_JOINT_SCENARIO_NOT_INCLUDED", source_KKT1=False,
        conditions="fixed tree, fitted parameters, geographic coding, time-stratified M1 DEC and R4",
        j_fixed=0, source_area_assignment="epoch-specific weighted imputation; not unique observed direction",
        interval_definition=INTERVAL, quantile_method="linear, equivalent to R type 7", sd_ddof=1,
        top5_rule="descending mean; nonzero_fraction>=0.1; ties source then destination internal code; first five",
        pooled_share="sum of route counts across maps / sum of all d counts across maps",
        mean_map_share="mean over maps of route count / d count; zero-event period denominators are NA",
        temporal_denominator="total fixed-tree lineage-Myr, not occupied-area or source-opportunity exposure",
        root_age_ma=root_age, macroregions=macro, source_hashes=SOURCES,
        script={"path": str(Path(__file__).resolve().relative_to(JOB)), "sha256": sha(Path(__file__))},
        versions={"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__},
        checks_passed=len(QA), checks_failed=0,
        outputs={name: {"sha256": sha(OUT / name), "rows": len(frame), "columns": list(frame.columns)} for name, frame in TABLES.items()},
        report={"path": str(REPORT.relative_to(JOB)), "sha256": sha(REPORT)})
    (OUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT), "qa_passed": len(QA), "top5": top5.route.tolist(),
                      "root_age_ma": root_age, "oldest_exposure_fraction": float(oldest.exposure_fraction),
                      "status": manifest["status"]}, ensure_ascii=False, indent=2))


def run(ctx):
    global CTX, SRC, N
    CTX = ctx; SRC = ctx["summary"]; N = ctx["n"]
    main()

if __name__ == "__main__":
    import importlib.util
    spec = importlib.util.spec_from_file_location("r4_bsm_audit", Path(__file__).with_name("30_audit_R4_BSM_20260914.py"))
    audit = importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)
    run(audit.context())
