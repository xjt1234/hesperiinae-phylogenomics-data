# R4 M1 bounded conditional BSM: numerical and scope note

Date: 2026-09-14. User authorization: “R4按照计划继续分析”.

All M0/M1/M2 R4 fit and postfit stages have completed. The predetermined main model remains M1 (conservative, five epochs, DEC, d/e free, j=0, w=1). This BSM continuation does not fit or optimize any model and does not operate on R6.

M1 source fit: `04_runs/M1_R4_final_v1/fit_result.rds`, SHA-256 `bb855fe8bfc2919d5eed94224b177463d517dc78fbe958cc34faf44fcca1bff4`.
M1 fixed-parameter postfit: `04_runs/M1_R4_postfit_v1/postfit_recalculated_ancestral_states.rds`, SHA-256 `566b86c007d123be0fccf1064abf076ad26e6c082ceac8133cbed447ef52e08d`.

The new fit has finite log likelihood -1124.49742109455; convergence code 0, KKT1 FALSE and KKT2 TRUE. Its three screens closely reproduce the likelihood, and the independent postfit has zero likelihood difference and normalized ancestral probabilities. These observations support bounded conditional follow-up, but are not a certificate of strict KKT satisfaction or a global optimum. No claims from old geography-specific fixed-point or double-precision tests are transferred to this run.

Input identity: same 417-tip tree; geographic coding amended except that Pelopidas_mathias uses its explicitly authorized legacy code. The comparison with R6 changes both this species coding and maximum range size; it is a joint scenario comparison, not a pure range-cap test. Historical R4/R6 caches and histories are forbidden as inputs.

The event-counting, transaction, seed, warning, and stopping routines are imported from a byte-pinned vendored source. The private branch sampler reuses the proven uniformization implementation. It simulates the native waiting rates and normalized destination law (Qeff), while the native endpoint preparation continues to use Qraw. They differ by small documented single-precision rounding residuals; the helper rejects discrepancies exceeding 1e-7 absolute or 4*2^-23 relative, and controls the endpoint-relative Poisson tail at 1e-12. No installed namespace is modified.

Evidence is layered. The newly executed generic selftest is synthetic algorithm evidence only, explicitly not a current-fit or whole-tree validation. The new runner's preflight verifies this actual fit, postfit, package bodies, 562 states, 5 strata, 1891-row stratified arrays, 833-node probabilities and exact input paths. It then builds a fresh native cache from the new postfit. Each pilot/production map must have 1473 audited branch calls, finite tail bounds, exactly 416 unique true cladogenetic nodes, one-area DEC gains/losses, zero j events, exact independent/native event totals, conserved area/epoch totals, zero warnings and no force-fit history. Pilot passage is not production completion.

Sampling is finite: one-map pilot; first 100 histories; then complete 100-history checkpoints through 500 maximum, with 1000 global attempts and stop after 10 consecutive failed attempts. Stop at the first checkpoint >=200 passing the unchanged 5% MCSE, cumulative-change and top-three set rules. Report INCOMPLETE if the cap is reached without stability; do not create an endless retry chain.

The per-history 2.5% and 97.5% quantiles describe variability among conditional histories, not confidence intervals for their means. MCSE estimates Monte Carlo precision, not phylogenetic/model/geographic uncertainty. Dispersal source areas are imputed with epoch-specific weights using an independent recorded seed; route counts are not independently demonstrated colonizations. With j fixed at zero, absence of sampled j events is a model invariant, not an empirical rejection of founder events.

Scientific acceptance remains NONE / conditional review required. The final main/supplement figures and Windows handoff must retain these limitations.
