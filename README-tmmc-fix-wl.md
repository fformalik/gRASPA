# Branch `tmmc-fix-wl` (fork fformalik/gRASPA)

TMMC changes on top of upstream snurr-group/gRASPA `main` (e6da4e3).
This branch includes everything from `tmmc-dev`.

## New input keywords (component block, after `RunTMMC yes`)

| Keyword | Default | Meaning |
|---|---|---|
| `WriteTMMCEvery N` (global) | 0 (off) | Writes `TMMC/System_X/TMMC_Statistics_<PHASE>_<cycle>.data` every N cycles. The file format is the same as the final `TMMC_Statistics.data`. |
| `RezeroAfterInitialization yes` | no | Clears the C matrix, histogram and bias at the start of equilibration. |
| `UseWLInEquilibration yes` | no | Uses a Wang-Landau bias during equilibration, then switches to the TMMC bias. Requires `UseBiasOnMacrostate yes`. |
| `WLInitialLnF x` | 1.0 | Initial ln f. |
| `WLFlatness x` | 0.8 | The WL histogram counts as flat if min > x * mean over the whole window. |
| `WLFlatnessCheckEvery N` | 10000 | Flatness is checked every N cycles. |
| `WLStagesBeforeSwitch N` | 1 | The bias switches to TMMC after N flat stages, or at the end of equilibration at the latest. |
| `WLCollectAfterStages K` | 0 | The C matrix is collected only after K flat WL stages; at stage K it is zeroed. Must be < `WLStagesBeforeSwitch`. 0 = collect from the start of WL. |

Without the WL keywords, old inputs run as before. The only differences come from the corrected bias (see below).

## Phases with `UseWLInEquilibration yes`

1. **Initialization.** No bias. The C matrix is collected.
2. **Start of equilibration.**
   - With `RezeroAfterInitialization yes`, the C matrix, histogram and TM bias are cleared first.
   - Then WL starts: ln g, the WL bias and the WL histogram are set to 0, and ln f = `WLInitialLnF`.
3. **Equilibration, WL phase.**
   - After every MC step, `ln g[N] += ln f` for the post-move state N, and the bias is −ln g.
   - Every `WLFlatnessCheckEvery` cycles, the WL histogram is checked. If it is flat, ln f is halved and the WL histogram is cleared.
   - After `WLCollectAfterStages` flat stages, the C matrix is zeroed and collected from then on.
   - After `WLStagesBeforeSwitch` flat stages, the bias switches to TMMC.
4. **Production.**
   - If WL is still active, it is switched off.
   - If the WL histogram was never flat, `output.txt` gets a warning with the N range known in the C matrix.
   - From then on the bias is TMBias = −ln Π from the C matrix, updated every `UpdateTMMCEvery` cycles.

The C matrix uses the unbiased Pacc. With `WLCollectAfterStages 0` it is collected from the start of WL. It is never cleared at the switch or at the start of production.

**Use `WLCollectAfterStages` > 0.** An unbiased Pacc is not enough for C to be correct.
- While ln f is large, the WL bias changes every step and pushes the system out of each N before the configuration equilibrates. The transition probabilities in C are then sampled on non-equilibrium configurations.
- In the bulk-water test, C from the first WL stages gave P(2→1) = 0.154 instead of 0.104, and lnΠ that was far too steep.
- The TM bias built from that C pushed every window to its upper edge.
- FEASST WLTM avoids this with `collect_flatness` (18 of 22 stages in its TIP4P tutorial).

The WL log lines in `output.txt` start with `TMMC-WL`.

## Forbidden combinations (error when the input is read)

- `RestartTMMC yes` + `RezeroAfterInitialization yes`: the rezero would erase the restarted statistics.
- `RestartTMMC yes` + `UseWLInEquilibration yes`: not supported yet.
- `UseWLInEquilibration yes` without `UseBiasOnMacrostate yes`.

## Fixes in the TMMC bias (`AdjustTMBias`, `src_clean/data_struct.h`)

- **Range comes from the C matrix.** The link a ↔ a+1 is used only if C(a→a+1) > 0 and C(a+1→a) > 0.
  - The old code took the range from the histogram of move-attempt targets.
  - Its last step was one-sided (P_r = 1), so the bias pushed the system into the first unoccupied state. In the bulk-water W24 window the push was ×1e5, a ratchet of about one state per update.
- **Integration starts at the most visited state** and goes both up and down. It stops at the first unknown link in each direction.
  - Outside the known range, ln Π is held constant, so the bias there is neutral.
- **ln Π is cleared before every update.** Previously a gap in the C matrix brought in the value from the previous update.

## Meaning of the columns in `TMMC_Statistics*.data` (format unchanged)

- `Histogram`: counts of the state AFTER each MC step, one entry per step.
  - Before: the target of every insertion/deletion attempt (accepted or not) and the current state for other moves.
  - Old and new histograms cannot be compared 1:1. The sum is still the number of MC steps since the last rezero.
- `WLBias`, `ln_g`: values from the WL phase, indexed by the relative bin. They are 0 if WL is not used.
  - Before: per-move counts indexed by the absolute N, which read out of bounds when TMMCMin > 0.

## Other fixes

- **Deletion bound check in `TMMC::Update()`.** It compared the relative bin index with the absolute `TMMCMin`. This affected only histogram bookkeeping, not the C matrix.
- **`WLFactor` from `TMMCInitial`** is now read with `stod`. It used to be `stoi`, which turned 0.5 into 0.

## Known issues, not fixed here

- `ApplyTMBiasCBCF` checks `UseWLBias` instead of `UseTMBias`, so with CBCFC + TMMC the TM bias is never applied.
- `WriteTMMCEvery` also writes snapshots during initialization. With `RezeroAfterInitialization` they mean nothing.

## Build and test on Perlmutter

```bash
Cluster-Setup/NERSC_Installation_Nov2025/build_vanilla_perlmutter.sh <build_dir>   # writes <build_dir>/BUILD_INFO with the commit
module load PrgEnv-nvidia nvidia/24.5 cudatoolkit/12.4
nvc++ -std=c++20 -cuda -Isrc_clean tests/test_tmmc_adjust_bias.cpp -lstdc++fs -o test_tmmc_adjust_bias && ./test_tmmc_adjust_bias
```
