// Unit test for TMMC::AdjustTMBias / Wang-Landau helpers (src_clean/data_struct.h).
// Build and run on Perlmutter (from the repository root):
//   module load PrgEnv-nvidia nvidia/24.5 cudatoolkit/12.4
//   nvc++ -std=c++20 -cuda -Isrc_clean tests/test_tmmc_adjust_bias.cpp -lstdc++fs -o test_tmmc_adjust_bias && ./test_tmmc_adjust_bias
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "data_struct.h"

static TMMC MakeTMMC(size_t Min, size_t Max)
{
  TMMC t;
  t.DoTMMC = true; t.DoUseBias = true; t.UseTMBias = true;
  t.MinMacrostate = Min; t.MaxMacrostate = Max; t.nbinPerMacrostate = 1;
  size_t NBin = Max - Min + 1;
  t.CMatrix.assign(NBin, {0.0, 0.0, 0.0});
  t.WLBias.assign(NBin, 0.0); t.TMBias.assign(NBin, 1.0); t.ln_g.assign(NBin, 0.0);
  t.lnpi.assign(NBin, 0.0); t.forward_lnpi.assign(NBin, 0.0); t.reverse_lnpi.assign(NBin, 0.0);
  t.Histogram.assign(NBin, 0.0); t.WLHistogram.assign(NBin, 0.0);
  return t;
}

// C matrix obeying detailed balance for a target lnpi: per step P(a->a+1) and P(a+1->a)//
static void FillDetailedBalance(TMMC& t, const std::vector<double>& target)
{
  for(size_t a = 0; a + 1 < target.size(); a++)
  {
    double d = target[a+1] - target[a];
    t.CMatrix[a].z   += 0.1 * std::min(1.0, std::exp(d));
    t.CMatrix[a+1].x += 0.1 * std::min(1.0, std::exp(-d));
  }
  for(size_t a = 0; a < target.size(); a++)
  {
    t.CMatrix[a].y = 1.0 - t.CMatrix[a].x - t.CMatrix[a].z;
    t.Histogram[a] = 100.0;
  }
}

static bool Near(double a, double b, double tol = 1e-12) { return std::abs(a - b) < tol; }

static double LnpiDiff(const TMMC& t, size_t a, size_t b) { return t.lnpi[b] - t.lnpi[a]; }

int main()
{
  const std::vector<double> target = {0.0, -1.0, -3.0, -2.5, 0.5, 1.0, -0.5, -4.0, -6.0, -5.0};

  // 1. Complete C matrix: exact recovery of lnpi differences, normalization, bias = -lnpi//
  {
    TMMC t = MakeTMMC(0, 9);
    FillDetailedBalance(t, target);
    t.Histogram[5] = 1000.0; // reference = most visited bin//
    t.AdjustTMBias();
    double norm = 0.0;
    for(size_t a = 0; a < target.size(); a++)
    {
      assert(Near(LnpiDiff(t, 0, a), target[a] - target[0], 1e-10));
      assert(Near(t.TMBias[a], -t.lnpi[a]));
      norm += std::exp(t.lnpi[a]);
    }
    assert(Near(norm, 1.0, 1e-12));
    printf("[ok] complete C matrix: lnpi recovered, sum exp(lnpi) = 1\n");
  }

  // 2. Gap in the middle: C(3->4) = 0 but C(4->3) > 0 -> integration stops at 3//
  //    (reference at 1), states 4..9 get the constant lnpi[3]//
  {
    TMMC t = MakeTMMC(0, 9);
    FillDetailedBalance(t, target);
    t.CMatrix[3].y += t.CMatrix[3].z; t.CMatrix[3].z = 0.0;
    t.Histogram[1] = 1000.0;
    t.AdjustTMBias();
    for(size_t a = 0; a <= 3; a++) assert(Near(LnpiDiff(t, 0, a), target[a] - target[0], 1e-10));
    for(size_t a = 4; a <= 9; a++) assert(Near(t.lnpi[a], t.lnpi[3]));
    printf("[ok] gap C(3->4)=0: known range 0..3, constant beyond\n");

    //Same, reference on the other side of the gap (most visited = 7): known range 4..9//
    TMMC u = MakeTMMC(0, 9);
    FillDetailedBalance(u, target);
    u.CMatrix[3].y += u.CMatrix[3].z; u.CMatrix[3].z = 0.0;
    u.Histogram[7] = 1000.0;
    u.AdjustTMBias();
    for(size_t a = 4; a <= 9; a++) assert(Near(LnpiDiff(u, 4, a), target[a] - target[4], 1e-10));
    for(size_t a = 0; a <= 3; a++) assert(Near(u.lnpi[a], u.lnpi[4]));
    printf("[ok] gap C(3->4)=0, reference above the gap: known range 4..9\n");
  }

  // 3. Missing C(Max->Max-1): the last link is unknown, no one-sided step//
  {
    TMMC t = MakeTMMC(0, 9);
    FillDetailedBalance(t, target);
    t.CMatrix[9].y += t.CMatrix[9].x; t.CMatrix[9].x = 0.0;
    t.Histogram[4] = 1000.0;
    t.AdjustTMBias();
    for(size_t a = 0; a <= 8; a++) assert(Near(LnpiDiff(t, 0, a), target[a] - target[0], 1e-10));
    assert(Near(t.lnpi[9], t.lnpi[8]));
    //Insertion 8 -> 9 is neither pushed nor hindered by the bias//
    double preFactor = 1.0;
    t.ApplyTMBias(preFactor, 8, INSERTION);
    assert(Near(preFactor, 1.0));
    printf("[ok] missing C(Max->Max-1): lnpi[Max] = lnpi[Max-1], neutral bias into Max\n");
  }

  // 4. Missing C(Min->Min+1), window with MinMacrostate > 0 (relative indexing)//
  {
    TMMC t = MakeTMMC(479, 488);
    FillDetailedBalance(t, target);
    t.CMatrix[0].y += t.CMatrix[0].z; t.CMatrix[0].z = 0.0;
    t.Histogram[6] = 1000.0;
    t.AdjustTMBias();
    for(size_t a = 1; a <= 9; a++) assert(Near(LnpiDiff(t, 1, a), target[a] - target[1], 1e-10));
    assert(Near(t.lnpi[0], t.lnpi[1]));
    double preFactor = 1.0;
    t.ApplyTMBias(preFactor, 480, DELETION); //480 -> 479//
    assert(Near(preFactor, 1.0));
    printf("[ok] missing C(Min->Min+1) with TMMCMin=479: lnpi[Min] = lnpi[Min+1], neutral bias\n");
  }

  // 5. No dependence on the previous lnpi: garbage in lnpi/forward/reverse before the update//
  {
    TMMC clean = MakeTMMC(0, 9);
    FillDetailedBalance(clean, target);
    clean.CMatrix[3].y += clean.CMatrix[3].z; clean.CMatrix[3].z = 0.0; //with a gap//
    clean.CMatrix[9].y += clean.CMatrix[9].x; clean.CMatrix[9].x = 0.0;
    clean.Histogram[2] = 1000.0;
    TMMC dirty = clean;
    for(size_t a = 0; a < dirty.lnpi.size(); a++)
    {
      dirty.lnpi[a] = 123.0 + 7.0 * a; dirty.forward_lnpi[a] = -55.0; dirty.reverse_lnpi[a] = 1e30;
    }
    clean.AdjustTMBias();
    dirty.AdjustTMBias();
    //second update on the same data must give the same result as the first//
    TMMC twice = clean; twice.AdjustTMBias();
    for(size_t a = 0; a < clean.lnpi.size(); a++)
    {
      assert(clean.lnpi[a] == dirty.lnpi[a] && clean.TMBias[a] == dirty.TMBias[a]);
      assert(clean.forward_lnpi[a] == dirty.forward_lnpi[a] && clean.reverse_lnpi[a] == dirty.reverse_lnpi[a]);
      assert(clean.lnpi[a] == twice.lnpi[a]);
    }
    printf("[ok] result independent of previous lnpi/forward/reverse values (bitwise)\n");
  }

  // 6. Old upper-edge ratchet case (W24-like): visited 0,1; insertion attempts from 1 recorded//
  //    (C(1->2) > 0) but state 2 never occupied (C(2->1) = 0) -> neutral, not a push//
  {
    TMMC t = MakeTMMC(479, 501);
    t.CMatrix[0] = {0.17, 713.8, 0.0171};
    t.CMatrix[1] = {2.07, 1785.0, 0.0167};
    t.Histogram[0] = 586; t.Histogram[1] = 1545;
    t.AdjustTMBias();
    double up = 1.0;   t.ApplyTMBias(up, 480, INSERTION);  //480 -> 481//
    double down = 1.0; t.ApplyTMBias(down, 480, DELETION); //480 -> 479 (known link)//
    assert(Near(up, 1.0));
    double expected_down = std::exp(LnpiDiff(t, 0, 1)); //Pi(480)/Pi(479)//
    assert(Near(down, expected_down, 1e-12));
    printf("[ok] W24-like edge: bias 480->481 = %.3g (old code: ~1.1e5), bias 480->479 = %.3g\n", up, down);
  }

  // 7. Nothing sampled yet: bias left unchanged//
  {
    TMMC t = MakeTMMC(0, 9);
    t.AdjustTMBias();
    for(double b : t.TMBias) assert(b == 1.0);
    assert(t.TMUpdateTimes == 0);
    printf("[ok] empty histogram: bias unchanged\n");
  }

  // 8. Wang-Landau helpers: relative indexing (no out-of-bounds), flatness, switch//
  {
    TMMC t = MakeTMMC(479, 481);
    t.UseWLInEquilibration = true; t.WLInitialLnF = 1.0; t.WLFlatness = 0.8;
    t.ln_g[1] = 99.0; t.WLHistogram[2] = 5.0; //leftovers must be cleared by StartWL//
    t.StartWL();
    assert(t.UseWLBias && !t.UseTMBias && t.WLFactor == 1.0 && t.ln_g[1] == 0.0 && t.WLHistogram[2] == 0.0);
    t.RecordState(480); t.RecordState(480); t.RecordState(479);
    assert(t.ln_g[1] == 2.0 && t.WLBias[1] == -2.0 && t.ln_g[0] == 1.0 && t.WLHistogram[1] == 2.0);
    t.RecordState(500); //outside the window: ignored//
    assert(t.Histogram[0] + t.Histogram[1] + t.Histogram[2] == 3.0);
    assert(!t.CheckWLFlatness()); //state 481 never visited//
    t.RecordState(481); t.RecordState(481); t.RecordState(479);
    assert(t.CheckWLFlatness()); //2,2,2 -> flat//
    assert(t.WLFactor == 0.5 && t.WLStagesDone == 1 && t.WLHistogram[0] == 0.0);
    double preFactor = 1.0; //WL bias: insertion 479 -> 480 uses exp(ln g(479) - ln g(480))//
    t.ApplyWLBias(preFactor, 479, INSERTION);
    assert(Near(preFactor, std::exp(t.ln_g[0] - t.ln_g[1])));
    t.SwitchWLToTM();
    assert(!t.UseWLBias && t.UseTMBias);
    printf("[ok] WL: StartWL reset, relative ln_g, flatness -> ln f/2, switch to TM\n");
  }

  // 9. WLCollectAfterStages: no C collection until K flat stages, then C reset and collected//
  {
    TMMC t = MakeTMMC(0, 2);
    t.UseWLInEquilibration = true; t.WLCollectAfterStages = 2; t.WLStagesBeforeSwitch = 3;
    t.CMatrix[1] = {5.0, 5.0, 5.0}; //data from initialization, must be discarded//
    t.StartWL();
    assert(!t.CollectCMatrix);
    t.Update(0.5, 1, INSERTION);
    assert(t.CMatrix[1].z == 5.0); //not collected//
    for(int stage = 0; stage < 2; stage++)
    {
      assert(!t.StartCollectionIfDue());
      for(size_t N = 0; N <= 2; N++) t.RecordState(N);
      assert(t.CheckWLFlatness());
    }
    assert(t.StartCollectionIfDue() && t.CollectCMatrix);
    assert(t.CMatrix[1].x == 0.0 && t.CMatrix[1].y == 0.0 && t.CMatrix[1].z == 0.0);
    t.Update(0.5, 1, INSERTION);
    assert(t.CMatrix[1].z == 0.5 && t.CMatrix[1].y == 0.5);
    assert(t.SwitchWLToTM()); //was collecting//

    //forced switch before K stages: C zeroed, collection starts, returns false//
    TMMC u = MakeTMMC(0, 2);
    u.UseWLInEquilibration = true; u.WLCollectAfterStages = 2; u.WLStagesBeforeSwitch = 3;
    u.CMatrix[0] = {1.0, 1.0, 1.0};
    u.StartWL();
    assert(!u.SwitchWLToTM());
    assert(u.CollectCMatrix && u.UseTMBias && !u.UseWLBias && u.CMatrix[0].z == 0.0);
    printf("[ok] WLCollectAfterStages: no C before K stages, reset at K, forced switch zeroes C\n");
  }

  printf("all TMMC tests passed\n");
  return 0;
}
