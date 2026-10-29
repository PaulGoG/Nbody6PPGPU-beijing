# Notice

This repository is a fork of [nbody6ppgpu/Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing),
the NBODY6++GPU code maintained by Rainer Spurzem's group at ARI Heidelberg and KIAA, Peking University.
The branch `main` carries the corrections I apply for the simulations driven by
[Nbody6Dynamics.jl](https://github.com/PaulGoG/Nbody6Dynamics.jl), one commit per correction on top of the
upstream release v2026.07 (commit 618d7a4). Releases are tagged `v2026.07.N` from N = 6 on; the earlier ones were tagged
`v2026.07-campaign.N` (N = 1–6) on the branch's former name `campaign`, and those tags remain because earlier package
releases pin them (`v2026.07-campaign.6` and `v2026.07.6` hold the same source). The other branches are upstream's as of the fork, unchanged
and not kept in sync. Pull requests upstream follow once the corrections have been validated in production.

## Corrections on this branch

Each is one commit on `main`; the tag that first contains it is given in brackets (`v2026.07-campaign.6` holds
the same source as `v2026.07.6`).

- `escape.F`: a step with more than 1000 escapers no longer ends the run with a format error [v2026.07-campaign.1].
- `energy.F`, `energy_mpi.F`: under a linearised external field (KZ(14) = 1, 2) the potential energy of the
  field enters the total energy and the energy check again [v2026.07-campaign.1].
- `extra_inc/cuda`: helper headers from cuda-samples v13.0, so GPU builds compile with CUDA 13 [v2026.07-campaign.1].
- `configure` outputs untracked, so a configured checkout of a tag is clean [v2026.07-campaign.2].
- `string_left.f`: the digit count of an output interval is bounded; intervals without a short exact binary
  expansion no longer hang the run, and a ten-decimal interval no longer builds an invalid format [v2026.07-campaign.3].
- `mydump.F`: unit-1 COMMON saves (end of run, STOP file, CPU limit, halt) are written whatever `NCOMM`
  [v2026.07-campaign.3]; dumps are written on a fixed unit instead of one numbered by the call counter
  [v2026.07-campaign.3].
- `intgrt.F`: the `COMMON SAVED` line prints each value under its label [v2026.07-campaign.3].
- `configure.ac`, `configure`: `nvcc` is found under `--with-cuda=PREFIX` when it is not on `PATH` [v2026.07-campaign.3].
- `reg.avx.cpp`, `reg.sse.cpp`, `gpunb.velocity.cu`: thread limits raised from 32 and 96 to 512 and 1024
  [v2026.07-campaign.3].
- `xtrnlv.F`: in the Milky-Way potential mode (KZ(14) = 5) an escaper or a mass change no longer carries
  the absolute Galactic potential into the energy bookkeeping; the mode no longer halts at the first
  energy check [v2026.07-campaign.4].
- `fnfw.f`, `fpowercut.f`: jerk of the NFW halo and of the cut-off power-law bulge of MWPotential2014
  [v2026.07-campaign.4].
- `ksint.f`, `unpert.f`, `kstide.f`, `impact.f`, `mdot.F` and the pair tables (`common6.h`, `mydump.F`,
  `ksinit.F`, `ksin2.f`, `remove.f`, `remove_ks.F`, `subint.F`): under KZ(27) > 2 the orbit-averaged
  gravitational-wave decay of a compact KS pair (Peters 1964, DOI 10.1103/PhysRev.136.B1224) is accumulated
  from a per-pair clock written to the dump, so decay below a relative change of 10⁻⁴ is no longer discarded and the KS slow-down is included;
  it is integrated over the accumulated interval in second-order steps and applied at any orbital phase by
  rescaling the KS pericentre and apocentre vectors, so an unperturbed pair no longer gains eccentricity at
  each application; the check is no longer skipped on the unperturbed path without perturbers, and
  unperturbed intervals resolve the decay time. The type −25 is set only while `KSTIDE` applies the decay,
  so a decaying pair is no longer treated as a chaotic binary; the hard-binary exception of `IMPACT` uses its
  intended 50 R☉; the BSE look-up no longer applies the same decay a second time; the per-pair interval is
  broadcast as a real in `mydump.F` (it was sent as an integer) and `mdot.F` guards the decay times against
  a vanishing rate [v2026.07-campaign.5].
- `units.f`, `chain.f`, `qpmod.f`, `rpmin.f`, `rpmax.f`: the speed of light is 299 792.458 km/s instead of
  3×10⁵ km/s [v2026.07-campaign.5].
- `impact.f`, `chterm.f`, `induce.f`, `mwpotinit.F`, `verify.f`: five FORMAT statements carry the comma
  between descriptors that gfortran's run-time format parser requires; with KZ(15) > 0 the first
  binary–binary chain ended every run [v2026.07.6].
- `adjust.F`, `chterm.f`: a run does not end while a triple, quad or chain is active, whose END RUN dump
  could not be restarted; the subsystem is terminated and the end criterion tested again at the next
  adjustment (`END RUN DEFERRED`). A dump skipped during a chain is written at its termination also
  under KZ(2) = 1 [v2026.07.6].
- `chain.f`: steps below 10⁻¹⁴ are counted per chain, and the limit of 10⁶ terminates the chain instead
  of stopping the run without a dump [v2026.07.6].
- `chmod_chain.f`: KZ(26) is restored after the absorption of an inert binary [v2026.07.6].
- `mix.f`: the `&INSSE` choice of `nsflag` is no longer overwritten with 4 (the delayed supernova
  prescription of Fryer et al. 2012, DOI 10.1088/0004-637X/749/1/91) [v2026.07.6].
- `search.f`, `intgrt.F`: a body inside a perturbed KS pair, close to its c.m., whose step has collapsed below
  10⁻⁴ of the c.m. step no longer stalls the block time: SEARCH requests the termination of that pair
  (`KS TERMINATION INSIDE PAIR`), which INTGRT performs at the end of the block [v2026.07.7].
- `mdot.F`: a mass change of a KS component below 10⁻⁸ N-body mass units is applied with the KS orbit and the
  c.m. mass corrected and its energy booked (`HCORR`), as larger ones were; before, the binding energy of a hard
  pair changed by H Δμ without a counterpart at each such step [v2026.07.7].
- `ksint.f`: under KZ(27) > 2 a hyperbolic KS pair of two neutron stars or black holes loses at each
  pericentre the quadrupole energy of one Keplerian passage (Hansen 1972, DOI 10.1103/PhysRevD.5.1021;
  Turner 1977, DOI 10.1086/155501), applied at the pericentre at fixed pericentre distance, and is captured
  when it exceeds the orbital energy (`GR CAPTURE`, `GR FLYBY` lines); before, these pairs received the
  orbit-averaged rate of a bound orbit at every step, and the capture branch used stellar radii, under which
  two black holes were never captured [v2026.07.8].
- `mdot.F`: a mass change below 10⁻⁸ N-body units of a component of a KS pair whose centre of mass is a
  merged hierarchy is applied with the KS orbit corrected and its energy booked, as in other KS pairs
  [v2026.07.8].
- `intgrt.F`: a merged hierarchy terminated by SEARCH (a stalled body inside it) is restored by `RESET`, as at
  the KSINT terminations of such pairs; before, `KSTERM` orphaned the merger entry and the inner ghost [v2026.07.9].
- `offset.f`, `ksint.f`, `unpert.f`: the time offset (KZ(35) > 0) also shifts the reference epochs of the KS
  gravitational-wave decay, which otherwise ran backwards over a negative interval after the first offset; the
  interval is clipped at zero [v2026.07.9].
- `ksint.f`: the KS slow-down factor is applied once in a step whose polynomials the GW decay re-initialised
  (it was applied twice) [v2026.07.9].
- `roche.f`: under KZ(27) > 2 a pair of compact components in mass transfer does not receive the GW decay of the
  look-up in addition to the one integrated in KSINT and UNPERT [v2026.07.9].
- `chain.f`: the tiny-step diagnostic format takes the ten reals of a ten-member chain [v2026.07.9].
- `hrplot.F`: the second component of a KS pair in `bev.82` is evolved with its own main-sequence time (it took
  the first component's) [v2026.07.9].
- `ksint.f`, `unpert.f`, `kstide.f`: under KZ(27) > 2 the Peters decay of a bound compact KS pair is accumulated at
  every check as the decrements of a and e over the interval since the previous check, integrated at the elements of
  that interval (`GWACCUM`), and applied once the sums exceed 10⁻⁴ of a or e; before, the whole sub-threshold
  interval was integrated at the elements of the application, so that an encounter raising e charged the interval
  before it at the new rate (up to 10³ times larger) and brought the merger forward by up to that interval. Three
  per-pair arrays (`GWDA`, `GWDE`, `TGWCHK`) join `COMMON/PAIRS/`, the restart dump and the MPI exchange: a dump
  written by an earlier tag cannot be restarted with this one [v2026.07.10].
- `rmesc_tlist.f`: at the removal of an escaper the last index of every step level is moved into the hole left in
  the block-step list unconditionally; the former guard compared the hole position with the level number and, when
  the two coincided (the escaper among the smallest-step particles, a freshly kicked remnant), left a stale copy of a
  moved index in the list. The duplicated particle was then corrected twice per regular block, once with a zero
  step, its force derivatives became infinite and the run spun in `CORE` without output. `regcor_gpu.F` now stops
  on a zero regular step and `core.f` after 64 cutoff expansions or a non-finite cutoff, so that a corrupted state
  ends the run at once [v2026.07.10].
- `imbhinit.F`, `start.F`, `lagr.f`, `bhplot.f`, `mydump.F`, `ksint.f`, `intgrt.F`, `nbody6.F`, `file_init.F`,
  `common6.h`: the flag route of the massive black holes, `KZ(24) = 3`. After `&INSCALE` (and `&INXTRNL0`) the
  engine reads `NIMBHS DTBH` and a line of up to 16 body names (`dat.10` indices); those bodies keep their `dat.10`
  mass, position and velocity, become stellar type 14, are excluded from the initial single-star mass and from the
  Lagrangian radii (the exclusion of the KS centre of mass of a pair holding one compared a name with an index and
  never acted; it acts now), and each gets its record in `mbh.45`, whose rows end with a new NAME column under both
  routes. The list lives in `COMMON/IMBHS/`, written as an extra record of the restart dump and read back when
  present: dumps of earlier tags restart with the list rebuilt from `NIMBH`, and a black hole missing after a restart
  is dropped with a warning instead of ending the run [v2026.07.11].
- `tstep2.f`, `stepk.f`: a body at exact rest (or with a speed below 10⁻² F/Ḟ in N-body units) gets its first step
  from the force criterion alone; before, its zero velocity time scale gave a zero step, `STEPK` returned a level
  beyond the step table and the run stopped in `SORT_TLIST` with a step of 2×10⁻³¹⁴ [v2026.07.12].
- `kspn.f` (new), `kspn.h` (new), `ksint.f`, `kscorr.f`, `input.F`: post-Newtonian terms of compact KS pairs,
  `KZ(43)`. Bits 1, 2 and 4 add the 1PN, 2PN and 2.5PN relative accelerations (Blanchet 2024, DOI
  10.1007/s41114-024-00050-z; harmonic coordinates, centre-of-mass frame) to the perturbation of every perturbed KS
  step of a bound pair of two neutron stars or black holes whose pericentre compactness m/(c²q) is at least 10⁻⁶;
  bit 8 admits any bound pair (a test mode) and bit 16 writes a `PN KSINT` line every 64 PN steps. The terms are
  added after the relative perturbation γ is formed, so termination, slow-down and perturber selection stay
  Newtonian; such a step has no slow-down, includes the PN strength in the step criterion, and skips the
  orbit-averaged Peters decay of `KZ(27) > 2`, whose osculating elements oscillate under 1PN; unperturbed pairs
  keep that decay and have no periastron advance. The work of the 1PN and 2PN terms is booked as the change of the
  conservative energy (Kidder 1995, DOI 10.1103/PhysRevD.52.821, Eq. 2.7), the 2.5PN work by quadrature, both into
  `ECOLL`, the 2.5PN part also into `EGRAV`. `KZ(43) = 7` with `KZ(27) = 3` is the intended setting (the input
  warns when 2.5PN is selected without `KZ(27) > 2`, which supplies the merger test) [v2026.07.12].
- `mbh.F` (new), `mbh.h` (new), `start.F`, `intgrt.F`, `nbint_cor.f`, `fpoly0.F`, `fpoly1.f`, `fpoly1_ks.f`,
  `kspert.f`, `kspoly.f`, `search.f`, `energy.F`, `energy_mpi.F`, `mydump.F`: a static massive black hole at the
  coordinate origin, `KZ(44) = 1`, read from the namelist `&INMBH MBH, EPSPC, RTMODE, RACCRPC, FACC, BETAF, IPART`
  after the IMBH input (M☉, pc; isolated clusters, `KZ(14) = 0`). Its Plummer-softened force enters the irregular
  force and the initial polynomials of every body and the perturbation of KS pairs as a tidal difference, also in
  the perturbation estimate of the KS candidate search (without it a pair formed deep in the black hole's field
  was integrated as unperturbed); its potential enters the external energy. After every block step a single star
  whose two-body orbit about the black hole has its pericentre within the step and inside r_t/β, or inside
  4 G M/c², is removed at the passage: its energy in the cluster and black-hole potentials goes into `ECOLL`, its
  mass leaves the cluster, it becomes a ghost as in `coal.f`, and a stamped line `TDE` or `TDE PLUNGE` is written.
  r_t is `RACCRPC` (`RTMODE = 2`) or R★ (M/m★)^{1/3} with the stellar radius (`RTMODE = 1`, `KZ(19) ≥ 3`).
  `IPART = 1, 2` replaces the full-disruption threshold by the stripped fraction of Guillochon & Ramirez-Ruiz
  (2013, DOI 10.1088/0004-637X/767/1/25) for γ = 4/3 or 5/3: full disruption above β = 1.85 or 0.9, below that
  down to β = 0.6 or 0.5 the star loses that fraction of its mass with the matching energy and keeps its orbit
  (`TDE PARTIAL`, at most one per pericentre passage). The black-hole mass stays fixed (`FACC` is read and
  checked only). Restart dumps carry the black hole as an optional record [v2026.07.13].
- `kspn.f`, `kspn.h`, `kscorr.f`, `resolv.f`, `energy.F`, `energy_mpi.F`, `remove.f`, `mydump.F`: the work of the
  PN terms of `KZ(43)` is booked as their part of the two-point Hermite rule by which `KSCORR` advances h (start values
  stored per pair, keyed on HDOT and HDOT2, moved with the pair, dumped), the PN part of a prediction that ends a step
  early is booked by `RESOLV`, and `ENERGY` leaves the PN part of the predicted h out of `EBIN`; the energy check then
  closes exactly. The booking of v2026.07.12 (the change of the conservative PN energy and a separate quadrature of
  the 2.5PN work) left about 4×10⁻⁴ of the radiated energy, which halted an embedded run whose binary radiated 10⁷
  times the energy of the cluster [v2026.07.14].
- `fpoly0.F`: the start-up neighbour radius is bisected once radii with too few and too many neighbours are known; a
  star far outside an untruncated Plummer sphere could alternate between the two for ever, and the run never started
  [v2026.07.14].
- `reg.avx.cpp`, `irr.avx.cpp`, `irr.sse.cpp`, `util_gpu.F`, `fpoly0.F`: the massless padding particle of the CPU
  force libraries sat at (255, 255, 255) N-body units with index NTOT + 1; an escaper kept in a long run that came
  within its neighbour radius of that point received it as a neighbour, the list comparison of `REGCOR_GPU` ran past
  its sentinel and beyond LMAX, and the run died of SIGSEGV. The dummy now sits at 10¹⁸ and list indices beyond NTOT
  are dropped [v2026.07.15].
- From upstream `dev`, cherry-picked with their authorship: supernova ghosts removed from neighbour and
  perturber lists in `mdot.F`, ghost guards in `coal.f`, `chterm.f`, `expel.f`, `expel2.f`; `lagr.f`
  cumulative statistics and initialised radii under KZ(8) = 0; the Eggleton IMF exponent in `imf2.f`; HDF5
  handle leaks in `custom_output_facility.F`, with the upstream test `tests/test_lagr_cumulative.py`
  [v2026.07-campaign.3].

## Credit

NBODY6++GPU is the work of Sverre Aarseth (NBODY6), Rainer Spurzem (NBODY6++), Long Wang (the GPU, SIMD and
MPI version), Keigo Nitadori (GPU kernels), the authors of the stellar-evolution routines (Jarrod Hurley,
Onno Pols, Christopher Tout; Sambaran Banerjee; Albrecht Kamlah and collaborators), and the contributors
recorded in the upstream history. Results obtained with this fork cite the upstream code papers:

- Aarseth 1999, PASP 111, 1333, DOI 10.1086/316455
- Spurzem 1999, J. Comput. Appl. Math. 109, 407, DOI 10.1016/S0377-0427(99)00166-1
- Nitadori & Aarseth 2012, MNRAS 424, 545, DOI 10.1111/j.1365-2966.2012.21227.x
- Wang et al. 2015, MNRAS 450, 4070, DOI 10.1093/mnras/stv817
- Kamlah et al. 2022, MNRAS 511, 4060, DOI 10.1093/mnras/stab3748
- Spurzem & Kamlah 2023, Living Rev. Comput. Astrophys. 9, 3, DOI 10.1007/s41115-023-00018-w

## Licensing

The upstream repository carries no licence file. The earlier public repository of the same code,
[nbodyx/Nbody6ppGPU](https://github.com/nbodyx/Nbody6ppGPU), is distributed under the MIT licence
("Copyright (c) 2015 NBODY6++GPU development team", added on 2023-07-28). I treat the inherited code on
those terms, with the attribution above, until the upstream maintainers state terms of their own.

My modifications on the branch `main` (the commits after v2026.07) are released under the MIT licence:

    Copyright (c) 2026 Paul-Adrian Gogîță

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
    associated documentation files (the "Software"), to deal in the Software without restriction,
    including without limitation the rights to use, copy, modify, merge, publish, distribute,
    sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or
    substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
    NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
    OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Third-party files keep their own terms: `extra_inc/cuda/helper_cuda.h` and `helper_string.h` come from
NVIDIA's cuda-samples v13.0 under the BSD-3-Clause licence (`extra_inc/cuda/LICENSE.cuda-samples`);
`extra_inc/cuda/cutil.h`, `extra_inc/nompi/mpif.h` and the other bundled files carry their original notices.
