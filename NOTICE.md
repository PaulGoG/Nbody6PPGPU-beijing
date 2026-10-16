# Notice

This repository is a fork of [nbody6ppgpu/Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing),
the NBODY6++GPU code maintained by Rainer Spurzem's group at ARI Heidelberg and KIAA, Peking University.
The branch `campaign` carries the corrections I apply for the simulation campaigns driven by
[Nbody6Dynamics.jl](https://github.com/PaulGoG/Nbody6Dynamics.jl), one commit per correction on top of the
upstream release v2026.07 (commit 618d7a4). The other branches are upstream's as of the fork, unchanged
and not kept in sync. Pull requests upstream follow once the corrections have been validated in production.

## Corrections on this branch

Each is one commit on `campaign`; the tag that first contains it is given in brackets.

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
- From upstream `dev`, cherry-picked with their authorship: supernova ghosts removed from neighbour and
  perturber lists in `mdot.F`, ghost guards in `coal.f`, `chterm.f`, `expel.f`, `expel2.f`; `lagr.f`
  cumulative statistics and initialised radii under KZ(8) = 0; the Eggleton IMF exponent in `imf2.f`; HDF5
  handle leaks in `custom_output_facility.F` [v2026.07-campaign.3].

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

My modifications on the branch `campaign` (the commits after v2026.07) are released under the MIT licence:

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
