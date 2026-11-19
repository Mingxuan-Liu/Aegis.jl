# Measuring native-block analysis

A reduction reads the requested field and, for a mass-weighted statistic, the selected density. Geometry validation also reads the stored faces and centers of each block. For repeated analyses inside the same open snapshot, construct `CartesianMesh(s)` once. Its validated axes are shared across blocks at the same logical position. Each later selection visits the block metadata and searches the axis segments for the requested centers.

## Choosing a memory budget

For one selected block, let $P$ be the number of cells in its rectangular slab, and $S$ the sum of the stored element sizes of the requested distinct fields. The serial path's largest field payload is $PS$ bytes. The threaded path has one reader and up to $W$ concurrent compute tasks, each assigned at most $K$ blocks. Its maximum in-flight field payload is bounded by $WK P S$ using the largest selected slab for $P$. The reader owns the file handle. Workers receive Julia arrays and return scalar partial sums and moments. They never issue HDF5 calls. `workers` and `chunk_blocks` set $W$ and $K$; the effective number of workers is also bounded by the available default-pool threads and blocks. This path requires no MPI setup.

For example, a full $32^3$ block of one Float32 field occupies 128 KiB. Four workers with two blocks per worker need at most 1 MiB of in-flight field arrays, or 2 MiB for a distinct Float32 density and field. Allocate additional space for Julia, HDF5 buffers, coordinate arrays, descriptors and allocator headroom. These budgets describe live application data. The runtime can retain previously allocated memory until collection.

The snapshot stores the field catalog and three logical indices and a level per block. Geometry construction reads six one-dimensional coordinate vectors per block, retains axis segments proportional to the sum of root-axis lengths, and temporarily uses a coverage bit per block. A reduction visits each selected cell once. A separate reduction call scans the block list again. On chunked or compressed HDF5 datasets, the storage library may fetch and decompress complete underlying storage chunks for a smaller hyperslab. HDF5.jl's MPI backend has separate build requirements. [HDF5.jl parallel I/O](https://juliaio.github.io/HDF5.jl/stable/mpi/) · [Julia threading](https://docs.julialang.org/en/v1/manual/multi-threading/)

## Combining block results

For selected values $f_i$ and nonnegative weights $w_i$, the reported mean is $N/W$, where $N=\sum_i w_i f_i$ and $W=\sum_i w_i$. The code carries compensated sums for these terms, the represented volume and the volume integral. This matters for signed velocity components and for small volumes or masses added to a much larger one. [Neumaier (1974)](https://doi.org/10.1002/zamm.19740540106)

The variance uses a central second moment. Let $M_{2,A}$ and $M_{2,B}$ be the unnormalized central moments of two disjoint selections, with positive weights $W_A,W_B$ and means $\mu_A,\mu_B$. Their merged moment is

$$
M_2 = M_{2,A} + M_{2,B} +
      (\mu_B-\mu_A)^2 \frac{W_A W_B}{W_A+W_B}.
$$

Dividing by the combined weight gives the population variance even when blocks carry very different masses. A zero-weight partial contributes no moment. Signed cancellation, large constant offsets, tiny-mass outliers and full block permutations have dedicated regression cases. Floating-point products and central moments still have finite dynamic range. [Pébay et al. (2016)](https://doi.org/10.1007/s00180-015-0637-z)

## Reproduce the measurement

From the repository root, set `snap` to the target snapshot. Each invocation below runs in a fresh Julia process. The last argument selects an eight-cell slab, a read-only pass over one field, geometry validation, an in-memory arithmetic kernel, or an end-to-end summary. The first two numeric arguments select worker count and blocks per compute task.

```sh
snap=test/test_data/3D_Athena/256_eta5e-3/256_eta5e-3.out2.00010.athdf
julia -t 4 --project=. benchmark/native.jl "$snap" rho 1 1 partial
julia -t 4 --project=. benchmark/native.jl "$snap" rho 1 1 read
julia -t 4 --project=. benchmark/native.jl "$snap" rho 1 1 mesh
julia -t 4 --project=. benchmark/native.jl "$snap" rho 1 1 kernel
julia -t 4 --project=. benchmark/native.jl "$snap" rho 1 1 summary
julia -t 4 --project=. benchmark/native.jl "$snap" rho 4 2 summary
```

The benchmark reports the synthetic warm-up (including compilation), the first access to the real format, and three steady trials. Each I/O trial opens the snapshot and closes it. The kernel trial runs the serial arithmetic on the same preloaded block 32 times. Garbage collection runs before the timer starts. The OS file cache is left intact. Julia's cumulative allocated bytes include arrays reclaimed or reused during the computation. Process high-water RSS includes runtime, compiler, allocations and resident storage-library buffers up to that point. The before/after readings show an incremental high-water change; they cannot isolate an instantaneous peak for a single call.


## Measured 256³ snapshot

The 1,074,162,072-byte local MHD file contains 512 blocks of $32^3$ cells and 16 contiguous Float32 fields. Selecting density reads 64 MiB of field values for the full domain. The slab contains 32 bytes of field values. The test machine was a 32 GB Mac Studio with an Apple M2 Max (eight performance and four efficiency cores), Julia 1.12.0, HDF5.jl 0.17.3 and HDF5 2.1.1. All processes had four Julia compute threads available. The summary worker count below controls how many were used for the reduction.

| Operation | Median steady time (ms) | Cumulative Julia allocations per trial | Process high-water RSS, after warm-up → through last trial (MiB) |
| --- | ---: | ---: | ---: |
| Read a $2 \times 2 \times 2$ slab | 1.0 | 74.8 KiB | 467.6 → 468.7 |
| Read the whole density field blockwise | 20.0 | 65.04 MiB | 475.0 → 558.4 |
| Verify Cartesian geometry | 43.5 | 3.33 MiB | 475.0 → 478.9 |
| Full summary, one worker | 217.1 | 68.55 MiB | 480.2 → 547.8 |
| Full summary, four workers, two blocks per task | 113.2 | 69.27 MiB | 475.5 → 547.6 |

The full-domain density integral was 512.0000000019579 in both summary configurations. On this cached file, four workers reduced the median elapsed time by a factor of 1.92. The synthetic summary warm-up took 1.10 s with one worker and 1.26 s with four. The first real summary took 0.39 s and 0.29 s respectively, including format-specific compilation and I/O. The separate in-memory kernel processed 1,048,576 cell visits in a median 9.11 ms with 4 KiB of Julia allocations for 32 per-block accumulators. Storage throughput, compilation and allocator behavior can differ substantially on other machines and filesystems.
