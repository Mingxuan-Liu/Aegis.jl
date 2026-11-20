# AGENTS.md — Aegis.jl

## Purpose

Aegis.jl is a Julia package for reusable, scalable, physically reliable analysis of Athena++ MHD-PIC simulations. Its scientific components should serve multiple investigations and remain understandable to researchers who did not write them.

These instructions apply to human and agent contributors throughout the repository. Keep them focused on enduring development principles.

## 1. Three requirements for every contribution

Evaluate every design against all three requirements:

- **Generalizability:** expose reusable scientific operations and composable interfaces. Put run parameters, research hypotheses, plotting choices, and case-specific orchestration in configuration, examples, or downstream projects.
- **Performance:** account for computational cost, peak memory, data movement, storage access, and parallel execution from the beginning. The intended workloads can be much larger than available memory.
- **Physical correctness:** make definitions, assumptions, validity limits, and numerical errors explicit, and validate them independently.

Physical correctness constrains optimization. Never silently change a scientific definition, normalization, precision requirement, or physical assumption to obtain speed. Resolve tradeoffs explicitly with evidence. Generalizability means clear extension points and honest capability limits; it does not require implementing every possible format, geometry, or backend in advance.

## 2. Build a scientific library with clear boundaries

- Use Julia as the primary implementation language for the public package and its computational work. Keep auxiliary languages optional and limited to justified interoperability or independent validation.
- Separate storage adapters, physical data representation, numerical operations, workflow composition, and presentation. A file-layout change should not require rewriting a physical diagnostic.
- Implement scientific work in reusable functions. Scripts and notebooks should orchestrate these functions rather than duplicate calculations or become the only implementation of a result.
- Pass dependencies and configuration explicitly. Avoid hidden global state, machine-specific paths, fixed resolutions, fixed snapshot numbers, and assumptions tied to one research project.
- Preserve a small, coherent public interface. Prefer existing, well-supported libraries when appropriate; evaluate their conventions, limitations, and dependencies before adopting them.
- Let demonstrated use cases guide abstraction. Add extensibility where needed without creating speculative frameworks or several unvalidated execution backends.
- Keep the simulation's CPU/GPU implementation distinct from the analysis backend. A common scientific interface may support several readers and execution strategies.

## 3. Treat physical meaning as part of the data

- Define required fields, units, coordinate systems, reference frames, time information, centering, weighting, and applicable physical models at interface boundaries. Preserve enough metadata to interpret derived results later.
- Make code-unit conversions, magnetic-field normalization, equation-of-state assumptions, and transport-coefficient conventions auditable. Do not infer missing physical parameters from familiar variable names alone.
- Respect cell volumes, coordinate metrics, active versus ghost cells, refinement information, and cell-centered versus face-centered fields. Count each represented physical region once in reductions.
- Distinguish operations that do not commute, such as averaging a nonlinear diagnostic and computing that diagnostic from averaged inputs. State which quantity an operation actually calculates.
- Keep volume-, mass-, sample-, and other weighted statistics distinct. Specify the measure for distributions and spectra, including bin widths and linear versus logarithmic coordinates.
- Treat identity and sampling semantics as data requirements where relevant. An array index is not a persistent physical identity; missing samples are not automatically zeros or evidence of a physical event.
- Reject unsupported or ambiguous combinations with actionable messages. Do not silently invent geometry, periodicity, units, missing fields, or complete coverage. Document any scientifically justified default and preserve it in result metadata.
- Distinguish measured simulation outputs from reconstructed estimates. A postprocessed approximation must not be presented as an exact solver quantity.

## 4. Establish the physics before coding it

Before implementing or materially changing a scientific calculation:

1. Write down its defining equation or algorithm, required inputs, output units, normalization, reference frame, and domain of validity.
2. Check those definitions against appropriate primary literature, authoritative documentation, and the actual simulation implementation when output semantics matter. Record stable references and relevant equations or sections in durable documentation or docstrings.
3. Separate physical assumptions, numerical approximations, and empirical choices. Identify which are user-selectable and which restrict applicability.
4. Specify independent expected behavior: an analytic limit, manufactured field, conservation or normalization identity, convergence result, or published benchmark.
5. Resolve consequential ambiguity before presenting the result as validated. Mark unverified behavior explicitly and keep it outside claims of supported scientific functionality.

Existing scripts, agent answers, and agreement with another package are useful evidence, but none replaces checking the physical definition. Do not fabricate references or imply that an uninspected source verifies a formula. Once a definition is documented and verified, reuse that evidence unless the method or assumptions change.

## 5. Design for large data and measurable performance

- State the expected input scale and estimate computational complexity, peak memory, temporary storage, and data movement before choosing an algorithm.
- Support selective access and bounded-memory processing where the operation permits it. Prefer reading only needed fields, regions, blocks, or times and combining partial results. Make expensive full-domain materialization explicit.
- Keep performance-critical Julia functions type-stable and inspect allocations and memory-access patterns. Reuse buffers and avoid unnecessary copies when measurements justify it; do not impose blanket bans on allocation or copying.
- Separate parallel computation from parallel I/O. Verify the concurrency guarantees of the actual storage libraries and file handles. Bound simultaneous reads and avoid loading a full dataset independently into every worker.
- Define ownership of mutable data, work buffers, and partial reductions. Prevent races and oversubscription across tasks, processes, numerical libraries, and devices.
- Maintain a verified serial path or independent small reference calculation. Parallel and accelerated paths must preserve the stated scientific meaning within justified numerical tolerances.
- Treat precision changes, fast-math settings, approximations, and changes in reduction order as numerical decisions requiring validation. Promise bitwise reproducibility only when the algorithm actually provides it.
- Benchmark representative end-to-end workloads as well as important kernels. Separate compilation/startup, steady execution, and I/O effects; report relevant hardware, input sizes, configuration, runtime, memory, and scaling.
- Preserve repeatable benchmark cases. Use stable correctness and resource checks in routine CI; avoid brittle timing gates on uncontrolled machines. A faster language or more threads is not evidence of a speedup.

## 6. Test the scientific contract rigorously

Build validation alongside implementation. For each meaningful scientific feature, choose tests that address its actual failure modes:

- Independent analytic or manufactured cases with known answers, including constant and zero-field limits where meaningful.
- Dimensional, normalization, conservation, symmetry, and limiting-case checks appropriate to the quantity and physical model.
- Geometry, field-centering, array-order, weighting, and mesh-coverage checks for every claimed input class.
- Invalid or incomplete input, missing metadata, non-finite values, empty selections, and degenerate cases with documented behavior.
- Numerical convergence and precision sensitivity for discretized or approximate methods.
- Agreement between serial, parallel, chunked, and accelerated implementations where supported, using justified absolute and relative tolerances.
- Small independent reader or calculation cross-checks, followed by representative real-output integration tests.
- Regressions that reproduce scientifically consequential failures and guard their causes.

Expected values must not simply call the production implementation or repeat its assumptions. Cross-code agreement alone cannot exclude a shared mistake. Explain tolerances using roundoff, truncation, sampling error, or the relevant physical scale; never relax a tolerance merely to make a failure pass.

Keep tests deterministic when possible and record random seeds when sampling is involved. Do not treat every evolving simulation as a closed conservative system: account for boundaries, sources, cooling, forcing, and missing terms before asserting a conservation law. Partition validation into small routine tests and reproducible larger checks so rigor remains practical.

## 7. Require adversarial scientific review

For substantive changes to physics, data interpretation, numerical algorithms, parallel execution, or architecture, use at least one independent adversarial agent when agents are available. Give the reviewer a concrete mandate to find counterexamples, incorrect assumptions, missing tests, and misleading performance or correctness claims. Routine prose or formatting edits do not need this gate.

The reviewer should examine the scientific contract and validation evidence independently of the implementer's explanation. Request separate scientific and performance reviews when the change spans both areas. Agent agreement is not proof: resolve substantive findings through derivations, authoritative sources, independent calculations, or targeted tests.

If independent agent review is unavailable, conduct a documented adversarial review and state the limitation honestly. Do not claim that an independent review occurred. Preserve durable scientific decisions and regression tests, rather than raw agent transcripts or work diaries.

## 8. Write documentation for users and maintainers

The GitHub `README.md` is a user guide to the package. It should explain:

- The scientific purpose, intended users, supported capabilities, and material limitations.
- How to install a working version and reproduce a small, verified example.
- What inputs are needed, which components or scripts to use, how those components relate, and what results they produce.
- The normal path from simulation output through selection and analysis to inspection or export.
- Where to find scientific definitions, further examples, tests, performance guidance, and citation information.

Use plain language and real, checked examples. Explain repository structure in terms of responsibilities and relationships, not a bare directory dump. Keep capabilities and installation instructions synchronized with implemented behavior; label proposals as proposals and do not invent runnable APIs.

Do not turn the README into a report of what an agent changed, removed, fixed, attempted, or left unaddressed. Put release changes in a changelog, task discussions in issue or review records, and scientific derivations in appropriate documentation.

Public scientific functions need docstrings covering their meaning, inputs, outputs, units, assumptions, limitations, and relevant references. Use enough inline comments to explain equations, conventions, indexing, algorithmic choices, and performance-sensitive decisions. Keep comments proportionate and accurate; avoid restating obvious syntax or accumulating redundant narration. Update comments and documentation when behavior changes.

## 9. Keep reproducibility without repository debris

- Record the provenance needed to interpret or reproduce scientific outputs: input identity, simulation and analysis versions where known, parameters, units, selections, weighting, algorithm choices, and execution settings that affect results.
- Preserve original inputs. Make generated outputs and caches identifiable, and document their lifecycle and invalidation rules when introduced.
- Keep dependencies and compatibility requirements explicit. Retain the environments and reference data needed to reproduce supported examples, validation, and benchmarks.
- After implementation and testing, inspect the working tree and remove temporary artifacts created by the work: debug dumps, scratch plots, abandoned variants, copied datasets, accidental logs, caches, and process-only notes.
- Preserve useful test fixtures, regression cases, benchmark scripts, concise scientific decisions, references, and validation evidence. These support future work and are not debris.
- Delete only artifacts attributable to the current work and known to be disposable. Never use broad cleanup to remove unrelated files, user data, or intentional changes. Keep large generated data outside versioned source and use an appropriate ignore policy.

## 10. Completion standard

A contribution is complete when its public behavior and scientific assumptions are documented; independent, relevant tests support its claims; significant performance implications have been measured or clearly bounded; substantive adversarial findings are resolved or explicitly restrict the supported scope; user examples remain accurate; and the repository contains only intentional, useful artifacts.

Report verified outcomes and remaining limitations precisely. Do not describe proposed capabilities, unrun tests, unavailable hardware checks, or unverified scientific conclusions as completed work.

## Authoritative starting points

Use these as entry points, and record the specific source version relevant to an implementation when necessary:

- [Julia performance guidance](https://docs.julialang.org/en/v1/manual/performance-tips/)
- [Julia multithreading and race-safety guidance](https://docs.julialang.org/en/v1/manual/multi-threading/)
- [HDF5.jl parallel I/O requirements](https://juliaio.github.io/HDF5.jl/stable/mpi/)
- [Athena++ reference Python readers](https://github.com/PrincetonUniversity/athena/blob/main/vis/python/athena_read.py)
- [Sun & Bai: the Athena++ MHD-PIC implementation](https://arxiv.org/abs/2304.10568)

These sources do not substitute for the primary literature and exact solver version needed by a particular diagnostic.
