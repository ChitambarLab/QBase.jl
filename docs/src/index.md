```@meta
CurrentModule = QBase
```
# QBase.jl

*A base library for quantum information.*

!!! note "Alpha Version"
    Your thoughts and opinions are valued and will help improve this software. To
    provide feedback or make feature requests, please contact [brian.d.doolittle@gmail.com](mailto:brian.d.doolittle@gmail.com).


## Features
  1. Base quantum types with a flexible tolerance for numerical error.
  2. Constructors for quantum state, measurement, and evolution operators.
  3. Methods for evolving and measuring quantum systems.
  4. Methods for calculating information-theoretic quantities.
  5. Mathematics utilities that support quantum mechanics.

## Base Quantum Types

QBase.jl provides a general framework for representing finite quantum
systems and their dynamics.
Quantum mechanics is simply an application of linear algebra where particular
constraints are applied to the vectors and matrices involved in the representation
of quantum systems.
There are three core data structures used to represent quantum systems and their behavior:
  1. **[Bras and Kets](@ref):**  Row and column vectors defined on a complex-valued Hilbert space.
  2. **[Operators](@ref):** Matrices defined on a complex-valued Hilbert space.
  3. **[Probabilities](@ref):** Real-valued vectors and matrices describing the probabilities of events.

Details regarding the definitions and constraints for each of these data structure
are provided in subsequent pages of this documentation.

Ideally, the constraints on quantum objects should be met exactly, however, numerical
errors are inherent to the computations involved.
Therefore, each type has an absolute tolerance parameter `atol` which specifies how
much error is allowed before the quantum object is deemed invalid.
By default, the `atol=1e-7` and is stored in the constant `QBase.ATOL`.
This tolerance is sufficient for most tasks, however, it can easily be relaxed or
tightened as needed.

## Citing

Please see [`CITATION.bib`](https://github.com/ChitambarLab/QBase.jl/blob/master/CITATION.bib)
for reference of this work.

## Licensing

QBase.jl is released under the MIT License.

## Acknowledgements

Development of QBase.jl was made possible by the advisory of Dr. Eric Chitambar
and general support from the Physics Department at the University of Illinois
Urbana-Champaign. Funding was provided by NSF Award 1914440.

## Contents

```@contents
Pages = ["base_types.md", "states.md", "evolution.md", "measurement.md", "informmation.md", "math_utilities.md", "examples.md"]
Depth = 1
```
