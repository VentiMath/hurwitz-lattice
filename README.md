# Hurwitz lattice verification

Exact-arithmetic Python verification accompanying John Clyde (VentiMath),
*The Hurwitz Lattice and the Two-Slice of Quaternionic Real Bordism*
(manuscript dated October 9, 2026).

Release **v1.0.0** verifies the finite degree-six calculation in Appendix B,
which proves Theorem 1.5: the underlying degree-six Brown–Peterson coefficient
is a sum of six induced sign lattices over the integers localized at 2.
The integral Hurewicz formulas and coefficient certificate are equations (9)
and (10) of that manuscript.

## Run the verification

Use Python 3.12 (the recorded run used Python 3.12.0). From this repository:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verification/degree-six.py
```

On Windows, activate with `.venv\Scripts\activate` instead. The pinned
packages are SymPy 1.14.0 and its dependency mpmath 1.3.0; the verification
itself uses exact integer, rational, and mod-two arithmetic.

A successful run exits with code 0 and ends with:

```text
All checks passed: the degree-six lattice is L_A^6; boundary is injective.
```

A failed assertion exits with an error. Python's `-O` and `-OO` options are
rejected because they would disable assertions. No network access is needed
after installing the dependencies.

## What is checked

- The 24-dimensional Hurewicz basis and its determinant, `2^16`.
- Integral action matrices for the quaternion generators, including
  `i^2 = j^2 = -I` and `ij = -ji`.
- Mod-two fixed-space dimension 6, with zero image in the indecomposable quotient.
- Both integral identities for `Omega S_2` and `Omega S_3`, before reduction mod 2.
- Every entry of the six-by-six coefficient matrix in equation (10),
  and its invertibility over the field with two elements.
- The determinant `729 = 3^6` of the 24 lifted quaternionic orbit vectors,
  a unit after localization at 2.

This is a check of the finite algebraic calculation. It does not compute a
genuine six-slice, determine the degree-six MU coefficient, or certify the
equivariant arguments in Section 6. The manuscript's proof is self-contained.

## Files and recorded outputs

- [`verification/degree-six.py`](verification/degree-six.py): executable verification.
- [`verification/README.md`](verification/README.md): export instructions.
- [`verification/degree-six-check.txt`](verification/degree-six-check.txt): complete recorded output.
- [`verification/degree-six-matrices.txt`](verification/degree-six-matrices.txt): exact matrices.
- [`CITATION.cff`](CITATION.cff): citation metadata for this code release.

To regenerate the recorded outputs after installing dependencies:

```sh
python verification/degree-six.py --matrices verification/degree-six-matrices.txt > verification/degree-six-check.txt
```

Without `--matrices`, the script writes only to standard output. Matrix export
occurs only after every mathematical check passes. GitHub Actions reruns the
verification and compares its output and matrices with the recorded files.

## Citation and provenance

Cite the tagged release
[v1.0.0](https://github.com/VentiMath/hurwitz-lattice/releases/tag/v1.0.0)
when referring to these computations. Later changes will receive new versions.

Author: **John Clyde, VentiMath**. AI tools assisted with implementation and
checking under the author's direction. Responsibility for the mathematical
claims and their presentation rests with the author.

## License

The verification code and documentation are released under the
[MIT License](LICENSE), following the convention of VentiMath's other
verification repository. This repository does not redistribute the manuscript.
