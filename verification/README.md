# Degree-six verification

`degree-six.py` verifies the degree-six computation in Appendix B of
`quaternionic-slice-counterexample.tex` using exact integer, rational, and
mod-two arithmetic. It requires Python 3 and SymPy; the recorded run used
Python 3.12.0 and SymPy 1.14.0.

If needed, install SymPy for the interpreter you will use with
`python3 -m pip install sympy==1.14.0 mpmath==1.3.0`, or use the repository's
`requirements.txt`.

From the directory containing `verification/`, run:

```sh
python3 verification/degree-six.py
```

The script checks integrality of the quaternionic action matrices, the quaternion
relations, the mod-two fixed space, every entry of the six-by-six coefficient
certificate, the determinant `729 = 3^6` of the 24 orbit vectors, and both
displayed integral identities. A failed check exits with an error. Do not use
Python's `-O` or `-OO` options; the script rejects execution with assertions disabled.

Matrix export is optional. To refresh both recorded outputs:

```sh
python3 verification/degree-six.py --matrices verification/degree-six-matrices.txt > verification/degree-six-check.txt
```

The matrix file is overwritten only after every check passes. Without
`--matrices`, the script writes only to standard output. The mathematical proof
is self-contained; this script checks its finite algebraic calculations and does
not verify or establish a genuine six-slice decomposition.
