# How to Find Ternary LWE Keys Using Locality Sensitive Hashing

Scripts accompanying the paper 
**[How to Find Ternary LWE Keys Using Locality Sensitive Hashing](https://eprint.iacr.org/2021/1255.pdf)**

## Requirements

* [Sage](https://www.sagemath.org/)

## Examples

To call May'21 estimates (implemented only for depth 3) on NTRU parameters, run

```bash
sage may21.sage
```

For other parameters, modify the script manually.

To call KM21 estimames (implemented up tp depth 4), run

```bash
sage no_guessing.sage
```

To compare against Drop&Solve attack from the [LatticeEstimator](https://github.com/malb/lattice-estimator) run

```bash
sage lattice_attack_complexity.sage
```
