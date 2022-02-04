# bezier-reduction
Optimal degree reduction matrix for Bézier curves, based on a [paper by H. Sunwoo](https://doi.org/10.1016/j.cagd.2004.12.002).

Both in Julia and in plain C. The computation is exact (using bignum rationals), which is then converted to doubles.
