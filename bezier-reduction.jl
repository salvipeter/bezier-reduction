module Reduction

# Based on:
#   H. Sunwoo: Matrix representation for multi-degree reduction of Bezier curves.
#     CAGD 22(3):261-273, 2005. https://doi.org/10.1016/j.cagd.2004.12.002

using LinearAlgebra

"""
    reduction_matrix(n, m, r, s)

Computes a matrix Q for reducing a Bézier curve of degree `n` to degree `m`,
while retaining the first `r` (`s`) derivatives at the 0 (1) parameter.
It is assumed that `r + s < m < n - 1`.

The new control points can be computed as `Q * P`,
where `P` is the vector of control points.

All computations are exact, using rational numbers with big integers.
"""
function reduction_matrix(n, m, r, s)
    (n, m, r, s) = BigInt.([n, m, r, s])
    N, M = n - (r + s + 2), m - (r + s + 2)
    L = computeL(M, 2r + 2, 2s + 2)
    E = computeE(N, 2r + 2, 2s + 2)
    A = computeA(n, m, max(r, s))
    D = diagm([binomial(n, r + 1 + i) // binomial(N, i) for i in 0:N]) # Lemma 4

    # Lemma 3 (Eq. 41)
    C = zeros(Rational{BigInt}, N + 1, n + 1)
    for j in r+1:n+r-m, k in 0:r
        d = -binomial(n, k) // binomial(n, j)
        C[j-r,k+1] = d * sum(i -> binomial(n - m, j - i) * A[i-k+1], max(k, j - (n - m)):r)
    end
    for j in r+1:n-s-1
        C[j-r,j+1] = 1
    end
    for j in m-s:n-s-1, k in n-s:n
        d = -binomial(n, k) // binomial(n, j)
        C[j-r,k+1] = d * sum(i -> binomial(n - m, j - m + i) * A[i+k-n+1], max(n - k, m - j):s)
    end

    # Lemma 7 (Eq. 45)
    D1 = zeros(Rational{BigInt}, m + 1, M + 1)
    for i in r+1:m-s-1
        D1[i+1,i-r] = binomial(M, i - r - 1) // binomial(m, i)
    end

    # Theorem 2 (Eq. 36)
    Q1 = zeros(Rational{BigInt}, m + 1, n + 1)
    for j in 0:r, k in 0:j
        Q1[j+1,k+1] = binomial(n, k) // binomial(m, j) * A[j-k+1]
    end
    for j in m-s:m, k in j+(n-m):n
        Q1[j+1,k+1] = binomial(n, k) // binomial(m, j) * A[k-j-(n-m)+1]
    end

    # Theorem 3 (Eq. 46)
    Q2 = D1 * L' * Matrix(I, M + 1, N + 1) * E' * D * C

    # Theorem 4 (Eq. 47)
    Q1 + Q2
end

# Lemma 1 (Eq. 18)
function computeL(n, r, s)
    L = zeros(Rational{BigInt}, n + 1, n + 1)
    for k in 0:n, j in 0:n
        for i in max(0, j + k - n):min(j, k)
            b = binomial(k, i) * binomial(n - k, j - i) // binomial(n, j) # Eq. 10
            sign = (k + i) % 2 == 0 ? 1 : -1
            d = binomial(k + r, i) * binomial(k + s, k - i) // binomial(k, i)
            L[k+1,j+1] += sign * b * d
        end
    end
    L
end

# Lemma 2 (Eq. 19)
function computeE(n, r, s)
    E = zeros(Rational{BigInt}, n + 1, n + 1)
    for k in 0:n, j in 0:n
        for i in 0:j
            sign = (j + i) % 2 == 0 ? 1 : -1
            d = binomial(j + r, i) * binomial(j + s, j - i) // binomial(n + r + s + j, k + s + i)
            E[k+1,j+1] += sign * d
        end
        d = binomial(j + r + s, r) * binomial(n, k) // binomial(j + r, r)
        E[k+1,j+1] *= (2j + r + s + 1) // (n + r + s + j + 1) * d
    end
    E
end

# Eq. 32
function computeA(n, m, l)
    A = zeros(Rational{BigInt}, l + 1)
    A[1] = 1
    for k in 1:l
        A[k+1] = -sum(i -> binomial(n - m, k - i) * A[i+1], 0:k-1)
    end
    A
end

end
