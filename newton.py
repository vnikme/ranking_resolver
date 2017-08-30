#!/usr/bin/python
# coding: utf-8

import random


def mul_ij(X, t, eps, i, j):
    n, m = len(X), len(X[0])
    Y = [[0.0 for v in xrange(m)] for u in xrange(n)]
    for u in xrange(n):
        for v in xrange(m):
            Y[u][v] = t * (X[u][i] - X[u][j]) * (X[i][v] - X[j][v]) + eps * (X[u][i] * X[i][v] + X[u][j] * X[j][v])
    return Y


def add(a, b):
    n, m = len(a), len(a[0])
    c = [[0.0 for j in xrange(m)] for i in xrange(n)]
    for i in xrange(n):
        for j in xrange(m):
            c[i][j] = a[i][j] + b[i][j]
    return c


def sub(a, b):
    n, m = len(a), len(a[0])
    c = [[0.0 for j in xrange(m)] for i in xrange(n)]
    for i in xrange(n):
        for j in xrange(m):
            c[i][j] = a[i][j] - b[i][j]
    return c


def mul(a, b):
    n = len(a)
    c = [[0.0 for j in xrange(n)] for i in xrange(n)]
    for i in xrange(n):
        for j in xrange(n):
            for k in xrange(n):
                c[i][j] += (a[i][k] * b[k][j])
    return c


def mulc(a, b):
    n = len(a)
    c = [[0.0 for j in xrange(n)] for i in xrange(n)]
    for i in xrange(n):
        for j in xrange(n):
            c[i][j] = (a[i][j] * b)
    return c


def Newton(A):
    n = len(A)
    #X = [[1.0 if i == j else 0.0 for j in xrange(n)] for i in xrange(n)]
    #X = [[random.random() * 0.001 if i == j else 0.0 for j in xrange(n)] for i in xrange(n)]
    X = [[0.001 if i == j else 0.0 for j in xrange(n)] for i in xrange(n)]
    for k in xrange(10):
        X = add(X, sub(X, mul(X, mul(A, X))))
    return X


def Print(A):
    n = len(A)
    for i in xrange(n):
        print "\t".join(map(lambda x: "%.4f" % x, A[i]))


def main():
    n, eps = 5, 0.01
    X = [[random.random() * 0.001 if i == j else 0.0 for j in xrange(n)] for i in xrange(n)]
    B = [[0.0 for j in xrange(n)] for i in xrange(n)]
    for k in xrange(5000):
        i, j = 0, 0
        while i == j:
            i, j = random.randint(0, n - 1), random.randint(0, n - 1)
        t = random.random()
        X = add(X, mul_ij(X, t, eps, i, j))
        B[i][i] += (t + eps)
        B[j][j] += (t + eps)
        B[i][j] -= t
        B[j][i] -= t
        #print k
        #Print(mul(X, B))
        #print
    #exit()
    print "Newton:"
    Print(B)
    print
    b1 = Newton(B)
    Print(b1)
    print
    Print(mul(b1, B))
    #for j in xrange(n):
    #    B[j][j] += 0.5


if __name__ == "__main__":
    main()

