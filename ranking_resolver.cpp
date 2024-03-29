#include "ranking_resolver.h"
#include <algorithm>
#include <cmath>
#include <iostream>


TRankingResolver::TRankingResolver(size_t height)
    : Gradient(1 << (height + 1))
    , Hessian(1 << (height + 1), std::vector<double>(1 << (height + 1)))
    , SumWeights(1e-38)
{
}

void TRankingResolver::Add(size_t targetBin, size_t otherBin, double scoreDiffSigma, double weight) {
    SumWeights += weight;
    bool same = (targetBin == otherBin);
    targetBin = targetBin * 2 + 1;
    otherBin = otherBin * 2 + 1;
    double t = scoreDiffSigma * weight, t2 = scoreDiffSigma * (1 - scoreDiffSigma) * weight;
    Gradient[targetBin] -= t;
    Gradient[otherBin] += t;
    size_t a = std::max(targetBin, otherBin), b = std::min(targetBin, otherBin);
    if (same)
        Hessian[a][b] += (2 * t2);
    else
        Hessian[a][b] += t2;
}

void TRankingResolver::MoveTarget(size_t targetBin, size_t otherBin, double scoreDiffSigma, bool otherMoved, double weight) {
    bool same = (targetBin == otherBin);
    targetBin = targetBin * 2;
    otherBin = otherBin * 2 + (otherMoved ? 0 : 1);
    double t = scoreDiffSigma * weight, t2 = scoreDiffSigma * (1 - scoreDiffSigma) * weight;
    Gradient[targetBin] -= t;
    if (same) {
        if (!otherMoved) {
            size_t a = targetBin + 1, b = targetBin;
            Hessian[a][b] += t2;
        } else {
            size_t a = targetBin, b = targetBin + 1;
            Hessian[a][b] += t2;
        }
    } else {
        if (!otherMoved) {
            size_t a = std::max(targetBin, otherBin), b = std::min(targetBin, otherBin);
            Hessian[a][b] += t2;
        } else {
            size_t a = std::max(targetBin + 1, otherBin), b = std::min(targetBin + 1, otherBin);
            Hessian[b][a] += t2;
        }
    }
}

void TRankingResolver::MoveOther(size_t targetBin, size_t otherBin, double scoreDiffSigma, bool targetMoved, double weight) {
    bool same = (targetBin == otherBin);
    targetBin = targetBin * 2 + (targetMoved ? 0 : 1);
    otherBin = otherBin * 2;
    double t = scoreDiffSigma * weight, t2 = scoreDiffSigma * (1 - scoreDiffSigma) * weight;
    Gradient[otherBin] += t;
    if (same) {
        if (!targetMoved) {
            size_t a = otherBin + 1, b = otherBin;
            Hessian[a][b] += t2;
        } else {
            size_t a = otherBin, b = otherBin + 1;
            Hessian[a][b] += t2;
        }
    } else {
        if (!targetMoved) {
            size_t a = std::max(targetBin, otherBin), b = std::min(targetBin, otherBin);
            Hessian[a][b] += t2;
        } else {
            size_t a = std::max(targetBin, otherBin + 1), b = std::min(targetBin, otherBin + 1);
            Hessian[b][a] += t2;
        }
    }
}

namespace {

    using TVector = std::vector<double>;
    using TMatrix = std::vector<TVector>;

    TMatrix Mul(const TMatrix &a, const TMatrix &b) {
        size_t m = a.size(), n = b.size(), v = b.front().size();
        TMatrix c(m, TVector(v));
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < v; ++j)
                for (size_t k = 0; k < n; ++k)
                    c[i][j] += (a[i][k] * b[k][j]);
        return c;
    }

    TVector Mul(const TMatrix &a, const TVector &b) {
        size_t m = a.size(), n = a.front().size();
        TVector c(m);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                c[i] += (a[i][j] * b[j]);
        return c;
    }

    TMatrix Add(const TMatrix &a, const TMatrix &b) {
        size_t m = a.size(), n = a.front().size();
        TMatrix c(m, TVector(n));
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                c[i][j] = a[i][j] + b[i][j];
        return c;
    }

    TMatrix Sub(const TMatrix &a, const TMatrix &b) {
        size_t m = a.size(), n = a.front().size();
        TMatrix c(m, TVector(n));
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                c[i][j] = a[i][j] - b[i][j];
        return c;
    }

    TMatrix Mul(const TMatrix &a, double b) {
        TMatrix c(a);
        for (TVector &row : c)
            for (double &val : row)
                val *= b;
        return c;
    }

    TVector Mul(const TVector &a, double b) {
        TVector c(a);
        for (double &val : c)
            val *= b;
        return c;
    }

    double MaxDiff(const TMatrix &a, const TMatrix &b) {
        double res = 0.0;
        size_t m = a.size(), n = a.front().size();;
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                res = std::max(res, std::abs(a[i][j] - b[i][j]));
        return res;
    }

    double MaxAbs(const TMatrix &a) {
        double res = 0.0;
        size_t m = a.size(), n = a.front().size();
        for (size_t i = 0; i < m; ++i)
            for(size_t j = 0; j < n; ++j)
                res = std::max(res, std::abs(a[i][j]));
        return res;
    }

    bool IsE(const TMatrix &e, double eps) {
        size_t n = e.size();
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                if (std::abs(e[i][j] - (i == j ? 1.0 : 0.0)) > eps)
                    return false;
        return true;
    }

    bool DoInverse(TMatrix a, double add, double seed, double eps, double maxValue, size_t maxSteps, TMatrix &result) {
        size_t n = a.size();
        TMatrix x(n, TVector(n));
        for (size_t i = 0; i < n; ++i) {
            x[i][i] = seed;
            a[i][i] += add;
        }
        for (size_t step = 0; step < maxSteps; ++step) {
            TMatrix e = Mul(a, x);
            if (IsE(e, eps)) {
                result.swap(x);
                return true;
            }
            x = Sub(Mul(x, 2), Mul(x, e));
            if (MaxAbs(x) > maxValue)
                return false;
        }
        return false;
    }

    void Print(const TMatrix &a) {
        for (const TVector &row : a) {
            for (const double &val : row) {
                std::cout << val << '\t';
            }
            std::cout << std::endl;
        }
    }

    void Print(const TVector &a) {
        for (double val : a)
            std::cout << val << '\t';
        std::cout << std::endl;
    }

    void Cholesky(const TMatrix &a, TMatrix &l, TVector &d) {
        size_t n = a.size();
        for (size_t j = 0; j < n; ++j) {
            l[j][j] = 1.0;
            d[j] = a[j][j];
            for (size_t k = 0; k < j; ++k)
                d[j] -= (l[j][k] * l[j][k] * d[k]);
            double id = std::abs(d[j]) > 1e-5 ? 1/d[j] : 0.0;
            for (size_t i = j + 1; i < n; ++i) {
                l[i][j] = a[i][j] * id;
                for (size_t k = 0; k < j; ++k)
                    l[i][j] -= (l[i][k] * l[j][k] * d[k] * id);
                l[j][i] = 0.0;
            }
        }
    }

    TMatrix Diag(const TVector &d) {
        size_t n = d.size();
        TMatrix a(n, TVector(n));
        for (size_t i = 0; i < n; ++i)
            a[i][i] = d[i];
        return a;
    }

    TMatrix Transpose(TMatrix a) {
        size_t n = a.size();
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < i; ++j)
                std::swap(a[i][j], a[j][i]);
        return a;
    }

    void AddRow(TMatrix &a, size_t dest, size_t src, double mul) {
        size_t m = a.front().size();
        for (size_t j = 0; j < m; ++j)
            a[dest][j] += (mul * a[src][j]);
    }

    TMatrix InverseL(TMatrix a) {
        size_t n = a.size();
        for (size_t i = 0; i < n; ++i) {
            a[i].resize(2 * n);
            a[i][i + n] = 1.0;
        }
        for (size_t i = n; i > 0; --i) {
            for (size_t k = i; k < n; ++k)
                AddRow(a, k, i - 1, -a[k][i - 1]);
        }
        for (size_t i = 0; i < n; ++i)
            a[i].erase(a[i].begin(), a[i].begin() + n);
        return a;
    }

    TVector InverseDiag(TVector d) {
        for (double &val : d)
            if (std::abs(val) > 1e-10)
                val = 1 / val;
        return d;
    }

    TMatrix Inverse(const TMatrix &a) {
        size_t n = a.size();
        TMatrix l(n, TVector(n));
        TVector d(n);
        Cholesky(a, l, d);
        l = InverseL(l);
        d = InverseDiag(d);
        return Mul(Transpose(l), Mul(Diag(d), l));
    }

} // namespace

std::vector<double> TRankingResolver::MakeGradient() const {
    TVector gradient(Gradient);
    for (size_t i = 0, n = Gradient.size(); i < n; i += 2)
        gradient[i + 1] -= gradient[i];
    for (double &val : gradient)
        val /= SumWeights;
    return gradient;
}

std::vector<std::vector<double>> TRankingResolver::MakeHessian() const {
    size_t n = Gradient.size();
    TMatrix hessian(n, TVector(n));
    for (size_t i = 0; i < n; i += 2) {
        double v = Hessian[i + 1][i + 1], s = Hessian[i + 1][i], t = Hessian[i][i + 1];
        v -= (s + t);
        double u = (s + t);
        s -= t;
        hessian[i][i] += u;
        hessian[i + 1][i + 1] += v;
        hessian[i + 1][i] -= s;
        hessian[i][i + 1] -= s;
        for (size_t j = 0; j < i; j += 2) {
            double v = Hessian[i + 1][j + 1];
            double s = Hessian[i + 1][j], t = Hessian[i][j + 1];
            double x = Hessian[j][i + 1], y = Hessian[j + 1][i];
            v -= (s + t);
            double u = x + y;
            s -= x;
            t -= y;
            hessian[i][i] += (u + t);
            hessian[i + 1][i + 1] += (s + v);
            hessian[j][j] += (u + s);
            hessian[j + 1][j + 1] += (t + v);
            hessian[i][j] -= u;
            hessian[j][i] -= u;
            hessian[i + 1][j + 1] -= v;
            hessian[j + 1][i + 1] -= v;
            hessian[i + 1][j] -= s;
            hessian[j][i + 1] -= s;
            hessian[i][j + 1] -= t;
            hessian[j + 1][i] -= t;
        }
    }
    for (TVector &row : hessian)
        for (double &val : row)
            val /= SumWeights;
    return hessian;

}

std::vector<double> TRankingResolver::NewtonStep(bool lite) const {
    size_t n = Gradient.size();
    std::vector<double> gradient = MakeGradient();
    std::vector<std::vector<double>> hessian = MakeHessian();
    //Print(Hessian);
    //std::cout << std::endl;
    //Print(hessian);
    //std::cout << std::endl;
    if (!lite) {
        TMatrix inverse = Inverse(hessian);
        //std::cout << "Max diff: " << MaxDiff(hessian, Mul(hessian, Mul(inverse, hessian))) << std::endl;
        return Mul(Mul(inverse, gradient), -1.0);
    }
    for (size_t i = 0; i < n; ++i) {
        if (std::abs(hessian[i][i]) > 1e-10)
            gradient[i] /= -hessian[i][i];
        else
            gradient[i] = 0;
    }
    return gradient;
}

double TRankingResolver::Approx(const std::vector<double> &dx) const {
    TVector gradient = MakeGradient();
    TMatrix hessian = MakeHessian();
    double res = 0.0;
    size_t n = dx.size();
    for (size_t i = 0; i < n; ++i) {
        res += (dx[i] * gradient[i]);
        for (size_t j = 0; j < n; ++j) {
            res += (dx[i] * hessian[i][j] * dx[j] / 2.0);
        }
    }
    return res;
}

