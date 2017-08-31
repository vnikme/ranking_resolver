#pragma once

#include <vector>
#include <cstddef>


class TRankingResolver {
    public:
        TRankingResolver(size_t height);
        void Add(size_t targetBin, size_t otherBin, double scoreDiffSigma, double weight);
        void MoveTarget(size_t targetBin, size_t otherBin, double scoreDiffSigma, bool otherMoved, double weight);
        void MoveOther(size_t targetBin, size_t otherBin, double scoreDiffSigma, bool targetMoved, double weight);
        std::vector<double> MakeGradient() const;
        std::vector<std::vector<double>> MakeHessian() const;
        std::vector<double> NewtonStep(bool lite) const;
        double Approx(const std::vector<double> &dx) const;

    private:
        std::vector<double> Gradient;               // Size is 2 ^ (height + 1), +1 because of splitting for next level
        std::vector<std::vector<double>> Hessian;   // Size is 2 ^ (height + 1) by 2 ^ (height + 1)
        double SumWeights;
};

