#include "ranking_resolver.h"
#include <random>
#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
#include <iomanip>


struct TItem {
    size_t Bin;
    double Score;
    double Feature;
    bool IsTarget;
    size_t Group;
};

std::vector<TItem> GenerateData(size_t count, size_t n) {
    std::mt19937 rng;
    std::uniform_int_distribution<size_t> idist(0, (1 << n) - 1);
    std::uniform_real_distribution<> rdist(-1.0, 1.0);
    std::vector<TItem> data;
    data.reserve(count * 2);
    for (size_t i = 0; i < count; ++i) {
        data.push_back(TItem { idist(rng), rdist(rng), rdist(rng) + 5.0, true, i });
        data.push_back(TItem { idist(rng), rdist(rng), rdist(rng), false, i });
    }
    return data;
}

int main() {
    size_t n = 2;
    TRankingResolver resolver(n);
    std::vector<TItem> data = GenerateData(50000, n);
    std::vector<size_t> idx(data.size());
    for (size_t i = 0; i < idx.size(); ++i)
        idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&data] (size_t a, size_t b) { return data[a].Feature < data[b].Feature; });
    for (size_t i = 0; i < data.size(); i += 2)
        resolver.Add(data[i].Bin, data[i].Score, data[i + 1].Bin, data[i + 1].Score, 1.0);
    std::cout << std::setprecision(2) << std::fixed; 
    std::set<size_t> movedGroups;
    for (size_t i = 0; i < idx.size(); ++i) {
        if (i % 100 == 0) {
            std::cout << "Iteration " << i / 100 << ":" << std::endl;
            for (double val : resolver.MakeGradient())
                std::cout << val << "\t";
            std::cout << std::endl << std::endl;
            auto step = resolver.NewtonStep();
            std::cout << resolver.Approx(step) << std::endl;
            for (auto val : step)
                std::cout << val << ' ';
            std::cout << std::endl << std::endl << std::endl;
        }
        TItem &obj = data[idx[i]];
        bool wasMoved = (movedGroups.find(obj.Group) != movedGroups.end());
        movedGroups.insert(obj.Group);
        if (obj.IsTarget) {
            TItem &other = data[obj.Group * 2 + 1];
            resolver.MoveTarget(obj.Bin, obj.Score, other.Bin, other.Score, wasMoved, 1.0);
        } else {
            TItem &target = data[obj.Group * 2];
            resolver.MoveOther(target.Bin, target.Score, obj.Bin, obj.Score, wasMoved, 1.0);
        }
    }
    std::cout << "Final:" << std::endl;
    for (double val : resolver.MakeGradient())
        std::cout << val << "\t";
    std::cout << std::endl << std::endl;
    auto step = resolver.NewtonStep();
    std::cout << resolver.Approx(step) << std::endl;
    for (auto val : step)
        std::cout << val << ' ';
    std::cout << std::endl;
    std::cout << std::endl;
    return 0;
}

