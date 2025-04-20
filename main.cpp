#include "nets.h"
#include <numeric>

int main() {
    std::string filename;
    std::cout << "Enter CSV filename: ";
    std::cin >> filename;

    double thetaR;
    int thetaK, thetaW, slideSize, dimension, topK;
    std::cout << "Enter thetaR: ";
    std::cin >> thetaR;
    std::cout << "Enter thetaK: ";
    std::cin >> thetaK;
    std::cout << "Enter dimension: ";
    std::cin >> dimension;
    std::cout << "Enter window size: ";
    std::cin >> thetaW;
    std::cout << "Enter slide size: ";
    std::cin >> slideSize;
    std::cout << "Enter top-K dimensions (Dsub): ";
    std::cin >> topK;

    std::vector<DataPoint> allPoints = readCSV(filename, dimension);
    std::vector<int> Dfull(dimension);
    std::iota(Dfull.begin(), Dfull.end(), 0);
    std::vector<int> Dsub = getDsubByVMR(allPoints, topK);

    std::unordered_map<std::string, std::unordered_map<int, DataPoint>> gridSub, gridFull;

    for (size_t i = thetaW + slideSize; i <= allPoints.size(); i += slideSize) {
        size_t startExp = (i >= thetaW + slideSize) ? i - thetaW - slideSize : 0;
        size_t endExp = (i >= thetaW) ? i - thetaW : 0;

        std::vector<DataPoint> expiredSlide(allPoints.begin() + startExp, allPoints.begin() + endExp);
        std::vector<DataPoint> newSlide(allPoints.begin() + i - slideSize, allPoints.begin() + i);

        std::cout << "Processing window [" << i - thetaW << " to " << i - 1 << "]..." << std::endl;

        updateGridDelta(expiredSlide, newSlide, thetaR, Dsub, gridSub);
        updateGridDelta(expiredSlide, newSlide, thetaR, Dfull, gridFull);

        auto cellCatSub = categorizeCellsPoints(gridSub, thetaK);
        std::unordered_map<std::string, CellCategory> finalCat;
        for (const auto& [k, v] : cellCatSub)
            if (v != NON_DETERMINED) finalCat[k] = v;

        std::unordered_map<std::string, std::unordered_map<int, DataPoint>> undecidedFull;
        for (const auto& [k, v] : gridFull)
            if (!finalCat.count(k)) undecidedFull[k] = v;

        auto cellCatFull = categorizeCellsPoints(undecidedFull, thetaK);
        for (const auto& [k, v] : cellCatFull)
            finalCat[k] = v;

        auto outliers = pointLevelDetection(gridFull, finalCat, thetaR, thetaK);

        if (outliers.empty()) std::cout << "No outliers detected.\n";
        else {
            for (const auto& p : outliers) {
                std::cout << "Outlier ID" << p.id << ": (";
                for (size_t j = 0; j < p.coordinates.size(); ++j) {
                    std::cout << p.coordinates[j];
                    if (j != p.coordinates.size() - 1) std::cout << ", ";
                }
                std::cout << ")\n";
            }
        }
    }

    return 0;
}