#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <cmath>
#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <map>
#include <functional>
#include <limits>

struct DataPoint {
    std::vector<double> coordinates;
    int id;
};

enum CellCategory { INLIER, OUTLIER, NON_DETERMINED };

std::vector<DataPoint> readCSV(const std::string &filename, int expectedDimension) {
    std::vector<DataPoint> points;
    std::ifstream file(filename);
    std::string line;
    int idCounter = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string token;
        DataPoint dp;
        dp.id = idCounter++;
        while (std::getline(ss, token, ',')) {
            dp.coordinates.push_back(std::stod(token));
        }
        if (dp.coordinates.size() == static_cast<size_t>(expectedDimension)) {
            points.push_back(dp);
        }
    }
    return points;
}

std::vector<int> getDsubByVMR(const std::vector<DataPoint>& data, int topK) {
    int d = data[0].coordinates.size();
    std::vector<std::pair<double, int>> vmrs;
    for (int j = 0; j < d; ++j) {
        std::vector<double> values;
        for (const auto& p : data)
            values.push_back(p.coordinates[j]);
        double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        double var = 0.0;
        for (double v : values) var += (v - mean) * (v - mean);
        var /= values.size();
        double vmr = (mean == 0.0) ? std::numeric_limits<double>::infinity() : var / mean;
        vmrs.emplace_back(vmr, j);
    }
    std::sort(vmrs.begin(), vmrs.end());
    std::vector<int> topDims;
    for (int i = 0; i < topK && i < (int)vmrs.size(); ++i)
        topDims.push_back(vmrs[i].second);
    return topDims;
}

std::string getCellKey(const DataPoint& point, double cellSize, const std::vector<int>& dims) {
    std::ostringstream oss;
    for (size_t i = 0; i < dims.size(); ++i) {
        int index = static_cast<int>(std::floor(point.coordinates[dims[i]] / cellSize));
        oss << index;
        if (i != dims.size() - 1) oss << "_";
    }
    return oss.str();
}

void updateGridDelta(
    const std::vector<DataPoint>& expiredSlide,
    const std::vector<DataPoint>& newSlide,
    double cellSize,
    const std::vector<int>& dims,
    std::unordered_map<std::string, std::unordered_map<int, DataPoint>>& grid)
{
    for (const auto& p : expiredSlide) {
        std::string key = getCellKey(p, cellSize, dims);
        if (grid.count(key) && grid[key].count(p.id)) {
            grid[key].erase(p.id);
            if (grid[key].empty()) grid.erase(key);
        }
    }
    for (const auto& p : newSlide) {
        std::string key = getCellKey(p, cellSize, dims);
        grid[key][p.id] = p;
    }
}

std::unordered_map<std::string, CellCategory>
categorizeCellsPoints(const std::unordered_map<std::string, std::unordered_map<int, DataPoint>>& grid, int thetaK) {
    std::unordered_map<std::string, CellCategory> cellCategories;
    for (const auto& entry : grid) {
        const std::string& cellKey = entry.first;
        int count = entry.second.size();
        int lowerBound = count - 1;
        int upperBound = count - 1;

        std::vector<int> index;
        std::istringstream iss(cellKey);
        std::string token;`
        while (std::getline(iss, token, '_'))
            index.push_back(std::stoi(token));

        std::vector<std::vector<int>> offsets;
        std::vector<int> current(index.size(), 0);
        std::function<void(int)> generateOffsets = [&](int pos) {
            if (pos == index.size()) {
                offsets.push_back(current);
                return;
            }
            for (int i = -1; i <= 1; ++i) {
                current[pos] = index[pos] + i;
                generateOffsets(pos + 1);
            }
        };
        generateOffsets(0);

        for (const auto& offset : offsets) {
            if (offset == index) continue;
            std::ostringstream neighborKey;
            for (size_t i = 0; i < offset.size(); ++i) {
                neighborKey << offset[i];
                if (i != offset.size() - 1) neighborKey << "_";
            }
            if (grid.count(neighborKey.str()))
                upperBound += grid.at(neighborKey.str()).size();
        }

        if (lowerBound >= thetaK)
            cellCategories[cellKey] = INLIER;
        else if (upperBound < thetaK)
            cellCategories[cellKey] = OUTLIER;
        else
            cellCategories[cellKey] = NON_DETERMINED;
    }
    return cellCategories;
}

double euclideanDistance(const DataPoint& a, const DataPoint& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.coordinates.size(); ++i)
        sum += (a.coordinates[i] - b.coordinates[i]) * (a.coordinates[i] - b.coordinates[i]);
    return std::sqrt(sum);
}

std::vector<DataPoint> pointLevelDetection(
    const std::unordered_map<std::string, std::unordered_map<int, DataPoint>>& grid,
    const std::unordered_map<std::string, CellCategory>& cellCategories,
    double thetaR, int thetaK)
{
    std::vector<DataPoint> outliers;
    for (const auto& [cellKey, pointsMap] : grid) {
        if (cellCategories.find(cellKey) == cellCategories.end()) continue;
        if (cellCategories.at(cellKey) != NON_DETERMINED) continue;

        std::vector<int> index;
        std::istringstream iss(cellKey);
        std::string token;
        while (std::getline(iss, token, '_'))
            index.push_back(std::stoi(token));

        std::vector<std::vector<int>> offsets;
        std::vector<int> current(index.size(), 0);
        std::function<void(int)> generateOffsets = [&](int pos) {
            if (pos == index.size()) {
                offsets.push_back(current);
                return;
            }
            for (int i = -1; i <= 1; ++i) {
                current[pos] = index[pos] + i;
                generateOffsets(pos + 1);
            }
        };
        generateOffsets(0);

        std::vector<std::string> neighborKeys;
        for (const auto& offset : offsets) {
            std::ostringstream neighborKey;
            for (size_t i = 0; i < offset.size(); ++i) {
                neighborKey << offset[i];
                if (i != offset.size() - 1) neighborKey << "_";
            }
            neighborKeys.push_back(neighborKey.str());
        }

        for (const auto& [pid, p] : pointsMap) {
            int count = -1;
            for (const auto& nk : neighborKeys) {
                if (!grid.count(nk)) continue;
                for (const auto& [qid, q] : grid.at(nk)) {
                    if (pid == qid) continue;
                    if (euclideanDistance(p, q) <= thetaR)
                        count++;
                    if (count >= thetaK) break;
                }
                if (count >= thetaK) break;
            }
            if (count < thetaK) outliers.push_back(p);
        }
    }
    return outliers;
}

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