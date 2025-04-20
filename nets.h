// nets.h
#ifndef NETS_H
#define NETS_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>

struct DataPoint {
    std::vector<double> coordinates;
    int id;
};

enum CellCategory { INLIER, OUTLIER, NON_DETERMINED };

std::vector<DataPoint> readCSV(const std::string &filename, int expectedDimension);
std::vector<int> getDsubByVMR(const std::vector<DataPoint>& data, int topK);
std::string getCellKey(const DataPoint& point, double cellSize, const std::vector<int>& dims);
void updateGridDelta(
    const std::vector<DataPoint>& expiredSlide,
    const std::vector<DataPoint>& newSlide,
    double cellSize,
    const std::vector<int>& dims,
    std::unordered_map<std::string, std::unordered_map<int, DataPoint>>& grid);
std::unordered_map<std::string, CellCategory>
categorizeCellsPoints(const std::unordered_map<std::string, std::unordered_map<int, DataPoint>>& grid, int thetaK);
std::vector<DataPoint> pointLevelDetection(
    const std::unordered_map<std::string, std::unordered_map<int, DataPoint>>& grid,
    const std::unordered_map<std::string, CellCategory>& cellCategories,
    double thetaR, int thetaK);
double euclideanDistance(const DataPoint& a, const DataPoint& b);

#endif // NETS_H
