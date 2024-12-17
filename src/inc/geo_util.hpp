#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>

#include "reader.hpp"

// utilizing subroutines for tetrahedral computation
double findTetraAspectRatio(std::vector<std::vector<double>> &_verts);

double findCircumRadius(std::vector<std::vector<double>> &_verts);

std::vector<double> findCircumCenter(std::vector<std::vector<double>> &_verts);

double findInRadius(std::vector<std::vector<double>> &_verts);

std::vector<double> findInCenter(std::vector<std::vector<double>> &_verts);

std::vector<double> findformulaN(std::vector<std::vector<double>> &_verts);

double findformulaD(std::vector<std::vector<double>> &_verts);

double findLength(std::vector<double> &_p1, std::vector<double> &_p2);

std::vector<double> normalizeLocTetra(const std::vector<std::vector<double>> &_verts,
                                      const std::vector<double> &_point);

void outputAspectRatio(const std::string &data_dir,
     std::vector<int> &marked_elem, 
     std::vector<std::vector<int>> &cny, 
     std::vector<std::vector<double>> &coor, int nelem_marked);