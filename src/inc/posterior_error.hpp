#pragma once
#include <vector>

void poterior_error_estimation(const std::vector<std::vector<int>> &cny, 
                               const std::vector<std::vector<double>> &coor, 
                               const std::vector<std::vector<double>> &displacement,
                               const std::vector<std::vector<int>> &adj_elem);

double calc_hk(const std::vector<std::vector<double>> &xnode);

std::vector<double> calc_normal_vec(const std::vector<std::vector<double>> &xnode);

std::vector<std::vector<std::vector<double>>> calc_dndrdr();

std::vector<std::vector<double>> calc_dxdr(const std::vector<std::vector<double>> &xnode);

std::vector<std::vector<std::vector<double>>> calc_dndxdx(
    const std::vector<std::vector<std::vector<double>>> &dndrdr, 
    const std::vector<std::vector<double>> &dxdr);

double calc_volume(std::vector<std::vector<double>> &dxdr);

std::vector<std::vector<double>> calc_drdx(
    const std::vector<std::vector<double>> &dxdr);