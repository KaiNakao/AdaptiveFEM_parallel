#pragma once
#include <iostream>
#include <string>
#include <vector>

void write_eta(const int &myid, const std::vector<double> &eta_arr);

void write_marked_elem(const int &myid, const std::vector<int> &marked_elem);

void write_displacement_obs(
    const std::vector<std::vector<double>> &displacement_obs, const int &iload);