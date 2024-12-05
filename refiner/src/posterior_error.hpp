#pragma once
#include <vector>
#include <map>
#include <set>

void calc_rk(const std::vector<std::vector<int>> &cny, 
             const std::vector<std::vector<double>> &coor, 
             const std::vector<std::vector<double>> &displacement,
             const std::vector<std::vector<double>> &material,
             const std::vector<int> &matid_arr,
             std::vector<double> &rk_arr);

void calc_re(const std::vector<std::vector<int>> &cny, 
             const std::vector<std::vector<double>> &coor, 
             const std::vector<std::vector<double>> &displacement,
             const std::vector<std::vector<double>> &material,
             const std::vector<int> &matid_arr,
             const std::map<std::set<int>, std::vector<int>> &face_to_elems,
             std::map<std::set<int>, double> &re_map);

double calc_stress_jump(const std::vector<std::vector<int>> &cny, 
                        const std::vector<std::vector<double>> &coor, 
                        const std::vector<std::vector<double>> &displacement,
                        const std::vector<std::vector<double>> &material,
                        const std::vector<int> &matid_arr,
                        const std::set<int> &face,
                        const std::vector<int> &elems, 
                        const std::map<std::vector<int>, std::vector<std::vector<double>>> &gauss_point_dict);
                        
std::vector<std::vector<double>> calc_dndr(double l0, double l1, double l2, double l3);

std::vector<std::vector<double>> calc_dndx(const std::vector<std::vector<double>> &drdx, 
                                           const std::vector<std::vector<double>> &dndr);

std::map<std::vector<int>, std::vector<std::vector<double>>> set_gauss_point_dict();

double calc_hk(const std::vector<std::vector<double>> &xnode);

std::vector<double> calc_normal_vec(const std::vector<std::vector<double>> &coor, 
                                    const std::set<int> &face);

double calc_area(const std::vector<std::vector<double>> &coor, 
                 const std::set<int> &face);

double calc_he(const std::vector<std::vector<double>> &coor,
               const std::set<int> &face);
               
std::vector<std::vector<std::vector<double>>> calc_dndrdr();

std::vector<std::vector<double>> calc_dxdr(const std::vector<std::vector<double>> &xnode);

std::vector<std::vector<std::vector<double>>> calc_dndxdx(
    const std::vector<std::vector<std::vector<double>>> &dndrdr, 
    const std::vector<std::vector<double>> &drdx);

double calc_volume(std::vector<std::vector<double>> &dxdr);

std::vector<std::vector<double>> calc_drdx(
    const std::vector<std::vector<double>> &dxdr);

void posterior_error_estimation(const std::vector<std::vector<int>> &cny, 
                                const std::vector<std::vector<double>> &coor,
                                const std::vector<std::vector<double>> &displacement,
                                const std::vector<std::vector<double>> &material,
                                const std::vector<int> &matid_arr,
                                const std::map<std::set<int>, std::vector<int>> &face_to_elems,
                                std::vector<double> &eta_arr);