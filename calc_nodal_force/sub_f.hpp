#pragma once
extern "C" {
void read_nload_(int &nload);
void read_nobs_xnode_(int &nobs, double &xnode);
void read_centroid_(double &xc, double &mvec);
void calc_fvec_(double &fvec, double &xnode, double &xc, double &mvec);
}