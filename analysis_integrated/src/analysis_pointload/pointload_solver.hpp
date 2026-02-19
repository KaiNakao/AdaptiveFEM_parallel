#pragma once
#include <vector>
#include <iostream>
#include <mpi.h>

// interface to fortran programs
#ifdef __cplusplus
extern "C" {
#endif
    void initialize_(int *im, int *ncpu, int *np, int *kd, int *nfa, int *n2dvisl, int *n2dvisa, 
                     int *n_size, int *ne_size, int *n, int *ne, int *ib, int *nab, int *nad, int *nadto, 
                     int *n4_size, int *ne4_size, int *n4, int *ne4, int *ib4, int *nab4, int *nad4, int *nadto4, 
                     long *bufls, long *bufds, long *fidc, long *fidm, int *ierr);
    void main_compute_ni_(int *im, int *ncpu, int *np, int *kd, int *nfa, int *n2dvisl,
            int *n_size, int *ne_size, int *n, int *ne, int *ib, int *nab, int *nad, int *nadto, long *fidm,
            int *n4_size, int *ne4_size, int *n4, int *ne4, int *ib4, int *nab4, int *nad4, int *nadto4, long *fidc,
            long *bufds, long *bufls,
            double *plane, double *ds, int *ni_size, int *nei_size, int *ni, int *nei,
            int *ni4_size, int *nei4_size, int *ni4, int *nei4);
    void main_solver_(int* im,int* ncpu,int* np,int* kd,int* nfa,double* plane,double* ds,int* n2dvisl,
                      int* n_size,int* ne_size,int* ni_size,int* nei_size,int* n,int* ne,int* ni,int* nei,
                      int* ib,int* nab,int* nad,int* nadto,long* fidm,
                      int* n4_size,int* ne4_size,int* ni4_size,int* nei4_size,int* n4,int* ne4,int* ni4,int* nei4,
                      int* ib4,int* nab4,int* nad4,int* nadto4,long* fidc,
                      long* bufds,long* bufls,
                      double* coor,double* coori,double* dsi,double* uval,double* bcarray,int* cny,int* num,
                      int* cnyi,int* numi,int* nfla,double* duplixyz,
                      int* mpiadi,int* mpiadlisti,int* mpli,int* nadtoi,
                      double* coori4,double* dsi4,double* uval4,double* bcarray4,int* cny4,int* num4,
                      int* cnyi4,int* numi4,int* nfla4,double* duplixyz4,
                      int* mpiadi4,int* mpiadlisti4,int* mpli4,int* nadtoi4,
                      int* nird,int* niwd,int* pentawptrd,int* pentarptrd,int* pentawindd,int* pentarindd,
                      int* nir,int* niw,int* pentawptr,int* pentarptr,int* pentawind,int* pentarind,
                      int* nir4,int* niw4,int* pentawptr4,int* pentarptr4,int* pentawind4,int* pentarind4,
                      double* kcrsvald,double* kcrsval,double* kcrsval4,
                      int* map4to10,
                      int* cnysep_write,int* cnysep_write4,
                      int* nsep_ne_ptr,int* numsep,int* nsep_n_ptr,int* nsep_n_ptr4,int* nsep,
                      int* nsepnmap,int* nsep4,int* nsepnmap4,
                      int* ncolor,int* color_ind,
                      int* cnysep_readm,int* cnysep_readm4,
                      double* younglst,double* sig,double* up,double* rv,double* mpibuf,double* uvalcoe,
                      int* log_innerCG,int* log_innerCG_per_outite,
                      int* use_crs_d, int* overlap_comm_d,
                      int* overlap_comm_s, int* persistent_comm_s);
    void pointload_(int* im, int* n, double* coor, int* ne, int* cny, int* ni, double* rv,
                    double* load_arr, int* nload, int* iload, int* load_elem_arr,
                    int* surf_nelem, int* surf_cny);
    void sync_parallel_block_d_(int* n_size, int* n, double* b, int* nad, int* nadto, int* mpiad, int* mpiadlist, int* mpl, int* ierr, int* im, double* mpibuf);
    void ipcg_(int* im, int* np, int* kd, double* younglst, double* sig, double* mpibuf,
               int* n, int* ne, double* coor, int* cny, int* num, double* uval, double* duplixyz, double* bcarray,
               int* ni, int* nei, int* cnyi, double* dsi, int* numi, int* nfla,
               int* nad, int* nadto, int* mpiad, int* mpiadlist, int* mpl,
               int* nird, int* pentarptrd, int* pentarindd,
               int* niwd, int* pentawptrd, int* pentawindd, double* pentacrsvald,
               int* nir, int* pentarptr, int* pentarind,
               int* niw, int* pentawptr, int* pentawind, double* pentacrsval,
               int* nsep, int* numsep, int* nsep_n_ptr, int* nsep_ne_ptr, int* nsepnmap,
               int* cnysep_read, int* cnysep_write,
               int* ncolor, int* color_ind,
               int* n4, int* ne4, double* coor4, int* cny4, int* num4, double* uval4, double* duplixyz4, double* bcarray4,
               int* ni4, int* nei4, int* cnyi4, double* dsi4, int* numi4, int* nfla4,
               int* nad4, int* nadto4, int* mpiad4, int* mpiadlist4, int* mpl4,
               int* nir4, int* pentarptr4, int* pentarind4,
               int* niw4, int* pentawptr4, int* pentawind4, double* pentacrsval4,
               int* map4to10,
               int* cnysep_read4, int* cnysep_write4, int* nsep_n_ptr4, int* nsep4, int* nsepnmap4,
               double* uvalcoe, int* log_innerCG, int* log_innerCG_per_outite,
               int* use_crs_d, int* overlap_comm_d, int* overlap_comm_s, int* persistent_comm_s,
               double* b, double* u);
    void output_displacement_(int* im, int* n, int* ni, double* up, int* iload);
    void read_nload_(int* nload);
    void read_loads_(int* nload, double* load_arr);
    void read_surf_nelem_(int* im, int* surf_nelem);
    void read_surf_cny_(int* im, int* surf_nelem, int* surf_cny);
#ifdef __cplusplus
};
#endif

class Pointload_Solver {
public:
    Pointload_Solver(int myid, int numprocs);
    void prep();
    void exec(int iload, bool write_displacement);
    int im, ncpu, nload;
    std::vector<double> up;
    std::vector<int> load_elem_arr;
private:
    // necessary for exec()
    int n, ne, ni, surf_nelem, 
        nad, nadtoi, ierr, np, kd, nei, nird, niwd, nir, niw,
        nsep, ncolor, n4, ne4, ni4, nei4, nad4, nadtoi4, nir4,
        niw4, nsep4;
    double dsi, dsi4;
    std::vector<double> coor, rv, load_arr, mpibuf,
        younglst, sig, coori, uval, duplixyz, bcarray, 
        kcrsvald, kcrsval, coori4, uval4, duplixyz4, bcarray4,
        kcrsval4, uvalcoe;
    std::vector<int> cny, surf_cny, mpiadi, 
        mpiadlisti, mpli, num, numi, nfla, pentarptrd, pentarindd,
        pentawptrd, pentarptr, pentarind, pentawptr, pentawind, 
        numsep, nsep_n_ptr, nsep_ne_ptr, nsepnmap, cnysep_readm, 
        cnysep_write, color_ind, cny4, num4, cnyi4, numi4, nfla4,
        mpiadlisti4, mpli4, pentarptr4, pentarind4, pentawptr4,
        pentawind4, map4to10, cnysep_readm4, cnysep_write4, nsep_n_ptr4, 
        nsepnmap4, cnyi, pentawindd, mpiadi4;
    int log_innerCG, log_innerCG_per_outite, use_crs_d, overlap_comm_d,
        overlap_comm_s, persistent_comm_s;

    // other
    int nfa, n2dvisl, n2dvisa, n_size, ne_size, ib, nab, nadto,
        n4_size, ne4_size, ib4, nab4, nadto4, ni_size, nei_size, ni4_size, nei4_size;
    double ds;
    std::vector<double> plane;
    long fidc, fidm, bufls, bufds;
};
