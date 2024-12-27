// by Kohei Fujita 2015/2/11

// --- Inputs ---
//   data/tet4vox8model.h5 -- tet4+vox8 hybrid mesh
//   data/tet4vox8setting.dat -- read number of nodes and elements
//   data/para_setting.dat -- read number of mpi domains to partition mesh
//   data/microdomain.dat -- set micro domain (optional)

// --- Outputs ---
//   hdata/<000000>.data.h5 -- partitioned model
//   hdata/<000000>.part_mpi.vtk -- partitioned mesh with mpi settings
//   (optional, not used for analyses) hdata/<000000>.part_surf.vtk --
//   partitioned mesh with surface node & micro domain settings (optional, set
//   by OUTPUT_VTK_LOCAL) data/tet4vox8model.vtk -- global mesh (optional, set
//   by OUTPUT_VTK_LOCAL)

// --- Options ---
// Weights for tetra and voxel elements used for METIS_PartMeshDual
#define TETRA_WEIGHT 1
#define VOXEL_WEIGHT 6

// Output results on micro domain boundary
// reads configuration from data/microdomain.dat
// #define OUTPUT_MICRODOMAINBOUNDARY
// Number of ds that makes one block in coarsest mesh
// Boundary of micro domain is snapped to nearest block
#define DS_BLOCK_FACTOR 8

// Output vtk files of whole domain and partitioned domains
// Just for checking, not needed for analyses
#define OUTPUT_VTK_GLOBAL
// #define OUTPUT_VTK_LOCAL

// #define MODIFY_GEOMETRY
// #define INFINITE_BOUNDARY

// #define USE_PART
// for summit
// #define USE_2STAGE
// #define PART_LIQ

//-----------------------------------------------------------
// Other debug flags
// #define DEBUG_NO_HDF5
// #define DEBUG_NO_METIS
// #define DEBUG_NO_READMESH
// #define DEBUG_NO_SURFACE

//******************************************************************************
// include files
//
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef DEBUG_NO_HDF5
extern "C" {
#include "hdf_lib.c"
}
#endif

#ifdef DEBUG_NO_METIS
#define idx_t long
#else
#include "metis.h"
// #include "/home/fujita/Softwares/bin/include/metis.h"
#endif

//******************************************************************************
// macros
//
#ifndef Swap2Bytes
#define Swap2Bytes(data) ((((data) >> 8) & 0x00FF) | (((data) << 8) & 0xFF00))
#endif
#ifndef SWAP_SHORT
#define SWAP_SHORT(x) \
    (*(unsigned short*)&(x) = Swap4Bytes(*(unsigned short*)&(x)))
#endif
#ifndef Swap4Bytes
#define Swap4Bytes(data)                                            \
    ((((data) >> 24) & 0x000000FF) | (((data) >> 8) & 0x0000FF00) | \
     (((data) << 8) & 0x00FF0000) | (((data) << 24) & 0xFF000000))
#endif
#ifndef SWAP_INT
#define SWAP_INT(x) (*(unsigned int*)&(x) = Swap4Bytes(*(unsigned int*)&(x)))
#endif

#ifndef Swap8Bytes
#define Swap8Bytes(data)                     \
    ((((data) >> 56) & 0x00000000000000FF) | \
     (((data) >> 40) & 0x000000000000FF00) | \
     (((data) >> 24) & 0x0000000000FF0000) | \
     (((data) >> 8) & 0x00000000FF000000) |  \
     (((data) << 8) & 0x000000FF00000000) |  \
     (((data) << 24) & 0x0000FF0000000000) | \
     (((data) << 40) & 0x00FF000000000000) | \
     (((data) << 56) & 0xFF00000000000000))
#endif
#ifndef SWAP_DOUBLE
#define SWAP_DOUBLE(x) \
    (*(unsigned long*)&(x) = Swap8Bytes(*(unsigned long*)&(x)))
#endif
#ifndef SWAP_FLOAT
#define SWAP_FLOAT(x) (*(unsigned int*)&(x) = Swap4Bytes(*(unsigned int*)&(x)))
#endif

//******************************************************************************
// templates
//
template <class T>
void swap_bytes(T& in) {}
template <>
void swap_bytes<int>(int& in) {
    SWAP_INT(in);
}
template <>
void swap_bytes<float>(float& in) {
    SWAP_FLOAT(in);
}
template <class T>
void vec_swap_bytes(std::vector<T>& in) {
    for (long i = 0; i < in.size(); ++i) swap_bytes<T>(in[i]);
}
template <class T>
std::string ToString(const T in) {
    std::stringstream ss;
    ss << in;
    return ss.str();
}
template <class T>
std::string ToString(const T in, const long size, const char filling) {
    std::string strtmp(ToString<T>(in));
    std::string strout(size, filling);
    const long strsize(strtmp.size());
    if (strsize > size)
        return strtmp;
    else {
        for (long i = 0; i < strsize; ++i)
            strout[i + size - strsize] = strtmp[i];
        return strout;
    }
}
template <class T>
void Error(const T message) {
    std::cout << "ERROR: " << message << std::endl;
    exit(-1);
}
template <class T>
void check_open(const std::string filename, T& file) {
    file.open(filename.c_str());
    if (!file) Error("cannot find " + filename);
}

template <class T>
void check_open2(const std::string filename1, const std::string filename2,
                 T& file) {
    file.open(filename1.c_str());
    if (!file) {
        file.open(filename2.c_str());
        if (file) {
            std::cout << "cannot find " + filename1 + "; Use " + filename2 +
                             " instead."
                      << std::endl;
        } else {
            Error("cannot find " + filename1 + " nor " + filename2);
        }
    }
}

template <class T>
void check_open3(const std::string filename1, const std::string filename2,
                 T& file) {
    file.open(filename1.c_str());
    if (!file) {
        file.open(filename2.c_str());
        if (file) {
            std::cout << "cannot find " + filename1 + "; Use " + filename2 +
                             " instead."
                      << std::endl;
        }
    }
}

//******************************************************************************
// functions
//
void read_settings(long& n, long& ne_tet4, long& ne_vox8, long& nmat,
                   long& nproc, long& nthread) {
#ifdef DEBUG_NO_READMESH
    n = 9;
    ne_tet4 = 1;
    ne_vox8 = 1;
    nmat = 2;
    nproc = 2;
    nthread = 1;
#else
    std::string tmp;
    std::ifstream fin;
#ifdef INFINITE_BOUNDARY
    check_open("./data/tet4vox8setting_infinite.dat", fin);
#else
#ifdef MODIFY_GEOMETRY
    check_open("./data/tet4vox8setting_mod.dat", fin);
#else
    check_open("./data/tet4vox8setting.dat", fin);
#endif
#endif
    fin >> tmp;  // num
    fin >> tmp;  // of
    fin >> tmp;  // node
    fin >> n;
    fin >> tmp;  // num
    fin >> tmp;  // of
    fin >> tmp;  // tet
    fin >> tmp;  // elem
    fin >> ne_tet4;
    fin >> tmp;  // num
    fin >> tmp;  // of
    fin >> tmp;  // vox
    fin >> tmp;  // elem
    fin >> ne_vox8;
    fin >> tmp;  // num
    fin >> tmp;  // of
    fin >> tmp;  // material
    fin >> tmp;  // properties
    fin >> nmat;
    fin.close();

    check_open2("./data/para_setting.dat", "input/para_setting.dat", fin);
    fin >> tmp;  // num
    fin >> tmp;  // of
    fin >> tmp;  // MPI
    fin >> tmp;  // thread
    fin >> nproc;
    fin >> tmp;  // num
    fin >> tmp;  // of
    fin >> tmp;  // OpenMP
    fin >> tmp;  // thread
    fin >> nthread;
    fin.close();
#endif
}

void check_mesh_quality(long n, long ne_tet4, const std::vector<double>& coor,
                        const std::vector<long>& conn_tet4) {
    double minl, maxl, gmaxrat(0.0), tmp;
    long lconn[6][2], j, i1, i2;
    lconn[0][0] = 0;
    lconn[0][1] = 1;
    lconn[1][0] = 0;
    lconn[1][1] = 2;
    lconn[2][0] = 0;
    lconn[2][1] = 3;
    lconn[3][0] = 1;
    lconn[3][1] = 2;
    lconn[4][0] = 1;
    lconn[4][1] = 3;
    lconn[5][0] = 2;
    lconn[5][1] = 3;
    for (long ie = 0; ie < ne_tet4; ++ie) {
        j = 0;
        i1 = conn_tet4[5 * ie + lconn[j][0]];
        i2 = conn_tet4[5 * ie + lconn[j][1]];
        tmp = (coor[3 * i1 + 0] - coor[3 * i2 + 0]) *
                  (coor[3 * i1 + 0] - coor[3 * i2 + 0]) +
              (coor[3 * i1 + 1] - coor[3 * i2 + 1]) *
                  (coor[3 * i1 + 1] - coor[3 * i2 + 1]) +
              (coor[3 * i1 + 2] - coor[3 * i2 + 2]) *
                  (coor[3 * i1 + 2] - coor[3 * i2 + 2]);
        maxl = minl = tmp;
        for (j = 1; j < 6; ++j) {
            i1 = conn_tet4[5 * ie + lconn[j][0]];
            i2 = conn_tet4[5 * ie + lconn[j][1]];
            tmp = (coor[3 * i1 + 0] - coor[3 * i2 + 0]) *
                      (coor[3 * i1 + 0] - coor[3 * i2 + 0]) +
                  (coor[3 * i1 + 1] - coor[3 * i2 + 1]) *
                      (coor[3 * i1 + 1] - coor[3 * i2 + 1]) +
                  (coor[3 * i1 + 2] - coor[3 * i2 + 2]) *
                      (coor[3 * i1 + 2] - coor[3 * i2 + 2]);
            minl = std::min(minl, tmp);
            maxl = std::max(maxl, tmp);
        }
        tmp = sqrt(maxl / minl);
        gmaxrat = std::max(tmp, gmaxrat);
    }
    std::cout << "max length ratio " << gmaxrat << "\n";
    // exit(0);
}

void reorder_for_part(long n, long ne_tet4, long ne_vox8,
                      std::vector<double>& coor, std::vector<long>& conn_tet4,
                      std::vector<long>& conn_vox8, long& part_n,
                      long& part_ne_tet4) {
    if (conn_vox8.size() != 0)
        Error("reorder_for_part only for tetrahedral elements");
    // conn_tet4 // 0~4: node number, 5: material number, in cstyle
    // conn_vox8 // 0~8: node number, 9: material number, in cstyle
    double range[2][3];
    std::vector<long> conn_tet4_tmp(conn_tet4.size(), 0);
    std::vector<long> conn_vox8_tmp(conn_vox8.size(), 0);
    std::vector<int> renumberflag(conn_tet4.size() / 5 + conn_vox8.size() / 9);
    std::vector<int> mapnode(n, -1);
    long psize(0), j, ii, jj, kk;
#ifdef PART_LIQ
    std::fstream fin;
    std::vector<int> part_flag;
    std::string stmp;
    int kd, itemp;
    check_open3("./data/part_flag.dat", "input/part_flag.dat", fin);
    if (fin) {
        fin >> stmp;  // kd
        fin >> kd;
        fin >> stmp;  // part_flag
        for (int i = 0; i < kd; i++) {
            fin >> itemp;
            part_flag.push_back(itemp);
        }
        fin.close();
    } else {
        check_open("./data/tet4vox8setting.dat", fin);
        fin >> stmp >> stmp >> stmp >> stmp;
        fin >> stmp >> stmp >> stmp >> stmp >> stmp;
        fin >> stmp >> stmp >> stmp >> stmp >> stmp;
        fin >> stmp >> stmp >> stmp >> stmp;
        fin >> kd;
        std::cout << "kd " << kd << std::endl;
        for (int i = 0; i < kd; i++) {
            part_flag.push_back(0);
        }
    }

#endif

    for (long i = 0; i < ne_tet4; ++i) {
#ifdef PART_LIQ
        if (part_flag[conn_tet4[5 * i + 4]] == 1)
#else
        range[0][0] = coor[3 * conn_tet4[5 * i]];
        range[1][0] = coor[3 * conn_tet4[5 * i]];
        range[0][1] = coor[3 * conn_tet4[5 * i] + 1];
        range[1][1] = coor[3 * conn_tet4[5 * i] + 1];
        range[0][2] = coor[3 * conn_tet4[5 * i] + 2];
        range[1][2] = coor[3 * conn_tet4[5 * i] + 2];
        for (j = 0; j < 4; ++j) {
            range[0][0] = std::min(range[0][0], coor[3 * conn_tet4[5 * i + j]]);
            range[1][0] = std::max(range[1][0], coor[3 * conn_tet4[5 * i + j]]);
            range[0][1] =
                std::min(range[0][1], coor[3 * conn_tet4[5 * i + j] + 1]);
            range[1][1] =
                std::max(range[1][1], coor[3 * conn_tet4[5 * i + j] + 1]);
            range[0][2] =
                std::min(range[0][2], coor[3 * conn_tet4[5 * i + j] + 2]);
            range[1][2] =
                std::max(range[1][2], coor[3 * conn_tet4[5 * i + j] + 2]);
        }
        //    if(range[0][2] > 72.0 && range[1][2] < 77.0 )
        //    if(range[0][2] > 72.0-0.8899 && range[1][2] < 77.0+0.44 )
        if (range[0][2] > 73.77777 - 0.888888 * 2.0 - 0.1 &&
            range[1][2] < 75.55555 + 0.888888 * 2.0 + 0.1)
#endif
        {
            psize++;
            renumberflag[i] = 1;
        } else
            renumberflag[i] = 0;
    }
    ii = 0;
    jj = psize;
    for (long i = 0; i < ne_tet4; ++i) {
        if (renumberflag[i] == 1) {
            kk = ii;
            ii++;
        } else {
            kk = jj;
            jj++;
        }
        for (j = 0; j < 5; ++j)
            conn_tet4_tmp[5 * kk + j] = conn_tet4[5 * i + j];
    }
    conn_tet4 = conn_tet4_tmp;
    long pnodesize(0), itmp;
    part_n = 0;
    std::vector<double> coortmp;
    for (long i = 0; i < ne_tet4; ++i) {
        for (j = 0; j < 4; ++j) {
            itmp = conn_tet4[5 * i + j];
            if (mapnode[itmp] == -1) {
                mapnode[itmp] = pnodesize;
                coortmp.push_back(coor[3 * itmp]);
                coortmp.push_back(coor[3 * itmp + 1]);
                coortmp.push_back(coor[3 * itmp + 2]);
                pnodesize++;
                if (i < psize) part_n++;
            }
            conn_tet4[5 * i + j] = mapnode[itmp];
        }
    }
    coor = coortmp;
    part_ne_tet4 = psize;
}

void read_mesh(long n, long ne_tet4, long ne_vox8, std::string name,
               std::vector<double>& coor, std::vector<long>& conn_tet4,
               std::vector<long>& conn_vox8) {
    long fid(0), size, memsize;
    coor.resize(3 * n);
    conn_tet4.resize(
        5 * ne_tet4);  // 0~4: node number, 5: material number, in cstyle
    conn_vox8.resize(
        9 * ne_vox8);  // 0~8: node number, 9: material number, in cstyle
    size = 0;
#ifdef DEBUG_NO_READMESH
    coor[3 * 0 + 0] = 0.0;
    coor[3 * 0 + 1] = 0.0;
    coor[3 * 0 + 2] = 0.0;
    coor[3 * 1 + 0] = 1.0;
    coor[3 * 1 + 1] = 0.0;
    coor[3 * 1 + 2] = 0.0;
    coor[3 * 2 + 0] = 1.0;
    coor[3 * 2 + 1] = 1.0;
    coor[3 * 2 + 2] = 0.0;
    coor[3 * 3 + 0] = 0.0;
    coor[3 * 3 + 1] = 1.0;
    coor[3 * 3 + 2] = 0.0;
    coor[3 * 4 + 0] = 0.0;
    coor[3 * 4 + 1] = 0.0;
    coor[3 * 4 + 2] = 1.0;
    coor[3 * 5 + 0] = 1.0;
    coor[3 * 5 + 1] = 0.0;
    coor[3 * 5 + 2] = 1.0;
    coor[3 * 6 + 0] = 1.0;
    coor[3 * 6 + 1] = 1.0;
    coor[3 * 6 + 2] = 1.0;
    coor[3 * 7 + 0] = 0.0;
    coor[3 * 7 + 1] = 1.0;
    coor[3 * 7 + 2] = 1.0;
    coor[3 * 8 + 0] = 0.5;
    coor[3 * 8 + 1] = -1.0;
    coor[3 * 8 + 2] = 0.5;
    conn_vox8[0] = 0;
    conn_vox8[1] = 1;
    conn_vox8[2] = 2;
    conn_vox8[3] = 3;
    conn_vox8[4] = 4;
    conn_vox8[5] = 5;
    conn_vox8[6] = 6;
    conn_vox8[7] = 7;
    conn_vox8[8] = 0;
    conn_tet4[0] = 0;
    conn_tet4[1] = 1;
    conn_tet4[2] = 4;
    conn_tet4[3] = 8;
    conn_tet4[4] = 1;
#else
#ifndef DEBUG_NO_HDF5

#ifdef MODIFY_GEOMETRY
#ifdef INFINITE_BOUNDARY
    std::cout << "reading data/tet4vox8model_infinite.h5\n";
    hdf_open_file_("./data/tet4vox8model_infinite.h5", &fid);
#else
    std::cout << "reading data/tet4vox8model_mod.h5\n";
    hdf_open_file_("./data/tet4vox8model_mod.h5", &fid);
#endif
#else
    // std::cout << "reading data/tet4vox8model.h5\n";
    // hdf_open_file_("./data/tet4vox8model.h5", &fid);
    std::cout << "reading " << name << std::endl;
    hdf_open_file_(name.c_str(), &fid);
#endif

    memsize = 3 * n;
    hdf_read_double_array_(&fid, "/coor", &(coor[0]), &memsize, &size);
    // for (long i = 0; i < n; i++) {
    //     std::cout << coor[3 * i] << " " << coor[3 * i + 1] << " "
    //               << coor[3 * i + 2] << std::endl;
    // }
    // std::cout << " --------------- " << std::endl;

    memsize = 5 * ne_tet4;
    if (memsize > 0)
        hdf_read_long_array_(&fid, "/conn_tet4", &(conn_tet4[0]), &memsize,
                             &size);
    for (long i = 0; i < memsize; ++i) conn_tet4[i]--;

    // for (long i = 0; i < ne_tet4; i++) {
    //     std::cout << conn_tet4[5 * i] << " " << conn_tet4[5 * i + 1] << " "
    //               << conn_tet4[5 * i + 2] << " " << conn_tet4[5 * i + 3] << "
    //               "
    //               << conn_tet4[5 * i + 4] << std::endl;
    // }

    memsize = 9 * ne_vox8;
    if (memsize > 0)
        hdf_read_long_array_(&fid, "/conn_vox8", &(conn_vox8[0]), &memsize,
                             &size);
    for (long i = 0; i < memsize; ++i) conn_vox8[i]--;

    hdf_close_file_(&fid);
#endif
#endif
    check_mesh_quality(n, ne_tet4, coor, conn_tet4);
}

void create_modeldomain_dat(long n, std::vector<double>& coor) {
    std::ifstream fin;
    fin.open("./data/modeldomain.dat");
    if (fin) {
        // ./data/modeldomain.dat が存在する場合は何もしない
        fin.close();
        return;
    }
    std::cout << "./data/modeldomain.dat is not found. Creating..."
              << std::endl;
    double minx, maxx, miny, maxy, minz;
    double min2x;
    minx = maxx = coor[0];
    miny = maxy = coor[1];
    minz = coor[2];
    for (long i = 0; i < n; ++i) {
        minx = std::min(minx, coor[3 * i + 0]);
        maxx = std::max(maxx, coor[3 * i + 0]);
        miny = std::min(miny, coor[3 * i + 1]);
        maxy = std::max(maxy, coor[3 * i + 1]);
        minz = std::min(minz, coor[3 * i + 2]);
    }
    std::ofstream fout;
    check_open("./data/modeldomain.dat", fout);
    fout << std::fixed << std::setprecision(10);
    fout << "xmin-plane" << std::endl;
    fout << minx << std::endl;
    fout << "xmax-plane" << std::endl;
    fout << maxx << std::endl;
    fout << "ymin-plane" << std::endl;
    fout << miny << std::endl;
    fout << "ymax-plane" << std::endl;
    fout << maxy << std::endl;
    fout << "zmin-plane" << std::endl;
    fout << minz << std::endl;
    fout.close();
}

void partition_mesh_nodal(long nproc, long n,
                          const std::vector<long>& conn_tet4,
                          const std::vector<long>& conn_vox8,
                          std::vector<idx_t>& epart) {
    std::cout << " --- using METIS_PartMeshNodal\n";
    long ne_tet4(conn_tet4.size() / 5);
    long ne_vox8(conn_vox8.size() / 9);

    idx_t ne(ne_vox8 + ne_tet4);
    idx_t nn(n);
    epart.resize(ne);
#ifndef DEBUG_NO_METIS
    if (nproc == 1) {
        for (long i = 0; i < ne; ++i) epart[i] = 0;
    } else {
        std::vector<idx_t> eptr(ne + 1);
        std::vector<idx_t> eind(ne_tet4 * 4 + ne_vox8 * 8);
        idx_t* vwgt(NULL);
        idx_t* vsize(NULL);
        idx_t nparts(nproc);
        real_t* tpwgts(NULL);
        idx_t options[METIS_NOPTIONS];
        idx_t objval;
        std::vector<idx_t> npart(nn);

        METIS_SetDefaultOptions(options);
        options[0] = 1;   //! METIS_OPTION_PTYPE 0: recursive bisection, 1:Kway
        options[6] = 20;  //! METIS_OPTION_NITER default=10
        options[16] = 5;  //! 10 ! METIS_OPTION_UFACTOR default=30

        for (long itmp = 0; itmp < ne_vox8; ++itmp) {
            eptr[itmp] = itmp * 8;
            eind[8 * itmp] = conn_vox8[9 * itmp];
            eind[8 * itmp + 1] = conn_vox8[9 * itmp + 1];
            eind[8 * itmp + 2] = conn_vox8[9 * itmp + 2];
            eind[8 * itmp + 3] = conn_vox8[9 * itmp + 3];
            eind[8 * itmp + 4] = conn_vox8[9 * itmp + 4];
            eind[8 * itmp + 5] = conn_vox8[9 * itmp + 5];
            eind[8 * itmp + 6] = conn_vox8[9 * itmp + 6];
            eind[8 * itmp + 7] = conn_vox8[9 * itmp + 7];
        }
        long ii = 8 * ne_vox8;
        for (long itmp = 0; itmp < ne_tet4; ++itmp) {
            eptr[itmp + ne_vox8] = ii + 4 * itmp;
            eind[4 * itmp + ii] = conn_tet4[5 * itmp];
            eind[4 * itmp + ii + 1] = conn_tet4[5 * itmp + 1];
            eind[4 * itmp + ii + 2] = conn_tet4[5 * itmp + 2];
            eind[4 * itmp + ii + 3] = conn_tet4[5 * itmp + 3];
        }
        eptr[ne_tet4 + ne_vox8] = ii + ne_tet4 * 4;

        int errcheck = METIS_PartMeshNodal(
            &ne, &nn, &(eptr[0]), &(eind[0]), vwgt, vsize, &nparts, tpwgts,
            options, &objval, &(epart[0]), &(npart[0]));
        if (errcheck != METIS_OK) Error("Metis failed");
    }
#else
    for (long i = 0; i < ne; ++i) epart[i] = 0;
    epart[1] = 1;
#endif
}

void partition_mesh_dual(long nproc, long n, const std::vector<long>& conn_tet4,
                         const std::vector<long>& conn_vox8,
                         const std::vector<double>& coor,
                         std::vector<idx_t>& epart) {
    int weight_vox8(VOXEL_WEIGHT);
    int weight_tet4(TETRA_WEIGHT);
    std::cout << " --- using METIS_PartMeshDual\n";
    std::cout << " --- weight vox,tet: " << weight_vox8 << " " << weight_tet4
              << "\n";

    long ne_tet4(conn_tet4.size() / 5);
    long ne_vox8(conn_vox8.size() / 9);

    idx_t ne(ne_vox8 + ne_tet4);
    idx_t nn(n);
    epart.resize(ne);
#ifndef DEBUG_NO_METIS
    if (nproc == 1) {
        for (long i = 0; i < ne; ++i) epart[i] = 0;
    } else {
        std::vector<idx_t> eptr(ne + 1);
        std::vector<idx_t> eind(ne_tet4 * 4 + ne_vox8 * 8);
        std::vector<idx_t> vwgt(ne);
        idx_t* vsize(NULL);
        idx_t ncommon(3);  // number of nodes to be shared among elements to
                           // have connection in element-element graph
        idx_t nparts(nproc);
        real_t* tpwgts(NULL);
        idx_t options[METIS_NOPTIONS];
        idx_t objval;
        std::vector<idx_t> npart(nn);

        for (long itmp = 0; itmp < ne_vox8; ++itmp) {
            eptr[itmp] = itmp * 8;
            eind[8 * itmp] = conn_vox8[9 * itmp];
            eind[8 * itmp + 1] = conn_vox8[9 * itmp + 1];
            eind[8 * itmp + 2] = conn_vox8[9 * itmp + 2];
            eind[8 * itmp + 3] = conn_vox8[9 * itmp + 3];
            eind[8 * itmp + 4] = conn_vox8[9 * itmp + 4];
            eind[8 * itmp + 5] = conn_vox8[9 * itmp + 5];
            eind[8 * itmp + 6] = conn_vox8[9 * itmp + 6];
            eind[8 * itmp + 7] = conn_vox8[9 * itmp + 7];
        }
        long ii = 8 * ne_vox8;
#ifdef INFINITE_BOUNDARY
        long inode, iii;
        double maxx, maxy, minx, miny, limx, limy;
        limx = coor[0];
        limy = coor[1];
        for (inode = 0; inode < coor.size() / 3; ++inode) {
            limx = std::max(limx, coor[3 * inode]);
            limy = std::max(limy, coor[3 * inode + 1]);
        }
        limx *= 0.5;
        limy *= 0.5;
#endif
        for (long itmp = 0; itmp < ne_tet4; ++itmp) {
#ifdef INFINITE_BOUNDARY
            inode = conn_tet4[5 * itmp];
            maxx = minx = coor[3 * inode];
            maxy = miny = coor[3 * inode + 1];
            for (iii = 1; iii < 4; ++iii) {
                inode = conn_tet4[5 * itmp + iii];
                maxx = std::max(maxx, coor[3 * inode]);
                maxy = std::max(maxy, coor[3 * inode + 1]);
                minx = std::min(minx, coor[3 * inode]);
                miny = std::min(miny, coor[3 * inode + 1]);
            }
            if ((maxx - minx) > limx || (maxy - miny) > limy) {
                eptr[itmp + ne_vox8] = ii + 4 * itmp;
                eind[4 * itmp + ii] = conn_tet4[5 * itmp];
                eind[4 * itmp + ii + 1] = conn_tet4[5 * itmp];
                eind[4 * itmp + ii + 2] = conn_tet4[5 * itmp];
                eind[4 * itmp + ii + 3] = conn_tet4[5 * itmp];
            } else {
#endif
                eptr[itmp + ne_vox8] = ii + 4 * itmp;
                eind[4 * itmp + ii] = conn_tet4[5 * itmp];
                eind[4 * itmp + ii + 1] = conn_tet4[5 * itmp + 1];
                eind[4 * itmp + ii + 2] = conn_tet4[5 * itmp + 2];
                eind[4 * itmp + ii + 3] = conn_tet4[5 * itmp + 3];
#ifdef INFINITE_BOUNDARY
            }
#endif
        }
        eptr[ne_tet4 + ne_vox8] = ii + ne_tet4 * 4;

        for (long i = 0; i < ne_vox8; ++i) vwgt[i] = weight_vox8;
        for (long i = 0; i < ne_tet4; ++i) vwgt[i + ne_vox8] = weight_tet4;

        METIS_SetDefaultOptions(options);
        options[0] = 1;   //! METIS_OPTION_PTYPE 0: recursive bisection, 1:Kway
        options[6] = 20;  //! METIS_OPTION_NITER default=10
        options[16] = 5;  //! 10 ! METIS_OPTION_UFACTOR default=30

        int errcheck = METIS_PartMeshDual(
            &ne, &nn, &(eptr[0]), &(eind[0]), &(vwgt[0]), vsize, &ncommon,
            &nparts, tpwgts, options, &objval, &(epart[0]), &(npart[0]));
        if (errcheck != METIS_OK) Error("Metis failed");
    }
#else
    for (long i = 0; i < ne; ++i) epart[i] = 0;
    epart[1] = 1;
#endif
}

void partition_mesh_dual_2stage(long nproc, long n,
                                const std::vector<long>& conn_tet4,
                                const std::vector<long>& conn_vox8,
                                const std::vector<double>& coor,
                                std::vector<idx_t>& epart) {
    int blocksize;
    std::ifstream fin;
    check_open("data/blocksize.dat", fin);
    fin >> blocksize;
    fin.close();
    int nblock(nproc / blocksize);
    std::cout << " --- partition 2stage with blocksize " << blocksize << "\n";

    if (nproc % blocksize != 0) Error("nproc cannot be divided by blocksize");
    std::vector<idx_t> epart_block;
    partition_mesh_dual(nblock, n, conn_tet4, conn_vox8, coor, epart_block);

    std::vector<long> conn_tet4_part, conn_vox8_part;
    std::vector<double> coor_part;
    std::vector<idx_t> epart_part;
    std::vector<std::vector<long> > partition_map_tet4;
    std::vector<std::vector<long> > partition_map_vox8;
    long ne_tet4(conn_tet4.size() / 5);
    long ne_vox8(conn_vox8.size() / 9);
    long iproc, inode, ielem, itmp, jtmp, iblock, i, j, k, ii;
    bool flag;
    for (i = 0; i < ne_tet4 + ne_vox8; ++i) epart.push_back(-100000000);

    //  epart=epart_block; return;

    partition_map_vox8.resize(nblock);
    for (ielem = 0; ielem < ne_vox8; ++ielem)
        partition_map_vox8[epart_block[ielem]].push_back(ielem);

    partition_map_tet4.resize(nblock);
    for (ielem = ne_vox8; ielem < ne_tet4 + ne_vox8; ++ielem)
        partition_map_tet4[epart_block[ielem]].push_back(ielem - ne_vox8);

    std::vector<long> node_map_gtol(n, -1);
    std::vector<long> node_map_ltog;
    std::vector<long> elem_map_ltog;

    for (iblock = 0; iblock < nblock; ++iblock) {
        for (ii = 0; ii < partition_map_vox8[iblock].size(); ++ii) {
            i = partition_map_vox8[iblock][ii];
            for (j = 0; j < 8; ++j) {
                k = conn_vox8[9 * i + j];
                if (node_map_gtol[k] == -1) {
                    node_map_gtol[k] = node_map_ltog.size();
                    node_map_ltog.push_back(k);
                    coor_part.push_back(coor[3 * k]);
                    coor_part.push_back(coor[3 * k + 1]);
                    coor_part.push_back(coor[3 * k + 2]);
                }
                conn_vox8_part.push_back(node_map_gtol[k]);
            }
            conn_vox8_part.push_back(conn_vox8[9 * i + 8]);
            elem_map_ltog.push_back(i);
        }
        for (ii = 0; ii < partition_map_tet4[iblock].size(); ++ii) {
            i = partition_map_tet4[iblock][ii];
            for (j = 0; j < 4; ++j) {
                k = conn_tet4[5 * i + j];
                if (node_map_gtol[k] == -1) {
                    node_map_gtol[k] = node_map_ltog.size();
                    node_map_ltog.push_back(k);
                    coor_part.push_back(coor[3 * k]);
                    coor_part.push_back(coor[3 * k + 1]);
                    coor_part.push_back(coor[3 * k + 2]);
                }
                conn_tet4_part.push_back(node_map_gtol[k]);
            }
            conn_tet4_part.push_back(conn_tet4[5 * i + 4]);
            elem_map_ltog.push_back(i + ne_vox8);
        }
        partition_mesh_dual(blocksize, node_map_ltog.size(), conn_tet4_part,
                            conn_vox8_part, coor_part, epart_part);
        for (i = 0; i < epart_part.size(); ++i) {
            epart[elem_map_ltog[i]] = blocksize * iblock + epart_part[i];
        }
        for (ii = 0; ii < partition_map_vox8[iblock].size(); ++ii) {
            i = partition_map_vox8[iblock][ii];
            for (j = 0; j < 8; ++j) {
                k = conn_vox8[9 * i + j];
                node_map_gtol[k] = -1;
            }
        }
        for (ii = 0; ii < partition_map_tet4[iblock].size(); ++ii) {
            i = partition_map_tet4[iblock][ii];
            for (j = 0; j < 4; ++j) {
                k = conn_tet4[5 * i + j];
                node_map_gtol[k] = -1;
            }
        }
        conn_vox8_part.clear();
        conn_tet4_part.clear();
        coor_part.clear();
        node_map_ltog.clear();
        epart_part.clear();
        elem_map_ltog.clear();
    }
}

void partition_mesh(long nproc, long n, const std::vector<long>& conn_tet4,
                    const std::vector<long>& conn_vox8,
                    const std::vector<double>& coor, long part_n,
                    long part_ne_tet4,
                    std::vector<std::vector<long> >& partition_map_tet4,
                    std::vector<std::vector<long> >& partition_map_vox8,
                    std::vector<std::vector<int> >& node_duplicate_map) {
    std::vector<idx_t> epart;
#ifdef USE_PART
    std::vector<idx_t> epart2;
    std::vector<double> coorpart;
    std::vector<long> conn_tet4_part;
    std::vector<long> conn_vox8_part;

    for (long i = 0; i < 3 * part_n; ++i) coorpart.push_back(coor[i]);
    for (long i = 0; i < 5 * part_ne_tet4; ++i)
        conn_tet4_part.push_back(conn_tet4[i]);
#ifdef USE_2STAGE
    partition_mesh_dual_2stage(nproc, part_n, conn_tet4_part, conn_vox8_part,
                               coorpart, epart);
#else
    partition_mesh_dual(nproc, part_n, conn_tet4_part, conn_vox8_part, coorpart,
                        epart);
#endif

    std::vector<long> mapnode(n, -1);
    coorpart.clear();
    conn_tet4_part.clear();
    long jj, iitmp(0);
    for (long i = part_ne_tet4; i < conn_tet4.size() / 5; ++i) {
        for (long j = 0; j < 4; ++j) {
            jj = conn_tet4[5 * i + j];
            if (mapnode[jj] == -1) {
                mapnode[jj] = iitmp;
                ++iitmp;
                coorpart.push_back(coor[3 * jj + 0]);
                coorpart.push_back(coor[3 * jj + 1]);
                coorpart.push_back(coor[3 * jj + 2]);
            }
            conn_tet4_part.push_back(mapnode[jj]);
        }
        conn_tet4_part.push_back(conn_tet4[5 * i + 4]);
    }
#ifdef USE_2STAGE
    partition_mesh_dual_2stage(nproc, iitmp, conn_tet4_part, conn_vox8_part,
                               coorpart, epart2);
#else
    partition_mesh_dual(nproc, iitmp, conn_tet4_part, conn_vox8_part, coorpart,
                        epart2);
#endif

    for (long i = 0; i < conn_tet4.size() / 5 - part_ne_tet4; ++i)
        epart.push_back(epart2[i]);
//  partition_mesh_dual(nproc,n,conn_tet4,conn_vox8,coor,epart);
#else
    //  partition_mesh_nodal(nproc,n,conn_tet4,conn_vox8,epart);
    partition_mesh_dual(nproc, n, conn_tet4, conn_vox8, coor, epart);
#endif
    // partition_map_tet4
    // partition_map_vox8
    //  map[iproc][localelemid]-> globalelemid

    long ne_tet4(conn_tet4.size() / 5);
    long ne_vox8(conn_vox8.size() / 9);
    long iproc, inode, ielem, itmp, jtmp;
    bool flag;

    partition_map_vox8.resize(nproc);
    for (ielem = 0; ielem < ne_vox8; ++ielem)
        partition_map_vox8[epart[ielem]].push_back(ielem);

    partition_map_tet4.resize(nproc);
    for (ielem = ne_vox8; ielem < ne_tet4 + ne_vox8; ++ielem)
        partition_map_tet4[epart[ielem]].push_back(ielem - ne_vox8);

    // node_duplicate_map[globalnodeid].size()-> number of duplicates,
    // map[globalnodeid][]-> iproc
    node_duplicate_map.resize(n);
    for (long ielem = 0; ielem < ne_vox8; ++ielem) {
        iproc = epart[ielem];
        for (itmp = 0; itmp < 8; ++itmp) {
            inode = conn_vox8[9 * ielem + itmp];
            flag = true;
            for (jtmp = 0; jtmp < node_duplicate_map[inode].size(); ++jtmp) {
                if (node_duplicate_map[inode][jtmp] == iproc) flag = false;
            }
            if (flag) node_duplicate_map[inode].push_back(iproc);
        }
    }
    for (long ielem = 0; ielem < ne_tet4; ++ielem) {
        iproc = epart[ielem + ne_vox8];
        for (itmp = 0; itmp < 4; ++itmp) {
            inode = conn_tet4[5 * ielem + itmp];
            flag = true;
            for (jtmp = 0; jtmp < node_duplicate_map[inode].size(); ++jtmp) {
                if (node_duplicate_map[inode][jtmp] == iproc) flag = false;
            }
            if (flag) node_duplicate_map[inode].push_back(iproc);
        }
    }
}

template <class T1, class T2>
void write_vtk_hybrid_mesh(const std::vector<double>& coor_in,
                           const std::vector<T1>& intconn_tet4,
                           const std::vector<T1>& intconn_vox8,
                           const std::vector<T2>& node_intscalar,
                           const std::vector<double>& node_doublevec,
                           std::string name) {
    long nodes(coor_in.size() / 3);
    long elements_tet4(intconn_tet4.size() / 5);
    long elements_vox8(intconn_vox8.size() / 9);
    long elements(elements_tet4 + elements_vox8);
    long elementcomponents(5 * elements_tet4 + 9 * elements_vox8);
    long i, j;

    std::vector<float> coor(3 * nodes);
    std::vector<int> conn(elementcomponents);
    std::vector<int> inttype(elements);
    std::vector<int> matid(elements);
    std::vector<int> intscalar(nodes);
    std::vector<float> floatvec(3 * nodes);

    for (i = 0; i < nodes; ++i) {
        coor[3 * i] = coor_in[3 * i];
        coor[3 * i + 1] = coor_in[3 * i + 1];
        coor[3 * i + 2] = coor_in[3 * i + 2];
#ifdef USE_PART
        intscalar[i] = i;
#else
        intscalar[i] = node_intscalar[i];
#endif
        floatvec[3 * i] = node_doublevec[3 * i];
        floatvec[3 * i + 1] = node_doublevec[3 * i + 1];
        floatvec[3 * i + 2] = node_doublevec[3 * i + 2];
    }

    for (i = 0; i < elements_vox8; ++i) {
        conn[9 * i] = 8;
        conn[9 * i + 1] = intconn_vox8[9 * i];
        conn[9 * i + 2] = intconn_vox8[9 * i + 1];
        conn[9 * i + 3] = intconn_vox8[9 * i + 2];
        conn[9 * i + 4] = intconn_vox8[9 * i + 3];
        conn[9 * i + 5] = intconn_vox8[9 * i + 4];
        conn[9 * i + 6] = intconn_vox8[9 * i + 5];
        conn[9 * i + 7] = intconn_vox8[9 * i + 6];
        conn[9 * i + 8] = intconn_vox8[9 * i + 7];
        matid[i] = intconn_vox8[9 * i + 8];
        inttype[i] = 12;  // 12 for hexahedral elements
    }
    for (i = 0; i < elements_tet4; ++i) {
        conn[9 * elements_vox8 + 5 * i] = 4;
        for (j = 0; j < 4; ++j)
            conn[9 * elements_vox8 + 5 * i + j + 1] = intconn_tet4[5 * i + j];
#ifdef USE_PART
        matid[elements_vox8 + i] = i;
#else
        matid[elements_vox8 + i] = intconn_tet4[5 * i + 4];
#endif
        inttype[elements_vox8 + i] = 10;  // 10 for tet4 elements
    }

    vec_swap_bytes(coor);
    vec_swap_bytes(conn);
    vec_swap_bytes(intscalar);
    vec_swap_bytes(floatvec);
    vec_swap_bytes(matid);
    vec_swap_bytes(inttype);

    std::ofstream fout;
    check_open(name, fout);
    fout << "# vtk DataFile Version 2.0\n";
    fout << "file: " << name << std::endl;
    fout << "BINARY\n";
    fout << "DATASET UNSTRUCTURED_GRID\n";
    fout << "POINTS " << nodes << " float\n";
    fout.write((char*)&(coor[0]), sizeof(float) * 3 * nodes);
    fout << "CELLS " << elements << " " << elementcomponents << std::endl;
    fout.write((char*)&(conn[0]), sizeof(int) * elementcomponents);
    fout << "CELL_TYPES " << elements << std::endl;
    fout.write((char*)&(inttype[0]), sizeof(int) * elements);
    fout << "CELL_DATA " << elements << std::endl;
    fout << "SCALARS material int 1\n";
    fout << "LOOKUP_TABLE default\n";
    fout.write((char*)&(matid[0]), sizeof(int) * elements);
    fout << "POINT_DATA " << nodes << "\n";
    fout << "VECTORS floatvec float\n";
    fout.write((char*)&(floatvec[0]), sizeof(float) * 3 * nodes);
    fout << "SCALARS intvec int 1\n";
    fout << "LOOKUP_TABLE default\n";
    fout.write((char*)&(intscalar[0]), sizeof(float) * nodes);
    fout.close();
}

template <class T>
void write_vtk_nodes(const std::vector<double>& coor_in,
                     const std::vector<T>& node_intscalar, std::string name) {
    // output nodes that have intscalar[i] >= 0

    long nodes(coor_in.size() / 3);
    std::vector<float> coor;
    for (long i = 0; i < nodes; ++i) {
        if (node_intscalar[i] >= 0) {
            coor.push_back(coor_in[3 * i]);
            coor.push_back(coor_in[3 * i + 1]);
            coor.push_back(coor_in[3 * i + 2]);
        }
    }
    nodes = coor.size() / 3;
    std::vector<int> conn(nodes * 2);
    std::vector<int> inttype(nodes, 1);
    long elements(nodes);
    long elementcomponents(2 * elements);
    long i, j;

    for (i = 0; i < nodes; ++i) {
        conn[2 * i] = 1;
        conn[2 * i + 1] = i;
    }

    vec_swap_bytes(coor);
    vec_swap_bytes(conn);
    vec_swap_bytes(inttype);

    std::ofstream fout;
    check_open(name, fout);
    fout << "# vtk DataFile Version 2.0\n";
    fout << "file: " << name << std::endl;
    fout << "BINARY\n";
    fout << "DATASET UNSTRUCTURED_GRID\n";
    fout << "POINTS " << nodes << " float\n";
    fout.write((char*)&(coor[0]), sizeof(float) * 3 * nodes);
    fout << "CELLS " << elements << " " << elementcomponents << std::endl;
    fout.write((char*)&(conn[0]), sizeof(int) * elementcomponents);
    fout << "CELL_TYPES " << elements << std::endl;
    fout.write((char*)&(inttype[0]), sizeof(int) * elements);
    fout.close();
}

void search_surface_nodes(const long n, const std::vector<double>& coor,
                          bool microdomainflag, long& nx, long& ny,
                          std::vector<long>& nodemap_surf2global,
                          std::vector<int>& nodemap_global2surf) {
#ifndef DEBUG_NO_SURFACE
    double minx, maxx, miny, maxy, minz;
    double min2x;
    minx = maxx = coor[0];
    miny = maxy = coor[1];
    minz = coor[2];
    for (long i = 0; i < n; ++i) {
        minx = std::min(minx, coor[3 * i + 0]);
        maxx = std::max(maxx, coor[3 * i + 0]);
        miny = std::min(miny, coor[3 * i + 1]);
        maxy = std::max(maxy, coor[3 * i + 1]);
        minz = std::min(minz, coor[3 * i + 2]);
    }
    min2x = maxx;
    for (long i = 0; i < n; ++i) {
        if (minx + 0.001 > coor[3 * i]) continue;
        min2x = std::min(min2x, coor[3 * i]);
    }
    // double ds(min2x-minx);
    double ds;

    std::ifstream fds;
    std::string stmp;
    fds.open("data/modeling_setting.dat");
    if (fds) {
        fds >> stmp >> stmp >> stmp >> stmp >> stmp;
        fds >> ds;
        fds.close();
    } else {
        check_open("data/ds_2Doutput.dat", fds);
        fds >> stmp;
        fds >> ds;
        fds.close();
    }

    nx = (round((maxx - minx) / ds) + 1);
    ny = (round((maxy - miny) / ds) + 1);
    std::cout << "--- start surface node info ---" << std::endl;
    std::cout << "ds: " << ds << std::endl;
    std::cout << "minx,maxx,nx: " << minx << " " << maxx << " " << nx
              << std::endl;
    std::cout << "miny,maxy,ny: " << miny << " " << maxy << " " << ny
              << std::endl;
    std::cout << "---  end  surface node info ---" << std::endl;

    nodemap_global2surf.resize(n);
    for (long i = 0; i < n; ++i) nodemap_global2surf[i] = -1;
    nodemap_surf2global.resize(nx * ny);
    for (long i = 0; i < nx * ny; ++i) nodemap_surf2global[i] = -1;
    std::vector<double> surfnodezval(nx * ny, -10000000.0);
    long indexx, indexy;
    double xtmp, ytmp;
    for (long i = 0; i < n; ++i) {
        indexx = round((coor[3 * i + 0] - minx) / ds);
        indexy = round((coor[3 * i + 1] - miny) / ds);
        xtmp = indexx * ds + minx;
        ytmp = indexy * ds + miny;
        if (fabs(xtmp - coor[3 * i + 0]) > 0.001 * ds ||
            fabs(ytmp - coor[3 * i + 1]) > 0.001 * ds)
            continue;
        // std::cout << indexx << " " << indexy << " "<< coor[3*i] << " " <<
        // coor[3*i+1] << std::endl;
        if (coor[3 * i + 2] > surfnodezval[indexy * nx + indexx]) {
            surfnodezval[indexy * nx + indexx] = coor[3 * i + 2];
            nodemap_surf2global[indexy * nx + indexx] = i;
        }
    }
    for (long i = 0; i < nx * ny; ++i) {
        if (nodemap_surf2global[i] == -1) {
            std::cout << "something wrong in search_surface_nodes_conn: "
                         "surfnodeindex\n";
            exit(-1);
        }
        nodemap_global2surf[nodemap_surf2global[i]] = i;
    }
    if (microdomainflag) {
        double md_minx, md_maxx, md_miny, md_maxy, md_minz;
        double ds_block(ds * DS_BLOCK_FACTOR);
        std::ifstream fin;
        std::string temp;
        check_open("data/microdomain.dat", fin);
        fin >> temp;  // microdomain-boundary
        fin >> temp;  // minx
        fin >> md_minx;
        fin >> temp;  // maxx
        fin >> md_maxx;
        fin >> temp;  // miny
        fin >> md_miny;
        fin >> temp;  // maxy
        fin >> md_maxy;
        fin >> temp;  // minz
        fin >> md_minz;
        fin.close();

        // fit to ds_block
        long md_nx, md_ny;
        md_minx = round((md_minx - minx) / ds_block) * ds_block + minx;
        md_miny = round((md_miny - miny) / ds_block) * ds_block + miny;
        md_minz = round((md_minz - minz) / ds_block) * ds_block + minz;
        md_nx = round((md_maxx - md_minx) / ds_block) + 1;
        md_ny = round((md_maxy - md_miny) / ds_block) + 1;
        md_maxx = (md_nx - 1) * ds_block + md_minx;
        md_maxy = (md_ny - 1) * ds_block + md_miny;
        std::cout << "--- start micro domain info ---" << std::endl;
        std::cout << "ds,ds_block,blockfactor: " << ds << " " << ds_block << " "
                  << DS_BLOCK_FACTOR << std::endl;
        std::cout << "minx,maxx,nx_block: " << md_minx << " " << md_maxx << " "
                  << md_nx << std::endl;
        std::cout << "miny,maxy,ny_block: " << md_miny << " " << md_maxy << " "
                  << md_ny << std::endl;
        std::cout << "minz: " << md_minz << std::endl;
        std::cout << "---  end  micro domain info ---" << std::endl;

        bool flag;
        double x, y, z;
        for (long i = 0; i < n; ++i) {
            flag = false;
            x = coor[3 * i + 0];
            y = coor[3 * i + 1];
            z = coor[3 * i + 2];
            if (x > md_minx - 0.01 * ds && x < md_maxx + 0.01 * ds &&
                z > md_minz - 0.01 * ds) {
                if (fabs(y - md_miny) < 0.01 * ds) flag = true;
                if (fabs(y - md_maxy) < 0.01 * ds) flag = true;
            }
            if (y > md_miny - 0.01 * ds && y < md_maxy + 0.01 * ds &&
                z > md_minz - 0.01 * ds) {
                if (fabs(x - md_minx) < 0.01 * ds) flag = true;
                if (fabs(x - md_maxx) < 0.01 * ds) flag = true;
            }
            if (y > md_miny - 0.01 * ds && y < md_maxy + 0.01 * ds &&
                x > md_minx - 0.01 * ds && x < md_maxx + 0.01 * ds &&
                fabs(z - md_minz) < 0.01 * ds)
                flag = true;
            if (flag) {
                nodemap_global2surf[i] = nodemap_surf2global.size();
                nodemap_surf2global.push_back(i);
            }
        }
    }
#else
    std::cout << "DEBUG:no surface nodes\n";
    nx = ny = 0;
    nodemap_global2surf.resize(n);
    for (long i = 0; i < n; ++i) nodemap_global2surf[i] = -1;
    nodemap_surf2global.clear();
#endif
}

void write_data_hdf(const std::vector<double>& coor_part,
                    const std::vector<int>& conn_tet4_part,
                    const std::vector<int>& conn_vox8_part,
                    const std::vector<int>& neighborrankid,
                    const std::vector<int>& mpinode_pointer,
                    const std::vector<int>& mpinode_index,
                    const std::vector<int>& surfnodes, std::string name) {
#ifndef DEBUG_NO_HDF5
    long fid(0), size;
    hdf_create_file_(name.c_str(), &fid);

    int settings[7];
    settings[0] = coor_part.size() / 3;       // number of nodes
    settings[1] = conn_vox8_part.size() / 9;  // number of vox8 elements
    settings[2] = conn_tet4_part.size() / 5;  // number of tet4 elements
    settings[3] = surfnodes.size();           // number of surface nodes
    settings[4] = neighborrankid.size();      // number of neighboring ranks
    settings[5] = mpinode_pointer.size();     // number of neighboring ranks
    settings[6] = mpinode_index.size();       // number of nodes to send/recv
    size = 7;
    hdf_write_int_array_(&fid, "/settings", settings, &size);

    size = coor_part.size();
    hdf_write_double_array_(&fid, "/coor", (double*)&(coor_part[0]), &size);

    size = conn_vox8_part.size();
    hdf_write_int_array_(&fid, "/conn_vox8", (int*)&(conn_vox8_part[0]), &size);

    size = conn_tet4_part.size();
    hdf_write_int_array_(&fid, "/conn_tet4", (int*)&(conn_tet4_part[0]), &size);

    size = surfnodes.size();
    hdf_write_int_array_(&fid, "/surfnodeid", (int*)&(surfnodes[0]), &size);

    size = neighborrankid.size();
    hdf_write_int_array_(&fid, "/neighborrankid", (int*)&(neighborrankid[0]),
                         &size);

    size = mpinode_pointer.size();
    hdf_write_int_array_(&fid, "/mpinode_pointer", (int*)&(mpinode_pointer[0]),
                         &size);

    size = mpinode_index.size();
    hdf_write_int_array_(&fid, "/mpinode_index", (int*)&(mpinode_index[0]),
                         &size);
    hdf_close_file_(&fid);
#endif
}

//******************************************************************************
// main function
//
int main() {
    double t1, t0;
    t0 = omp_get_wtime();
    // read settings
    long n, ne_tet4, ne_vox4, nmat;
    long nproc, nthread;
    read_settings(n, ne_tet4, ne_vox4, nmat, nproc,
                  nthread);  // nthread not used at present
    std::cout << "Start partitioning\n";
    std::cout << " --- n,ne_tet,ne_vox: " << n << " " << ne_tet4 << " "
              << ne_vox4 << "\n";
    std::cout << " --- nmat,nproc: " << nmat << " " << nproc << "\n";

    // read mesh
    t1 = omp_get_wtime();
    std::cout << "reading mesh...\n";
    std::vector<double> coor;
    std::vector<long> conn_tet4;
    std::vector<long> conn_vox8;

    read_mesh(n, ne_tet4, ne_vox4, "./data/tet4vox8model.h5", coor, conn_tet4,
              conn_vox8);
    std::cout << "...took " << omp_get_wtime() - t1 << "\n";

    // create ./data/modeldomain.dat (input for creation of cdata, mdata)
    create_modeldomain_dat(n, coor);

    long part_n(0), part_ne_tet4(0);
    std::vector<long> psizelocal(nproc * 2, 0);
#ifdef USE_PART
    reorder_for_part(n, ne_tet4, ne_vox4, coor, conn_tet4, conn_vox8, part_n,
                     part_ne_tet4);
    std::cout << "USE_PART part_n,n " << part_n << " " << n << " \n";
    std::cout << "USE_PART part_ne,ne " << part_ne_tet4 << " " << ne_tet4
              << " \n";
#endif

    // extract surface nodes
    t1 = omp_get_wtime();
    std::cout << "extracting surface nodes...\n";
    long nx, ny;  // number of nodes in x & y directions
    std::vector<long> nodemap_surf2global;  // map[iy*nx+ix]->global node id
    std::vector<int>
        nodemap_global2surf;  // map[global node id]->surf node id (iy*nx+ix)
    bool microdomainflag(false);
#ifdef OUTPUT_MICRODOMAINBOUNDARY
    microdomainflag = true;
#endif
    search_surface_nodes(n, coor, microdomainflag, nx, ny, nodemap_surf2global,
                         nodemap_global2surf);
    write_vtk_nodes(coor, nodemap_global2surf, "./data/surface_nodes.vtk");
    std::cout << "...took " << omp_get_wtime() - t1 << "\n";

#ifdef OUTPUT_VTK_GLOBAL
    // extract surface nodes
    t1 = omp_get_wtime();
    std::cout << "writing global vtk mesh...\n";
    std::vector<int> surfflag(n, 0);
    std::vector<double> nodevec(3 * n, 0.0);
    for (long itmp = 0; itmp < nodemap_surf2global.size(); ++itmp)
        surfflag[nodemap_surf2global[itmp]] = 1;
    write_vtk_hybrid_mesh(coor, conn_tet4, conn_vox8, surfflag, nodevec,
                          "data/tet4vox8model.vtk");
    std::cout << "...took " << omp_get_wtime() - t1 << "\n";
#endif

    // partition mesh using metis
    t1 = omp_get_wtime();
    std::cout << "partitioning...\n";
    std::vector<std::vector<long> >
        partition_map_tet4;  // map[iproc][localelemid]-> globalelemid
    std::vector<std::vector<long> > partition_map_vox8;
    std::vector<std::vector<int> >
        node_duplicate_map;  // map[globalnodeid].size()-> number of duplicates,
                             // map[globalnodeid][]-> iproc
    partition_mesh(nproc, n, conn_tet4, conn_vox8, coor, part_n, part_ne_tet4,
                   partition_map_tet4, partition_map_vox8, node_duplicate_map);
    std::cout << "...took " << omp_get_wtime() - t1 << "\n";

    std::vector<double> coor_part;
    std::vector<int> conn_tet4_part;
    std::vector<int> conn_vox8_part;
    std::vector<long> nodemap_local2global;  // map[localnodeid]->globalnodeid
    std::vector<int> nodemap_global2local(
        n, -1);  // map[globalnodeid]->localnodeid
    std::vector<int> mpimap(
        nproc);  // map[iproc]->1 if neighbor, map[iproc]->0 if not neighbor

    std::vector<int> neighborrankid;  // list of neighboring processes
    std::vector<int>
        mpinode_pointer;  // pointer to list of local nodes to be communicated
    std::vector<int> mpinode_index;  // list of local nodes to be communicated
    //  mpinode_index[neighborrank]->start of nodes to be communicated in
    //  mpinode_index for neighborrank,
    std::vector<long> mpinode_all;
    std::vector<long> tmpmpinode_pointer;
    std::vector<int> surfnodes;  // list of local nodes on surface

    // for each partition
    int logskip = 100;
    double t2, t3, t4;
    t2 = omp_get_wtime();
    std::cout << "Start partition loop\n";
    for (long iproc = 0; iproc < nproc; ++iproc) {
        long n_part;
        long itmp, jtmp, ktmp, ltmp;
        long ne_vox8_part(partition_map_vox8[iproc].size());
        long ne_tet4_part(partition_map_tet4[iproc].size());
        long ielem, inode, glbelemid, glbnodeid;
        long tmpproc, totcommnodes;
        t3 = omp_get_wtime();

        n_part = 0;  // number of nodes in partition
        nodemap_local2global.clear();

        // set flag for nodes inside domain
        //-- for nodes in voxel elements
        for (ielem = 0; ielem < ne_vox8_part; ++ielem) {
            glbelemid = partition_map_vox8[iproc][ielem];
            for (inode = 0; inode < 8; ++inode) {
                glbnodeid = conn_vox8[9 * glbelemid + inode];
                if (nodemap_global2local[glbnodeid] == -1) {
                    nodemap_global2local[glbnodeid] = n_part;
                    nodemap_local2global.push_back(glbnodeid);
                    n_part++;
                }
            }
        }
        //-- for nodes in tetra elements
        for (ielem = 0; ielem < ne_tet4_part; ++ielem) {
            glbelemid = partition_map_tet4[iproc][ielem];
            for (inode = 0; inode < 4; ++inode) {
                glbnodeid = conn_tet4[5 * glbelemid + inode];
                if (nodemap_global2local[glbnodeid] == -1) {
                    nodemap_global2local[glbnodeid] = n_part;
                    nodemap_local2global.push_back(glbnodeid);
                    n_part++;
                }
            }
        }
        if (iproc % logskip == 0)
            std::cout << " part " << iproc << "\n --- setnodestook "
                      << omp_get_wtime() - t3 << "\n";
        t4 = omp_get_wtime();

        // set coor_part
        coor_part.resize(3 * n_part);
        for (inode = 0; inode < n_part; ++inode) {
            glbnodeid = nodemap_local2global[inode];
#ifdef USE_PART
//      std::cout << iproc << " " << inode << " " << glbnodeid << "\n";
//      if(glbnodeid < part_n)
//        psizelocal[2*iproc] = inode+1;
#endif
            coor_part[3 * inode] = coor[3 * glbnodeid];
            coor_part[3 * inode + 1] = coor[3 * glbnodeid + 1];
            coor_part[3 * inode + 2] = coor[3 * glbnodeid + 2];
        }
        if (iproc % logskip == 0)
            std::cout << " --- setcoortook " << omp_get_wtime() - t4 << "\n";
        t4 = omp_get_wtime();

        // set conn_vox8_part
        conn_vox8_part.resize(9 * ne_vox8_part);
        for (ielem = 0; ielem < ne_vox8_part; ++ielem) {
            glbelemid = partition_map_vox8[iproc][ielem];
            for (inode = 0; inode < 8; ++inode)
                conn_vox8_part[9 * ielem + inode] =
                    nodemap_global2local[conn_vox8[9 * glbelemid + inode]];
            conn_vox8_part[9 * ielem + 8] = conn_vox8[9 * glbelemid + 8];
        }
        if (iproc % logskip == 0)
            std::cout << " --- setconnvoxtook " << omp_get_wtime() - t4 << "\n";
        t4 = omp_get_wtime();

        // set conn_tet4_part
        conn_tet4_part.resize(5 * ne_tet4_part);
        for (ielem = 0; ielem < ne_tet4_part; ++ielem) {
            glbelemid = partition_map_tet4[iproc][ielem];
            for (inode = 0; inode < 4; ++inode)
                conn_tet4_part[5 * ielem + inode] =
                    nodemap_global2local[conn_tet4[5 * glbelemid + inode]];
            conn_tet4_part[5 * ielem + 4] = conn_tet4[5 * glbelemid + 4];
#ifdef USE_PART
            if (glbelemid < part_ne_tet4) {
                psizelocal[2 * iproc + 1] = ielem + 1;
                for (inode = 0; inode < 4; ++inode)
                    psizelocal[2 * iproc] =
                        std::max<long>(psizelocal[2 * iproc],
                                       conn_tet4_part[5 * ielem + inode]);
            }
#endif
        }
        if (iproc % logskip == 0)
            std::cout << " --- setconntettook " << omp_get_wtime() - t4 << "\n";
        t4 = omp_get_wtime();

        // set mpi nodes
        mpimap.clear();
        mpimap.resize(nproc);
        neighborrankid.clear();
        mpinode_pointer.clear();
        tmpmpinode_pointer.clear();
        mpinode_index.clear();
        mpinode_all.clear();
        for (inode = 0; inode < n_part; ++inode) {
            glbnodeid = nodemap_local2global[inode];
            if (node_duplicate_map[glbnodeid].size() != 1) {
                mpinode_all.push_back(inode);
                for (itmp = 0; itmp < node_duplicate_map[glbnodeid].size();
                     ++itmp)
                    mpimap[node_duplicate_map[glbnodeid][itmp]]++;
            }
        }
        if (iproc % logskip == 0)
            std::cout << " --- setmpinodes_part1_took " << omp_get_wtime() - t4
                      << "\n";
        t4 = omp_get_wtime();
        totcommnodes = 0;
        mpimap[iproc] = 0;
        for (itmp = 0; itmp < nproc; ++itmp) {
            if (mpimap[itmp] != 0) {
                mpinode_pointer.push_back(totcommnodes);
                tmpmpinode_pointer.push_back(totcommnodes);
                neighborrankid.push_back(itmp);
                totcommnodes += mpimap[itmp];
            }
        }
        mpinode_pointer.push_back(totcommnodes);
        mpinode_index.resize(totcommnodes);
        if (iproc % logskip == 0)
            std::cout << " --- setmpinodes_part2_took " << omp_get_wtime() - t4
                      << "\n";
        t4 = omp_get_wtime();
        for (ltmp = 0; ltmp < mpinode_all.size(); ++ltmp) {
            inode = mpinode_all[ltmp];
            glbnodeid = nodemap_local2global[inode];
            for (itmp = 0; itmp < node_duplicate_map[glbnodeid].size();
                 ++itmp) {
                tmpproc = node_duplicate_map[glbnodeid][itmp];
                if (tmpproc == iproc) continue;

                // search for tmpproc
                for (jtmp = 0; jtmp < neighborrankid.size(); ++jtmp) {
                    if (neighborrankid[jtmp] == tmpproc) {
                        ktmp = jtmp;
                        break;
                    }
                }
                mpinode_index[tmpmpinode_pointer[ktmp]] = inode;
                tmpmpinode_pointer[ktmp]++;
            }
        }
        if (iproc % logskip == 0)
            std::cout << " --- setmpinodes_part3_took " << omp_get_wtime() - t4
                      << "\n";
        t4 = omp_get_wtime();

        // set surface nodes
        surfnodes.clear();
        for (inode = 0; inode < n_part; ++inode) {
            glbnodeid = nodemap_local2global[inode];
            if (nodemap_global2surf[glbnodeid] != -1)
                surfnodes.push_back(inode);
        }
        if (iproc % logskip == 0)
            std::cout << " --- setsurfacenodestook " << omp_get_wtime() - t4
                      << "\n";
        t4 = omp_get_wtime();

        // write data
        write_data_hdf(coor_part, conn_tet4_part, conn_vox8_part,
                       neighborrankid, mpinode_pointer, mpinode_index,
                       surfnodes,
                       "./hdata/" + ToString(iproc, 6, '0') + ".data.h5");
        if (iproc % logskip == 0)
            std::cout << " --- writedatatook " << omp_get_wtime() - t4 << "\n";
        t4 = omp_get_wtime();

#ifdef OUTPUT_VTK_LOCAL
        std::vector<int> node_intvec1(n_part, 0);
        std::vector<double> node_doublevec(3 * n_part, 0.0);
        for (itmp = 0; itmp < mpinode_all.size(); ++itmp)
            node_intvec1[mpinode_all[itmp]] = itmp;
        write_vtk_hybrid_mesh(
            coor_part, conn_tet4_part, conn_vox8_part, node_intvec1,
            node_doublevec,
            "./hdata/" + ToString(iproc, 6, '0') + ".part_mpi.vtk");

        std::vector<int> node_intvec2(n_part, 0);
        for (itmp = 0; itmp < surfnodes.size(); ++itmp)
            node_intvec2[surfnodes[itmp]] = 1;
        write_vtk_hybrid_mesh(
            coor_part, conn_tet4_part, conn_vox8_part, node_intvec2,
            node_doublevec,
            "./hdata/" + ToString(iproc, 6, '0') + ".part_surf.vtk");
        if (iproc % logskip == 0)
            std::cout << " --- writevtktook " << omp_get_wtime() - t4 << "\n";
        t4 = omp_get_wtime();
#endif

        // clear flag
        for (ielem = 0; ielem < ne_vox8_part; ++ielem) {
            glbelemid = partition_map_vox8[iproc][ielem];
            for (inode = 0; inode < 8; ++inode)
                nodemap_global2local[conn_vox8[9 * glbelemid + inode]] = -1;
        }
        for (ielem = 0; ielem < ne_tet4_part; ++ielem) {
            glbelemid = partition_map_tet4[iproc][ielem];
            for (inode = 0; inode < 4; ++inode)
                nodemap_global2local[conn_tet4[5 * glbelemid + inode]] = -1;
        }
        if (iproc % logskip == 0) {
            std::cout << " --- clearflagtook " << omp_get_wtime() - t4 << "\n";
            std::cout << " --- partition total " << omp_get_wtime() - t3
                      << "\n";
        }
    }
#ifdef USE_PART
    std::ofstream fout;
    check_open("data/partsize.dat", fout);
    fout << part_n << " " << part_ne_tet4 << "\n";
    for (long iproc = 0; iproc < nproc; ++iproc)
        fout << psizelocal[2 * iproc] << " " << psizelocal[2 * iproc + 1]
             << "\n";
    fout.close();
#endif
    std::cout << "partition loop took " << omp_get_wtime() - t2 << "\n";
    std::cout << "TOTAL_ELAPSED_TIME " << omp_get_wtime() - t0 << "\n";
    return 0;
}
