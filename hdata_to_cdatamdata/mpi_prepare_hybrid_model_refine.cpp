// by Kohei Fujita 2015/2/11

// Cross link with 4to10_grid.f, MPI_ABC2_mod_lib.f, MPI_ABC4_mod.f, MPI_mapping.f, read_tet10_geometry_hdf.f
// Note: Use 64 bit integer mode (i.e., "integer" in Fortran will be equivalent as "long" in c/c++, when otherwise specified) when compiling Fortran functions

// Execute with same number of mpi processes as analysis model
// (one thread is needed per mpi process)

// --- Inputs --- <R0> indicate I/O by rank #0, <R%> indicate I/O by all processes
//  <R0> data/para_setting.dat -- read number of MPIprocs and OpenMPthreads
//  <R0> data/modeling_setting.dat -- read grid size ds
//  <R%> data/modeldomain.dat -- read size of whole domain
//  <R%> data/material.dat -- read size of whole domain
//  <R%> hdata/<000000>.data.h5 -- partitioned model

// --- Outputs --- <R0> indicate I/O by rank #0, <R%> indicate I/O by all processes
//  <R0> data/setting.dat -- total number of nodes and elements in fine mesh
//  <R0> data/2Doutput.dat -- total number of surface nodes to be outputted during analysis (includes nodes on micro domain)
//  <R%> mdata/<000000>.2Doutputflag.dat -- number of surface nodes in each partition, local node number of surface nodes in fine mesh, and their coordinates
//  <R%> mdata/<000000>.data.h5 -- partitioned fine model with boundary conditions, communication table, and output node information
//  <R%> cdata/<000000>.data.h5 -- partitioned coarse model with boundary conditions, communication table, and mapping between fine mesh
//  <R%> hdata/<000000>.part_quad_surf.vtk -- partitioned fine model (optional, set by OUTPUT_VTK)

// --- Options ---
// output vtk file of quadratic mesh (optional)
// #define OUTPUT_VTK
#define IAMINTEL

// --- Options ---
// number of refinement of mesh
// one refinement will partition one voxel element into 8 voxel elements and one tetrahedral element to 8 tetrahedral elements
// 1: no refinement
#define REFINE_FACTOR 1

// #define RENUMBER_ELEMENTS
// #define RENUMBER_NODES

// #define MAKE_INFINITE_BOUNDARY
//  Switch for post processing (use one mpi process)
// #define MERGE_SURFACE_VTK
//    when MERGE_SURFACE_VTK is chosen:
//    Set NROUGH: number of nodes to stride in output surface mesh
//    For example, if NROUGH=1, all nodes on surface will be outputted (spacing ds/2, as quadratic mesh is used)
//    If NROUGH=4, nodes on surface with spacing of 2*ds will be outputted.
#define NROUGH 1
//   Outputs are
//    vtk_surface.geom -- Surface mesh
//    vtk_surface.rank -- Rank number of surface. To make vtk file, $cat vtk_surface.geom vtk_surface.rank > vtk_surface_rank.vtk
//    vtk_surface.displshead -- Header for making output file, $cat vtk_surface.geom vtk_surface.displshead <displs_data> > vtk_surface_displs.vtk
//    Here, <displs_data> is displacement data in big endian floats, in the order designated by mdata/<000000>.2Doutputflag.dat
//    displsx_node0_partition0,displsy_node0_partition0,
//    displsz_node0_partition0,displsx_node1_partition0,...,
//    displsz_node(n-1)_partition(m-1)

// #define OUTPUT_3DRESULT
#define USECRS4
#define USECRS10

// #define STRUCTURED_MESH_SET_MATERIAL

//-----------------------------------------------------------
// Other debug flags
// #define DEBUG_NO_MPI
// #define DEBUG_NO_HDF5
// #define DEBUG_NO_COMM

//******************************************************************************
// include files
//
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <cstdlib>
#include <map>
#include <climits>

#ifndef DEBUG_NO_MPI
#include "mpi.h"
#endif

#ifndef DEBUG_NO_HDF5
extern "C"
{
#include "hdf_lib.c"
}
#endif

//******************************************************************************
// macros
//
#ifndef Swap2Bytes
#define Swap2Bytes(data) ((((data) >> 8) & 0x00FF) | (((data) << 8) & 0xFF00))
#endif
#ifndef SWAP_SHORT
#define SWAP_SHORT(x) (*(unsigned short *)&(x) = Swap4Bytes(*(unsigned short *)&(x)))
#endif
#ifndef Swap4Bytes
#define Swap4Bytes(data) ((((data) >> 24) & 0x000000FF) | (((data) >> 8) & 0x0000FF00) | \
                          (((data) << 8) & 0x00FF0000) | (((data) << 24) & 0xFF000000))
#endif
#ifndef SWAP_INT
#define SWAP_INT(x) (*(unsigned int *)&(x) = Swap4Bytes(*(unsigned int *)&(x)))
#endif

#ifndef Swap8Bytes
#define Swap8Bytes(data) ((((data) >> 56) & 0x00000000000000FF) | (((data) >> 40) & 0x000000000000FF00) | \
                          (((data) >> 24) & 0x0000000000FF0000) | (((data) >> 8) & 0x00000000FF000000) |  \
                          (((data) << 8) & 0x000000FF00000000) | (((data) << 24) & 0x0000FF0000000000) |  \
                          (((data) << 40) & 0x00FF000000000000) | (((data) << 56) & 0xFF00000000000000))
#endif
#ifndef SWAP_DOUBLE
#define SWAP_DOUBLE(x) (*(unsigned long *)&(x) = Swap8Bytes(*(unsigned long *)&(x)))
#endif
#ifndef SWAP_FLOAT
#define SWAP_FLOAT(x) (*(unsigned int *)&(x) = Swap4Bytes(*(unsigned int *)&(x)))
#endif

//******************************************************************************
// fortan libraries
//
extern "C"
{
  void set_gridsize_(long *n4, long *coorsize, double *coor4, double *dsg, long *nxg, long *nyg, long *nzg, double *co);
  void convert_4to10_grid_(long *n4, long *ne, long *coorsize, double *coor, long *conn, long *nl, double *dsg,
                           long *nxg, long *nyg, long *nzg, double *co, long *n10, long *flag);
  void convert_4to10_grid_merge4and10_(long *n4, long *ne, long *coorsize, double *coor, long *conn, long *nl, double *dsg,
                                       long *nxg, long *nyg, long *nzg, double *co, long *n10, long *flag);
  void convert_4to10_confirming_mesh_(long *n4, long *ne, long *coorsize, double *coor, long *conn, long *n10);
  void abc4_lib_(int *n_linear, int *ne, int *nmat, int *ib_linear, long *fidc);
  void abc2_lib2_(long *n, long *ne, long *nmat, long *ib_quad, long *fidm, int *conn, long *connbufsize,
                  double *coor, long *coorbufsize, double *rmat, int *ibc, long *ibcsize,
                  long *bufsl, int *buf10, double *buf11, double *buf12, double *buf13, int *buf14);
  void mapping_pairs_lib_(int *n_quad, int *ne, int *n_linear, int *ne_linear, long *fidm, long *fidc);
  void abc2_lib3_(long *n_long, long *ne_long, long *nma_long, long *fidm, long *bufsl);
#ifdef USECRS4
  void crs_nodes4_lib_(int *n, int *ne, long *fid);
  void crs_connect_block4_lib_(int *n, int *ne, int *ib, long *fid);
#endif
#ifdef USECRS10
  void crs_nodes_lib_(int *n, int *ne, long *fid);
  void crs_connect_block_lib_(int *n, int *ne, int *ib, long *fid);
#endif
  void modify_mat(int rank, int nummat, const std::vector<double> &coor4, std::vector<int> &conn_tet4, std::vector<int> &conn_vox8, std::vector<int> &conn_tet10);
}
void convert_infinite_boundary_coor(int myrank, double dsorg, double gmaxx, double gmaxy, const std::vector<int> &conn_alltet10, const std::vector<int> &nodemap_linear2quad, std::vector<double> &coor_quad, std::vector<double> &coor);

//******************************************************************************
// templates
//
template <class T>
void swap_bytes(T &in) {}
template <>
void swap_bytes<int>(int &in) { SWAP_INT(in); }
template <>
void swap_bytes<float>(float &in) { SWAP_FLOAT(in); }
template <class T>
void vec_swap_bytes(std::vector<T> &in)
{
  for (long i = 0; i < in.size(); ++i)
    swap_bytes<T>(in[i]);
}
template <class T>
std::string ToString(const T in)
{
  std::stringstream ss;
  ss << in;
  return ss.str();
}
template <class T>
std::string ToString(const T in, const long size, const char filling)
{
  std::string strtmp(ToString<T>(in));
  std::string strout(size, filling);
  const long strsize(strtmp.size());
  if (strsize > size)
    return strtmp;
  else
  {
    for (long i = 0; i < strsize; ++i)
      strout[i + size - strsize] = strtmp[i];
    return strout;
  }
}
template <class T>
void Error(const T message)
{
  std::cout << "ERROR: " << message << std::endl;
  exit(-1);
}
template <class T>
void check_open(const std::string filename, T &file)
{
  file.open(filename.c_str());
  if (!file)
    Error("cannot find " + filename);
}

//******************************************************************************
// functions
//
void ToFortranStyle(const std::vector<int> &in, std::vector<int> &out)
{
  out.resize(in.size());
  for (long i = 0; i < in.size(); ++i)
    out[i] = in[i] + 1;
}

void read_global_settings(const int myrank, double &ds, double &gminx, double &gminy, double &gminz,
                          double &gmaxx, double &gmaxy, int &nmat, int &nproc, int &nthread)
{
  double dvec[6];
  if (myrank == 0)
  {
    std::string tmp;
    std::ifstream fin;
    check_open("./data/modeldomain.dat", fin);
    fin >> tmp; // xmin-plane
    fin >> dvec[0];
    fin >> tmp; // xmax-plane
    fin >> dvec[1];
    fin >> tmp; // ymin-plane
    fin >> dvec[2];
    fin >> tmp; // ymax-plane
    fin >> dvec[3];
    fin >> tmp; // zmin-plane
    fin >> dvec[4];
    fin.close();
    check_open("./data/modeling_setting.dat", fin);
    fin >> tmp; // nx,
    fin >> tmp; // ny
    fin >> tmp; //
    fin >> tmp; //
    fin >> tmp; // ds
    fin >> dvec[5];
    fin >> tmp; // num
    fin >> tmp; // of
    fin >> tmp; // layer
    fin >> nmat;
    fin.close();
    check_open("data/para_setting.dat", fin);
    fin >> tmp; // num
    fin >> tmp; // of
    fin >> tmp; // MPI
    fin >> tmp; // thread
    fin >> nproc;
    fin >> tmp; // num
    fin >> tmp; // of
    fin >> tmp; // OpenMP
    fin >> tmp; // thread
    fin >> nthread;
    fin.close();
  }
#ifndef DEBUG_NO_MPI
  MPI_Bcast(dvec, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nmat, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nthread, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  gminx = dvec[0];
  gminy = dvec[2];
  gmaxx = dvec[1];
  gmaxy = dvec[3];
  gminz = dvec[4];
  ds = dvec[5];
}

void read_local_data(std::string name, int &n, int &ne_tet4, int &ne_vox8, int &neighbors,
                     std::vector<double> &coor, std::vector<int> &conn_tet4, std::vector<int> &conn_vox8,
                     std::vector<int> &neighborrankid, std::vector<int> &mpinode_pointer,
                     std::vector<int> &mpinode_index, std::vector<int> &surfnodes,
                     std::vector<int> &global_node_id, std::vector<long> &global_node_id64)
{
#ifndef DEBUG_NO_HDF5
  long fid(0), count, size;
  int settings[7];

  hdf_open_file_(name.c_str(), &fid);

  size = 7;
  hdf_read_int_array_(&fid, "/settings", settings, &size, &count);
  n = settings[0]; // = coor_part.size()/3; // number of nodes
  coor.resize(3 * n);
  ne_vox8 = settings[1]; // = conn_vox8_part.size()/9; // number of vox8 elements
  conn_vox8.resize(9 * ne_vox8);
  ne_tet4 = settings[2]; //= conn_tet4_part.size()/5; // number of tet4 elements
  conn_tet4.resize(5 * ne_tet4);
  surfnodes.resize(settings[3]); // = surfnodes.size(); // number of surface nodes
  neighbors = settings[4];       //= neighborrankid.size(); // number of neighboring ranks
  neighborrankid.resize(neighbors);
  mpinode_pointer.resize(settings[5]); // = mpinode_pointer.size(); // number of neighboring ranks
  mpinode_index.resize(settings[6]);   // = mpinode_index.size(); // number of nodes to send/recv

  size = coor.size();
  hdf_read_double_array_(&fid, "/coor", &(coor[0]), &size, &count);

  size = conn_vox8.size();
  if (size != 0)
    hdf_read_int_array_(&fid, "/conn_vox8", &(conn_vox8[0]), &size, &count);

  size = conn_tet4.size();
  if (size != 0)
    hdf_read_int_array_(&fid, "/conn_tet4", &(conn_tet4[0]), &size, &count);

  size = surfnodes.size();
  if (size != 0)
    hdf_read_int_array_(&fid, "/surfnodeid", &(surfnodes[0]), &size, &count);

  size = neighborrankid.size();
  if (size != 0)
    hdf_read_int_array_(&fid, "/neighborrankid", &(neighborrankid[0]), &size, &count);

  size = mpinode_pointer.size();
  if (size != 0)
    hdf_read_int_array_(&fid, "/mpinode_pointer", &(mpinode_pointer[0]), &size, &count);

  size = mpinode_index.size();
  if (size != 0)
    hdf_read_int_array_(&fid, "/mpinode_index", &(mpinode_index[0]), &size, &count);

  // optional global node id (1-based)
  global_node_id.clear();
  global_node_id64.clear();
  if (hdf_has_dataset_(&fid, (char *)"/global_node_id")) {
    size = 0;
    hdf_get_array_size_(&fid, (char *)"/global_node_id", &size);
    if (size > 0) {
      global_node_id.resize(size);
      hdf_read_int_array_(&fid, (char *)"/global_node_id",
                          &(global_node_id[0]), &size, &count);
    }
  }
  if (hdf_has_dataset_(&fid, (char *)"/global_node_id64")) {
    size = 0;
    hdf_get_array_size_(&fid, (char *)"/global_node_id64", &size);
    if (size > 0) {
      global_node_id64.resize(size);
      hdf_read_long_array_(&fid, (char *)"/global_node_id64",
                           &(global_node_id64[0]), &size, &count);
    }
  }

  hdf_close_file_(&fid);

#endif
}

void convert_vox8tet4_to_tet4(const std::vector<int> &conn_vox8, const std::vector<int> &conn_tet4,
                              std::vector<int> &conn_alltet4)
{
  int ne_vox8(conn_vox8.size() / 9);
  int ne_tet4(conn_tet4.size() / 5);
  int ne(6 * ne_vox8 + ne_tet4);
  int connlist[6][4];
  // set conn list
  connlist[0][0] = 2;
  connlist[0][1] = 6;
  connlist[0][2] = 4;
  connlist[0][3] = 5;
  connlist[1][0] = 2;
  connlist[1][1] = 0;
  connlist[1][2] = 4;
  connlist[1][3] = 3;
  connlist[2][0] = 2;
  connlist[2][1] = 1;
  connlist[2][2] = 4;
  connlist[2][3] = 0;
  connlist[3][0] = 2;
  connlist[3][1] = 5;
  connlist[3][2] = 4;
  connlist[3][3] = 1;
  connlist[4][0] = 2;
  connlist[4][1] = 3;
  connlist[4][2] = 4;
  connlist[4][3] = 6;
  connlist[5][0] = 3;
  connlist[5][1] = 4;
  connlist[5][2] = 6;
  connlist[5][3] = 7;
  // HERE
  connlist[0][0] = 1;
  connlist[0][1] = 7;
  connlist[0][2] = 3;
  connlist[0][3] = 4;

  connlist[1][0] = 2;
  connlist[1][1] = 5;
  connlist[1][2] = 6;
  connlist[1][3] = 7;

  connlist[2][0] = 1;
  connlist[2][1] = 7;
  connlist[2][2] = 2;
  connlist[2][3] = 3;

  connlist[3][0] = 3;
  connlist[3][1] = 1;
  connlist[3][2] = 4;
  connlist[3][3] = 0;

  connlist[4][0] = 2;
  connlist[4][1] = 5;
  connlist[4][2] = 7;
  connlist[4][3] = 1;

  connlist[5][0] = 5;
  connlist[5][1] = 7;
  connlist[5][2] = 1;
  connlist[5][3] = 4;

  conn_alltet4.resize(5 * ne);

  int i, j, k, ns(5);
  for (i = 0; i < ne_vox8; ++i)
  {
    for (j = 0; j < 6; ++j) // elements in voxel
    {
      for (k = 0; k < 4; ++k) // nodes of tet element
        conn_alltet4[ns * 6 * i + ns * j + k] = conn_vox8[9 * i + connlist[j][k]];
      conn_alltet4[ns * 6 * i + ns * j + 4] = conn_vox8[9 * i + 8];
    }
  }
  k = ns * ne_vox8 * 6;
  for (i = 0; i < ne_tet4; ++i)
  {
    conn_alltet4[k + ns * i + 0] = conn_tet4[5 * i + 0];
    conn_alltet4[k + ns * i + 1] = conn_tet4[5 * i + 1];
    conn_alltet4[k + ns * i + 2] = conn_tet4[5 * i + 2];
    conn_alltet4[k + ns * i + 3] = conn_tet4[5 * i + 3];
    conn_alltet4[k + ns * i + 4] = conn_tet4[5 * i + 4];
  }
}

void convert_linear2quad(
    int myrank, double dsorg, double gmaxx, double gmaxy,
    bool merge4and10, const std::vector<double> &coor, const std::vector<int> &conn_vox8,
    const std::vector<int> &conn_tet4, double dsg, std::vector<double> &coor_quad, std::vector<int> &conn_alltet10,
    std::vector<int> &nodemap_quad2linear, std::vector<int> &nodemap_linear2quad, std::vector<int> &nodetype)
{
  long n4, coorsize, nxg, nyg, nzg, ne;
  double co[3];
  std::vector<int> conn_alltet4;
  std::vector<double> coor_quad_tmp;
  std::vector<long> conn_alltet10_tmp;

  convert_vox8tet4_to_tet4(conn_vox8, conn_tet4, conn_alltet4);
  ne = conn_alltet4.size() / 5;
  conn_alltet10_tmp.resize(11 * ne);
  for (long ie = 0; ie < ne; ++ie)
  {
    conn_alltet10_tmp[11 * ie + 0] = conn_alltet4[5 * ie + 0] + 1; // to Fortran style
    conn_alltet10_tmp[11 * ie + 1] = conn_alltet4[5 * ie + 1] + 1;
    conn_alltet10_tmp[11 * ie + 2] = conn_alltet4[5 * ie + 2] + 1;
    conn_alltet10_tmp[11 * ie + 3] = conn_alltet4[5 * ie + 3] + 1;
    conn_alltet10_tmp[11 * ie + 10] = conn_alltet4[5 * ie + 4] + 1;
  }

  //  long n10(ne*15); // can be reduced
  long n10(ne * 4); // can be reduced
  coor_quad_tmp.resize(3 * n10);
  n4 = coor.size() / 3;
  coorsize = n10;

  for (long i = 0; i < 3 * n4; ++i)
    coor_quad_tmp[i] = coor[i];

  set_gridsize_(&n4, &coorsize, &(coor_quad_tmp[0]), &dsg, &nxg, &nyg, &nzg, co);
  //  for( long nl = 30, flag = 0; flag == 0; nl += 10 )
  for (long nl = 300, flag = 0; flag == 0; nl += 20)
  {

    convert_4to10_confirming_mesh_(&n4, &ne, &coorsize, &(coor_quad_tmp[0]), &(conn_alltet10_tmp[0]), &n10);
    flag = 1;

    /*
        if(merge4and10)
        {
          Error("do not use convert_4to10_grid_merge4and10_");
          convert_4to10_grid_merge4and10_(&n4,&ne,&coorsize,&(coor_quad_tmp[0]),&(conn_alltet10_tmp[0]),
            &nl,&dsg,&nxg,&nyg,&nzg,co,&n10,&flag);
        }
        else
        {
          convert_4to10_grid_(&n4,&ne,&coorsize,&(coor_quad_tmp[0]),&(conn_alltet10_tmp[0]),
            &nl,&dsg,&nxg,&nyg,&nzg,co,&n10,&flag);
        }
      */
  }

  coor_quad.resize(3 * n10);
  for (long i = 0; i < 3 * n10; ++i)
    coor_quad[i] = coor_quad_tmp[i];

  conn_alltet10.resize(ne * 11);
  for (long i = 0; i < ne * 11; ++i)
    conn_alltet10[i] = conn_alltet10_tmp[i] - 1; // to C style

  nodemap_linear2quad.resize(n4);
  nodemap_quad2linear.resize(n10);

  for (long i = 0; i < n4; ++i)
  {
    nodemap_linear2quad[i] = i;
    nodemap_quad2linear[i] = i;
  }
  for (long i = n4; i < n10; ++i)
    nodemap_quad2linear[i] = -1;

  if (n4 != nodetype.size())
    Error("error in check nodetype");
  for (long i = 0; i < n4; ++i)
    nodetype[i]++;
  for (long i = 0; i < n10 - n4; ++i)
    nodetype.push_back(0);

//***
#ifdef MAKE_INFINITE_BOUNDARY
  std::vector<double> coor_tmp(coor);
  coor_quad_tmp.clear();
  coor_quad_tmp = coor_quad;
  convert_infinite_boundary_coor(myrank, dsorg, gmaxx, gmaxy, conn_alltet10, nodemap_linear2quad, coor_quad_tmp, coor_tmp);
  int j, index;
  for (long i = 0; i < ne; ++i)
  {
    for (j = 4; j < 10; ++j)
    {
      index = conn_alltet10[11 * i + j];
      coor_quad[3 * index + 0] = coor_quad_tmp[3 * index + 0];
      coor_quad[3 * index + 1] = coor_quad_tmp[3 * index + 1];
      coor_quad[3 * index + 2] = coor_quad_tmp[3 * index + 2];
    }
  }
#endif
  //***
}

template <class T1, class T2>
void write_vtk_hybrid_mesh(const std::vector<double> &coor_in, const std::vector<T1> &intconn_tet4,
                           const std::vector<T1> &intconn_vox8, const std::vector<T2> &node_intscalar, const std::vector<double> &node_doublevec,
                           std::string name)
{
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

  for (i = 0; i < nodes; ++i)
  {
    coor[3 * i] = coor_in[3 * i];
    coor[3 * i + 1] = coor_in[3 * i + 1];
    coor[3 * i + 2] = coor_in[3 * i + 2];
    intscalar[i] = node_intscalar[i];
    floatvec[3 * i] = node_doublevec[3 * i];
    floatvec[3 * i + 1] = node_doublevec[3 * i + 1];
    floatvec[3 * i + 2] = node_doublevec[3 * i + 2];
  }

  for (i = 0; i < elements_vox8; ++i)
  {
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
    inttype[i] = 12; // 12 for hexahedral elements
  }
  for (i = 0; i < elements_tet4; ++i)
  {
    conn[9 * elements_vox8 + 5 * i] = 4;
    for (j = 0; j < 4; ++j)
      conn[9 * elements_vox8 + 5 * i + j + 1] = intconn_tet4[5 * i + j];
    matid[elements_vox8 + i] = intconn_tet4[5 * i + 4];
    inttype[elements_vox8 + i] = 10; // 10 for tet4 elements
  }

#ifdef IAMINTEL
  vec_swap_bytes(coor);
  vec_swap_bytes(conn);
  vec_swap_bytes(intscalar);
  vec_swap_bytes(floatvec);
  vec_swap_bytes(matid);
  vec_swap_bytes(inttype);
#endif
  std::ofstream fout;
  check_open(name, fout);
  fout << "# vtk DataFile Version 2.0\n";
  fout << "file: " << name << std::endl;
  fout << "BINARY\n";
  fout << "DATASET UNSTRUCTURED_GRID\n";
  fout << "POINTS " << nodes << " float\n";
  fout.write((char *)&(coor[0]), sizeof(float) * 3 * nodes);
  fout << "CELLS " << elements << " " << elementcomponents << std::endl;
  fout.write((char *)&(conn[0]), sizeof(int) * elementcomponents);
  fout << "CELL_TYPES " << elements << std::endl;
  fout.write((char *)&(inttype[0]), sizeof(int) * elements);
  fout << "CELL_DATA " << elements << std::endl;
  fout << "SCALARS material int 1\n";
  fout << "LOOKUP_TABLE default\n";
  fout.write((char *)&(matid[0]), sizeof(int) * elements);
  fout << "POINT_DATA " << nodes << "\n";
  fout << "VECTORS floatvec float\n";
  fout.write((char *)&(floatvec[0]), sizeof(float) * 3 * nodes);
  fout << "SCALARS intvec int 1\n";
  fout << "LOOKUP_TABLE default\n";
  fout.write((char *)&(intscalar[0]), sizeof(float) * nodes);
  fout.close();
}

void write_vtk_tet10_mesh(const std::vector<double> &coor_in, const std::vector<int> &intconn_tet10,
                          const std::vector<int> &node_intscalar, const std::vector<double> &node_doublevec,
                          std::string name)
{
  long nodes(coor_in.size() / 3);
  long elements(intconn_tet10.size() / 11);
  long elementcomponents(intconn_tet10.size());
  long i, j;

  std::vector<float> coor(3 * nodes);
  std::vector<int> conn(elementcomponents);
  std::vector<int> inttype(elements);
  std::vector<int> matid(elements);
  std::vector<int> intscalar(nodes);
  std::vector<float> floatvec(3 * nodes);

  for (i = 0; i < nodes; ++i)
  {
    coor[3 * i] = coor_in[3 * i];
    coor[3 * i + 1] = coor_in[3 * i + 1];
    coor[3 * i + 2] = coor_in[3 * i + 2];
    intscalar[i] = node_intscalar[i];
    floatvec[3 * i] = node_doublevec[3 * i];
    floatvec[3 * i + 1] = node_doublevec[3 * i + 1];
    floatvec[3 * i + 2] = node_doublevec[3 * i + 2];
  }

  for (i = 0; i < elements; ++i)
  {
    conn[11 * i] = 10;
    for (j = 0; j < 10; ++j)
      conn[11 * i + 1 + j] = intconn_tet10[11 * i + j];
    matid[i] = intconn_tet10[11 * i + 10];
    inttype[i] = 24; // 24 for tet10 elements
  }
#ifdef IAMINTEL
  vec_swap_bytes(coor);
  vec_swap_bytes(conn);
  vec_swap_bytes(intscalar);
  vec_swap_bytes(floatvec);
  vec_swap_bytes(matid);
  vec_swap_bytes(inttype);
#endif
  std::ofstream fout;
  check_open(name, fout);
  fout << "# vtk DataFile Version 2.0\n";
  fout << "file: " << name << std::endl;
  fout << "BINARY\n";
  fout << "DATASET UNSTRUCTURED_GRID\n";
  fout << "POINTS " << nodes << " float\n";
  fout.write((char *)&(coor[0]), sizeof(float) * 3 * nodes);
  fout << "CELLS " << elements << " " << elementcomponents << std::endl;
  fout.write((char *)&(conn[0]), sizeof(int) * elementcomponents);
  fout << "CELL_TYPES " << elements << std::endl;
  fout.write((char *)&(inttype[0]), sizeof(int) * elements);
  fout << "CELL_DATA " << elements << std::endl;
  fout << "SCALARS material int 1\n";
  fout << "LOOKUP_TABLE default\n";
  fout.write((char *)&(matid[0]), sizeof(int) * elements);
  fout << "POINT_DATA " << nodes << "\n";
  fout << "VECTORS floatvec float\n";
  fout.write((char *)&(floatvec[0]), sizeof(float) * 3 * nodes);
  fout << "SCALARS intvec int 1\n";
  fout << "LOOKUP_TABLE default\n";
  fout.write((char *)&(intscalar[0]), sizeof(float) * nodes);
  fout.close();
}

void set_boundary(const std::vector<double> &coor, double ds,
                  double gminx, double gminy, double gminz, double gmaxx, double gmaxy, std::vector<int> &boundarynode)
{
  int n(coor.size() / 3);
  double eps(ds * 0.01);
  bool flag;
  boundarynode.clear();
  for (int i = 0; i < n; ++i)
  {
    flag = false;
#ifndef MAKE_INFINITE_BOUNDARY
    if (fabs(coor[3 * i] - gminx) < eps)
      flag = true;
    if (fabs(coor[3 * i + 1] - gminy) < eps)
      flag = true;
    if (fabs(coor[3 * i] - gmaxx) < eps)
      flag = true;
    if (fabs(coor[3 * i + 1] - gmaxy) < eps)
      flag = true;
#endif
    if (fabs(coor[3 * i + 2] - gminz) < eps)
      flag = true;
    if (flag)
      boundarynode.push_back(i);
  }
}

void MPI_check_comm_settings(const std::vector<double> &coor, double ds, const std::vector<int> &neighborrankid,
                             const std::vector<int> &mpinode_pointer, const std::vector<int> &mpinode_index,
                             int myrank, int totalrank, const std::vector<int> &nodetype)
{
#ifndef DEBUG_NO_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "START MPI CHECK\n";
  int nbrs(neighborrankid.size());
  std::vector<MPI_Request> req_send(nbrs);
  std::vector<MPI_Request> req_recv(nbrs);
  std::vector<MPI_Status> stat_send(nbrs);
  std::vector<MPI_Status> stat_recv(nbrs);
  std::vector<MPI_Request> reqi_send(nbrs);
  std::vector<MPI_Request> reqi_recv(nbrs);
  std::vector<MPI_Status> stati_send(nbrs);
  std::vector<MPI_Status> stati_recv(nbrs);
  std::vector<int> nodesize(nbrs);
  std::vector<int> nodesize_check(nbrs);

  int tagsend, tagrecv, nn, j;
  double eps(ds * 0.01);
  nn = 100000;

  // check number of nodes to send/recv
  for (int i = 0; i < nbrs; ++i)
  {
    nodesize[i] = mpinode_pointer[i + 1] - mpinode_pointer[i];
    tagsend = 1;
    (myrank + 1) * nn + neighborrankid[i]; // tag=(sendid+1)*nn+recvid
    MPI_Isend(&(nodesize[i]), 1, MPI_INT, neighborrankid[i],
              tagsend, MPI_COMM_WORLD, &(req_send[i]));
    tagrecv = 1;
    (neighborrankid[i] + 1) * nn + myrank; // tag=(sendid+1)*nn+recvid
    MPI_Irecv(&(nodesize_check[i]), 1, MPI_INT, neighborrankid[i],
              tagrecv, MPI_COMM_WORLD, &(req_recv[i]));
  }

  if (nbrs != 0)
  {
    MPI_Waitall(nbrs, &(req_recv[0]), &(stat_recv[0]));
    MPI_Waitall(nbrs, &(req_send[0]), &(stat_send[0]));
  }

  for (int i = 0; i < nbrs; ++i)
  {
    if (nodesize_check[i] != nodesize[i])
      Error("MPI_CHECK_FAILED: " + ToString(myrank) + " " + ToString(i) + " " +
            ToString(neighborrankid[i]) + " " + ToString(nodesize_check[i]) + " " + ToString(nodesize[i]));
  }

  std::vector<double> coorbuf(3 * mpinode_pointer[nbrs]);
  std::vector<double> coorbuf_check(3 * mpinode_pointer[nbrs]);
  std::vector<int> intbuf(mpinode_pointer[nbrs]);
  std::vector<int> intbuf_check(mpinode_pointer[nbrs]);
  // check coordinates of nodes to send/recv
  for (int i = 0; i < mpinode_pointer[nbrs]; ++i)
  {
    coorbuf[3 * i + 0] = coor[3 * mpinode_index[i] + 0];
    coorbuf[3 * i + 1] = coor[3 * mpinode_index[i] + 1];
    coorbuf[3 * i + 2] = coor[3 * mpinode_index[i] + 2];
    intbuf[i] = nodetype[mpinode_index[i]];
  }
  for (int i = 0; i < nbrs; ++i)
  {
    tagsend = 1;
    (myrank + 1) * nn + neighborrankid[i]; // tag=(sendid+1)*nn+recvid
    MPI_Isend(&(coorbuf[3 * mpinode_pointer[i]]), 3 * nodesize[i], MPI_DOUBLE, neighborrankid[i],
              tagsend, MPI_COMM_WORLD, &(req_send[i]));
    tagrecv = 1;
    (neighborrankid[i] + 1) * nn + myrank; // tag=(sendid+1)*nn+recvid
    MPI_Irecv(&(coorbuf_check[3 * mpinode_pointer[i]]), 3 * nodesize[i], MPI_DOUBLE, neighborrankid[i],
              tagrecv, MPI_COMM_WORLD, &(req_recv[i]));

    tagsend = 1;
    (myrank + 1) * nn + neighborrankid[i] + 99999; // tag=(sendid+1)*nn+recvid+99999
    MPI_Isend(&(intbuf[mpinode_pointer[i]]), nodesize[i], MPI_INT, neighborrankid[i],
              tagsend, MPI_COMM_WORLD, &(reqi_send[i]));
    tagrecv = 1;
    (neighborrankid[i] + 1) * nn + myrank + 99999; // tag=(sendid+1)*nn+recvid+99999
    MPI_Irecv(&(intbuf_check[mpinode_pointer[i]]), nodesize[i], MPI_INT, neighborrankid[i],
              tagrecv, MPI_COMM_WORLD, &(reqi_recv[i]));
  }
  if (nbrs != 0)
  {
    MPI_Waitall(nbrs, &(req_recv[0]), &(stat_recv[0]));
    MPI_Waitall(nbrs, &(req_send[0]), &(stat_send[0]));
    MPI_Waitall(nbrs, &(reqi_recv[0]), &(stati_recv[0]));
    MPI_Waitall(nbrs, &(reqi_send[0]), &(stati_send[0]));
  }
  for (int i = 0; i < nbrs; ++i)
  {
    for (j = mpinode_pointer[i]; j < mpinode_pointer[i + 1]; ++j)
    {
      if (fabs(coorbuf[3 * j] - coorbuf_check[3 * j]) > eps)
        Error("Error: MPI_CHECK_FAILED coor x: " + ToString(myrank) + " " + ToString(i) + " " +
              ToString(neighborrankid[i]) + " " + ToString(j) + " " + ToString(coorbuf[3 * j]) + " " + ToString(coorbuf_check[3 * j]));
      if (fabs(coorbuf[3 * j + 1] - coorbuf_check[3 * j + 1]) > eps)
        Error("Error: MPI_CHECK_FAILED coor y: " + ToString(myrank) + " " + ToString(i) + " " +
              ToString(neighborrankid[i]) + " " + ToString(j) + " " + ToString(coorbuf[3 * j + 1]) + " " + ToString(coorbuf_check[3 * j + 1]));
      if (fabs(coorbuf[3 * j + 2] - coorbuf_check[3 * j + 2]) > eps)
        Error("Error: MPI_CHECK_FAILED coor z: " + ToString(myrank) + " " + ToString(i) + " " +
              ToString(neighborrankid[i]) + " " + ToString(j) + " " + ToString(coorbuf[3 * j + 2]) + " " + ToString(coorbuf_check[3 * j + 2]));
      if (intbuf[j] != intbuf_check[j])
        Error("Error: MPI_CHECK_FAILED int: " + ToString(myrank) + " " + ToString(i) + " " +
              ToString(neighborrankid[i]) + " " + ToString(j) + " " + ToString(intbuf[j]) + " " + ToString(intbuf_check[j]));
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "MPI CHECK SUCCEED\n";
#endif
}

inline double dist2(double *x1, double *x2)
{
  return (x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[1] - x2[1]) * (x1[1] - x2[1]) + (x1[2] - x2[2]) * (x1[2] - x2[2]);
}

double global_coor_pair_eps;
class coor_pair
{
public:
  double x[3];
  int id;
  int type;
  coor_pair(int id_, double x0, double x1, double x2, int type_)
  {
    id = id_;
    x[0] = x0;
    x[1] = x1;
    x[2] = x2;
    type = type_;
  }
};

bool operator<(coor_pair a, coor_pair b)
{
  if (fabs(a.x[0] - b.x[0]) < global_coor_pair_eps)
  {
    if (fabs(a.x[1] - b.x[1]) < global_coor_pair_eps)
    {
      if (fabs(a.x[2] - b.x[2]) < global_coor_pair_eps)
        return a.type < b.type;
      else
        return a.x[2] < b.x[2];
    }
    else
      return a.x[1] < b.x[1];
  }
  else
    return a.x[0] < b.x[0];
}

void sort_by_coordinate(const std::vector<double> &coor, const std::vector<int> &nodetype, std::vector<int> &nodeids)
{
  std::vector<coor_pair> pairs;
  for (int i = 0; i < nodeids.size(); ++i)
    pairs.push_back(coor_pair(nodeids[i], coor[3 * nodeids[i]], coor[3 * nodeids[i] + 1], coor[3 * nodeids[i] + 2], nodetype[nodeids[i]]));
  std::sort(pairs.begin(), pairs.end());
  for (int i = 0; i < nodeids.size(); ++i)
    nodeids[i] = pairs[i].id;
}

void renumber_elements(const std::vector<int> &mpinode_linear_index, const std::vector<int> &nodemap_quad2linear, const std::vector<int> &nodemap_linear2quad,
                       std::vector<int> &conn_vox8, std::vector<int> &conn_tet4, std::vector<int> &conn_alltet10, int myrank)
{
  int netet(conn_tet4.size() / 5);
  int nevox(conn_vox8.size() / 9);
  std::vector<int> nodeflag(nodemap_linear2quad.size(), 0);
  std::vector<int> conn_vox8_new;
  std::vector<int> conn_tet4_new;
  std::vector<int> conn_vox27_new;
  std::vector<int> conn_tet10_new;

  for (int i = 0; i < mpinode_linear_index.size(); ++i)
    nodeflag[mpinode_linear_index[i]] = 1;

  int i, j, disp(6 * 11 * nevox);
  int nevox_outer(0);
  int netet_outer(0);
  bool flag;

  for (i = 0; i < nevox; ++i)
  {
    flag = false;
    for (j = 0; j < 8; ++j)
    {
      if (nodeflag[conn_vox8[9 * i + j]] == 1)
        flag = true;
    }
    if (flag)
    {
      nevox_outer++;
      for (j = 0; j < 9; ++j)
        conn_vox8_new.push_back(conn_vox8[9 * i + j]);
      for (j = 0; j < 66; ++j)
        conn_vox27_new.push_back(conn_alltet10[66 * i + j]);
    }
  }
  for (i = 0; i < netet; ++i)
  {
    flag = false;
    for (j = 0; j < 4; ++j)
    {
      if (nodeflag[conn_tet4[5 * i + j]] == 1)
        flag = true;
    }
    if (flag)
    {
      netet_outer++;
      for (j = 0; j < 5; ++j)
        conn_tet4_new.push_back(conn_tet4[5 * i + j]);
      for (j = 0; j < 11; ++j)
        conn_tet10_new.push_back(conn_alltet10[11 * i + j + disp]);
    }
  }

  for (i = 0; i < nevox; ++i)
  {
    flag = false;
    for (j = 0; j < 8; ++j)
    {
      if (nodeflag[conn_vox8[9 * i + j]] == 1)
        flag = true;
    }
    if (!flag)
    {
      for (j = 0; j < 9; ++j)
        conn_vox8_new.push_back(conn_vox8[9 * i + j]);
      for (j = 0; j < 66; ++j)
        conn_vox27_new.push_back(conn_alltet10[66 * i + j]);
    }
  }
  for (i = 0; i < netet; ++i)
  {
    flag = false;
    for (j = 0; j < 4; ++j)
    {
      if (nodeflag[conn_tet4[5 * i + j]] == 1)
        flag = true;
    }
    if (!flag)
    {
      for (j = 0; j < 5; ++j)
        conn_tet4_new.push_back(conn_tet4[5 * i + j]);
      for (j = 0; j < 11; ++j)
        conn_tet10_new.push_back(conn_alltet10[11 * i + j + disp]);
    }
  }

  if (myrank == 0)
  {
    std::cout << "sizes: " << conn_vox8.size() << " " << conn_tet4.size() << " " << conn_alltet10.size() << "\n";
    std::cout << "sizes: " << conn_vox8_new.size() << " " << conn_tet4_new.size();
    std::cout << " " << conn_vox27_new.size() << " " << conn_tet10_new.size() << "\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);

  conn_vox8 = conn_vox8_new;
  conn_tet4 = conn_tet4_new;
  conn_alltet10 = conn_vox27_new;
  for (int i = 0; i < conn_tet10_new.size(); ++i)
    conn_alltet10.push_back(conn_tet10_new[i]);
}

void renumber_nodes(std::vector<double> &coor_quad, std::vector<int> &nodetype, std::vector<int> &conn_alltet10,
                    std::vector<int> &nodemap_quad2linear, std::vector<int> &nodemap_linear2quad)
{
  //  std::cout << "renumbering quad mesh\n";
  int n_quad(nodemap_quad2linear.size());
  int n_linear(nodemap_linear2quad.size());
  int ne(conn_alltet10.size() / 11);
  std::vector<int> mapquad_new2old(n_quad, -1);
  std::vector<int> mapquad_old2new(n_quad, -1);
  int i, j, k, nt(0);
  for (i = 0; i < ne; ++i)
  {
    for (j = 0; j < 10; ++j)
    {
      k = conn_alltet10[11 * i + j];
      if (mapquad_old2new[k] == -1)
      {
        mapquad_old2new[k] = nt;
        mapquad_new2old[nt] = k;
        nt++;
      }
      conn_alltet10[11 * i + j] = mapquad_old2new[k];
    }
  }
  if (nt != n_quad)
    Error("renumber nodes");

  std::vector<double> coor_quad_(coor_quad);
  std::vector<int> nodetype_(nodetype);
  std::vector<int> nodemap_quad2linear_(nodemap_quad2linear);
  std::vector<int> nodemap_linear2quad_(nodemap_linear2quad);
  for (i = 0; i < n_quad; ++i)
  {
    coor_quad[3 * i] = coor_quad_[3 * mapquad_new2old[i]];
    coor_quad[3 * i + 1] = coor_quad_[3 * mapquad_new2old[i] + 1];
    coor_quad[3 * i + 2] = coor_quad_[3 * mapquad_new2old[i] + 2];
    nodetype[i] = nodetype_[mapquad_new2old[i]];
  }
  for (i = 0; i < n_quad; ++i)
    nodemap_quad2linear[i] = nodemap_quad2linear_[mapquad_new2old[i]];
  for (i = 0; i < n_linear; ++i)
    nodemap_linear2quad[i] = mapquad_old2new[nodemap_linear2quad_[i]];
}

void MPI_set_comm_quad_settings(const std::vector<double> &coor, const std::vector<double> &coor_quad,
                                double ds, const std::vector<int> &conn_alltet10, const std::vector<int> &neighborrankid,
                                const std::vector<int> &mpinode_linear_pointer, const std::vector<int> &mpinode_linear_index,
                                const std::vector<int> &nodemap_quad2linear, const std::vector<int> &nodemap_linear2quad,
                                std::vector<int> &mpinode_quad_pointer, std::vector<int> &mpinode_quad_index, int myrank, int totalrank, const std::vector<int> &nodetype)
{
#ifndef DEBUG_NO_MPI
  // -- check mpi settings - linear
  //  MPI_check_comm_settings(coor,ds,neighborrankid,
  //    mpinode_linear_pointer,mpinode_linear_index,myrank,totalrank,nodetype);

  int n_linear(coor.size() / 3);
  int n_quad(coor_quad.size() / 3);
  int ne(conn_alltet10.size() / 11);
  int nbrs(neighborrankid.size()), nbrrank, i, j, k;
  std::vector<int> nodeflag_linear;
  std::vector<int> nodeflag_quad;
  std::vector<int> elems;
  bool flag, tagrecv, tagsend;
  int nn(100000);
  double eps(0.01 * ds);
  double eps2(eps * eps);
  std::vector<std::vector<int>> mpinode_quad_test(nbrs);
  std::vector<MPI_Request> req_send(nbrs);
  std::vector<MPI_Request> req_recv(nbrs);
  std::vector<MPI_Status> stat_send(nbrs);
  std::vector<MPI_Status> stat_recv(nbrs);
  std::vector<int> nodesize_send(nbrs);
  std::vector<int> nodesize_recv(nbrs);

  std::vector<std::vector<double>> coorbuf_send(nbrs);
  std::vector<std::vector<double>> coorbuf_recv(nbrs);
  std::vector<std::vector<int>> typebuf_send(nbrs); // 1 for tet4 nodes, -1 for edge nodes
  std::vector<std::vector<int>> typebuf_recv(nbrs);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "-- 1\n"
              << std::flush;

  // extract quad nodes for each rank
  for (i = 0; i < nbrs; ++i)
  {
    //    nbrrank=neighborrankid[i];
    elems.clear();
    nodeflag_linear.clear();
    nodeflag_linear.resize(n_linear);
    for (j = 0; j < n_linear; ++j)
      nodeflag_linear[j] = 0;
    nodeflag_quad.clear();
    nodeflag_quad.resize(n_quad);
    for (j = 0; j < n_quad; ++j)
      nodeflag_quad[j] = -1;
    for (j = mpinode_linear_pointer[i]; j < mpinode_linear_pointer[i + 1]; ++j)
      nodeflag_linear[mpinode_linear_index[j]] = 1;

    for (j = 0; j < ne; ++j)
    {
      flag = false;
      for (k = 0; k < 4; ++k)
      {
        if (nodeflag_linear[nodemap_quad2linear[conn_alltet10[11 * j + k]]] == 1)
          flag = true;
      }
      if (flag)
      {
        elems.push_back(j);
        for (k = 0; k < 4; ++k)
          nodeflag_quad[conn_alltet10[11 * j + k]] = nodetype[conn_alltet10[11 * j + k]];
        for (k = 4; k < 10; ++k)
          nodeflag_quad[conn_alltet10[11 * j + k]] = nodetype[conn_alltet10[11 * j + k]];
      }
    }
    for (j = 0; j < n_quad; ++j)
    {
      if (nodeflag_quad[j] != -1)
      {
        mpinode_quad_test[i].push_back(j);
        coorbuf_send[i].push_back(coor_quad[3 * j]);
        coorbuf_send[i].push_back(coor_quad[3 * j + 1]);
        coorbuf_send[i].push_back(coor_quad[3 * j + 2]);
        typebuf_send[i].push_back(nodeflag_quad[j]);
      }
    }
    nodesize_send[i] = mpinode_quad_test[i].size();
    tagsend = 1;
    (myrank + 1) * nn + neighborrankid[i]; // tag=(sendid+1)*nn+recvid
    MPI_Isend(&(nodesize_send[i]), 1, MPI_INT, neighborrankid[i],
              tagsend, MPI_COMM_WORLD, &(req_send[i]));
    tagrecv = 1;
    (neighborrankid[i] + 1) * nn + myrank; // tag=(sendid+1)*nn+recvid
    MPI_Irecv(&(nodesize_recv[i]), 1, MPI_INT, neighborrankid[i],
              tagrecv, MPI_COMM_WORLD, &(req_recv[i]));
  }
  if (nbrs != 0)
  {
    MPI_Waitall(nbrs, &(req_recv[0]), &(stat_recv[0]));
    MPI_Waitall(nbrs, &(req_send[0]), &(stat_send[0]));
  }

  // allocate memory to recv coordinates
  int totnodes(0);
  for (i = 0; i < nbrs; ++i)
  {
    coorbuf_recv[i].resize(3 * nodesize_recv[i]);
    typebuf_recv[i].resize(nodesize_recv[i]);
    totnodes += nodesize_recv[i] + nodesize_send[i];
  }

  nn *= 2;
  // send recv coordinates
  for (i = 0; i < nbrs; ++i)
  {
    nodesize_send[i] = mpinode_quad_test[i].size();
    tagsend = 1;
    myrank *nn + neighborrankid[i]; // tag=sendid*nn+recvid
    tagsend = 1;
    (myrank + 1) * nn + neighborrankid[i]; // tag=(sendid+1)*nn+recvid
    MPI_Isend(&(coorbuf_send[i][0]), 3 * nodesize_send[i], MPI_DOUBLE, neighborrankid[i],
              tagsend, MPI_COMM_WORLD, &(req_send[i]));
    tagrecv = 1;
    (neighborrankid[i] + 1) * nn + myrank; // tag=(sendid+1)*nn+recvid
    MPI_Irecv(&(coorbuf_recv[i][0]), 3 * nodesize_recv[i], MPI_DOUBLE, neighborrankid[i],
              tagrecv, MPI_COMM_WORLD, &(req_recv[i]));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "-- 2\n"
              << std::flush;
  if (nbrs != 0)
  {
    MPI_Waitall(nbrs, &(req_recv[0]), &(stat_recv[0]));
    MPI_Waitall(nbrs, &(req_send[0]), &(stat_send[0]));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "-- 3\n"
              << std::flush;

  nn *= 2;
  // send recv types
  for (i = 0; i < nbrs; ++i)
  {
    nodesize_send[i] = mpinode_quad_test[i].size();
    tagsend = 1;
    myrank *nn + neighborrankid[i]; // tag=sendid*nn+recvid
    tagsend = 1;
    (myrank + 1) * nn + neighborrankid[i]; // tag=(sendid+1)*nn+recvid
    MPI_Isend(&(typebuf_send[i][0]), nodesize_send[i], MPI_INT, neighborrankid[i],
              tagsend, MPI_COMM_WORLD, &(req_send[i]));
    tagrecv = 1;
    (neighborrankid[i] + 1) * nn + myrank; // tag=(sendid+1)*nn+recvid
    MPI_Irecv(&(typebuf_recv[i][0]), nodesize_recv[i], MPI_INT, neighborrankid[i],
              tagrecv, MPI_COMM_WORLD, &(req_recv[i]));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "-- 4\n"
              << std::flush;
  if (nbrs != 0)
  {
    MPI_Waitall(nbrs, &(req_recv[0]), &(stat_recv[0]));
    MPI_Waitall(nbrs, &(req_send[0]), &(stat_send[0]));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "-- 5\n"
              << std::flush;
  if (myrank == 0)
    std::cout << "ds: " << ds << "\n"
              << std::flush;

  // check match coorbuf_send and coorbuf_recv
  int nsend, nrecv;
  std::vector<int> tmpmatch;
  std::vector<int> tmpmatch_index(totnodes);
  mpinode_quad_pointer.resize(nbrs + 1);
  totnodes = 0;
  int ncount, typeij;
  double tx, ty, tz;
  std::vector<int>::iterator is, ie, it;

  for (i = 0; i < nbrs; ++i)
  {

    nsend = nodesize_send[i];
    nrecv = nodesize_recv[i];
    tmpmatch.clear();

    std::map<long, std::vector<int>> map;
    for (k = 0; k < nrecv; ++k)
      (map[coorbuf_recv[i][3 * k] / ds]).push_back(k);

    for (j = 0; j < nsend; ++j)
    {
      ncount = 0;
      /*
            typeij=typebuf_send[i][j];
            tx=coorbuf_send[i][3*j+0];
            ty=coorbuf_send[i][3*j+1];
            tz=coorbuf_send[i][3*j+2];

            for( k=0; k<nrecv; ++k)
            {
              if( dist2(&(coorbuf_send[i][3*j]),&(coorbuf_recv[i][3*k]))<eps2 &&
                typebuf_send[i][j]==typebuf_recv[i][k] )
              {
                tmpmatch.push_back(mpinode_quad_test[i][j]);
                ncount++;
              }
            }
      */
      //      if( (map[coorbuf_send[i][3*j]/ds]).size() == 0 ) Error("something wrong in map");
      if ((map[coorbuf_send[i][3 * j] / ds]).size() != 0)
      {

        is = (map[coorbuf_send[i][3 * j] / ds]).begin();
        ie = (map[coorbuf_send[i][3 * j] / ds]).end();
        for (it = is; it != ie; ++it)
        {
          if (dist2(&(coorbuf_send[i][3 * j]), &(coorbuf_recv[i][3 * (*it)])) < eps2 &&
              typebuf_send[i][j] == typebuf_recv[i][(*it)])
          {
            tmpmatch.push_back(mpinode_quad_test[i][j]);
            ncount++;
          }
        }
      }

      if (ncount > 1)
      {
        std::cout << "ncount,isend,irecv: " << ncount << " " << myrank << " " << neighborrankid[i] << "\n";
        ncount = mpinode_quad_test[i][j];
        std::cout << ncount << " " << coor_quad[3 * ncount] << " " << coor_quad[3 * ncount + 1] << " " << coor_quad[3 * ncount + 2] << "\n";
      }
    }
    global_coor_pair_eps = eps;
    sort_by_coordinate(coor_quad, nodetype, tmpmatch);
    mpinode_quad_pointer[i] = totnodes;
    for (j = 0; j < tmpmatch.size(); ++j)
      tmpmatch_index[j + totnodes] = tmpmatch[j];
    totnodes += tmpmatch.size();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "-- 6\n"
              << std::flush;
  mpinode_quad_pointer[nbrs] = totnodes;
  mpinode_quad_index.resize(totnodes);
  for (j = 0; j < totnodes; ++j)
    mpinode_quad_index[j] = tmpmatch_index[j];

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "-- 7\n"
              << std::flush;
  // -- check mpi settings - quad
  MPI_check_comm_settings(coor_quad, ds, neighborrankid,
                          mpinode_quad_pointer, mpinode_quad_index, myrank, totalrank, nodetype);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "-- 8\n"
              << std::flush;
#endif
}

void MPI_set_comm_linear_settings(const std::vector<double> &coor, double ds, const std::vector<int> &neighborrankid,
                                  std::vector<int> &mpinode_linear_pointer, std::vector<int> &mpinode_linear_index,
                                  int myrank, int totalrank, const std::vector<int> &nodetype)
{
  double eps(ds * 0.01);
  std::vector<int> tmparr;
  for (int i = 0; i < neighborrankid.size(); ++i)
  {
    tmparr.resize(mpinode_linear_pointer[i + 1] - mpinode_linear_pointer[i]);
    std::copy(&(mpinode_linear_index[mpinode_linear_pointer[i]]),
              &(mpinode_linear_index[mpinode_linear_pointer[i + 1]]), tmparr.begin());
    global_coor_pair_eps = eps;
    sort_by_coordinate(coor, nodetype, tmparr);
    std::copy(tmparr.begin(), tmparr.end(), &(mpinode_linear_index[mpinode_linear_pointer[i]]));
  }
}

void MPI_extract_surface_nodes_quad(const std::vector<int> &surfnodes, const std::vector<int> &nodemap_linear2quad,
                                    double ds, const std::vector<double> &coor_quad, const std::vector<int> &conn_tet10, std::vector<int> &surfnodes_quad, int myrank,
                                    const std::vector<int> &neighborrankid, const std::vector<int> &mpinode_quad_pointer, const std::vector<int> &mpinode_quad_index)
{
  int i, j;
  int n_quad(coor_quad.size() / 3);
  int n_linear(nodemap_linear2quad.size());
  int ne(conn_tet10.size() / 11);
  int cnyele[6][2];
  cnyele[0][0] = 0;
  cnyele[0][1] = 1;
  cnyele[1][0] = 1;
  cnyele[1][1] = 2;
  cnyele[2][0] = 0;
  cnyele[2][1] = 2;
  cnyele[3][0] = 0;
  cnyele[3][1] = 3;
  cnyele[4][0] = 1;
  cnyele[4][1] = 3;
  cnyele[5][0] = 3;
  cnyele[5][1] = 2;

  std::vector<int> surfflag(n_quad, -1);
  for (i = 0; i < surfnodes.size(); ++i)
    surfflag[nodemap_linear2quad[surfnodes[i]]] = 1;
#define INCLUDE_SURF_QUAD
#ifdef INCLUDE_SURF_QUAD
  for (i = 0; i < ne; ++i)
  {
    for (j = 0; j < 6; ++j)
    {
      if (surfflag[conn_tet10[11 * i + cnyele[j][0]]] == 1 &&
          surfflag[conn_tet10[11 * i + cnyele[j][1]]] == 1)
        surfflag[conn_tet10[11 * i + j + 4]] = 1;
    }
  }
#endif

  // exclude overlapping mpi nodes
  int nbrs(neighborrankid.size());
  for (i = 0; i < nbrs; ++i)
  {
    if (neighborrankid[i] < myrank)
      continue;
    for (j = mpinode_quad_pointer[i]; j < mpinode_quad_pointer[i + 1]; ++j)
    {
      surfflag[mpinode_quad_index[j]] = -1;
    }
  }

  surfnodes_quad.clear();
  for (j = 0; j < n_quad; ++j)
  {
    if (surfflag[j] != -1)
      surfnodes_quad.push_back(j);
  }

  std::ofstream fout;

  long localsize(surfnodes_quad.size());
  check_open("./mdata/" + ToString(myrank, 6, '0') + ".2Doutputflag.dat", fout);
  fout << localsize << "\n";
  for (j = 0; j < localsize; ++j)
    fout << surfnodes_quad[j] + 1 << "\n"; // ToFortranStyle
  for (j = 0; j < localsize; ++j)
  {
    i = surfnodes_quad[j];
    fout << coor_quad[3 * i] << " " << coor_quad[3 * i + 1] << " " << coor_quad[3 * i + 2] << "\n";
  }
  fout.close();

  long glbsize;
#ifndef DEBUG_NO_MPI
  MPI_Reduce(&localsize, &glbsize, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  glbsize = localsize;
#endif
  if (myrank == 0)
  {
    check_open("./data/2Doutput.dat", fout);
    fout << glbsize << "\n";
    fout.close();
  }
}

void MPI_extract_3D_nodes_quad(const std::vector<double> &coor_quad, const std::vector<int> &conn_tet10,
                               std::vector<int> &surfnodes_quad, int myrank, double minelev)
{
  minelev = 300000;
  if (myrank == 0)
    std::cout << "output 3D node with surface and cut at y=" << minelev << "\n";
  //    std::cout << "output 3D node. minelev " << minelev << "\n";

  int i, j, k;
  int n_quad(coor_quad.size() / 3);

  std::vector<int> surfflag(n_quad, -1);
  for (i = 0; i < n_quad; ++i)
  {
    if (fabs(coor_quad[3 * i + 1] - minelev) < 1.0)
      surfflag[i] = 1;
  }

  // extract surface
  double minx, miny, maxx, maxy, ds, minx2;
  minx = maxx = coor_quad[0];
  miny = maxy = coor_quad[1];
  for (i = 0; i < n_quad; ++i)
  {
    minx = std::min(coor_quad[3 * i], minx);
    maxx = std::max(coor_quad[3 * i], maxx);
    miny = std::min(coor_quad[3 * i + 1], miny);
    maxy = std::max(coor_quad[3 * i + 1], maxy);
  }
  minx2 = maxx;
  for (i = 0; i < n_quad; ++i)
  {
    if (fabs(coor_quad[3 * i] - minx) < 0.001)
      continue;
    minx2 = std::min(coor_quad[3 * i], minx2);
  }
  ds = minx2 - minx;
  double dstot;
  MPI_Allreduce(&ds, &dstot, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  ds = dstot;
  if (myrank == 0)
    std::cout << "minx,maxx,miny,maxy,ds: " << minx << " " << maxx << " " << miny << " " << maxy << " " << ds << "\n";
  int nx((maxx - minx) / ds + 1);
  int ny((maxy - miny) / ds + 1);
  std::vector<double> surfelev(nx * ny, -10000.0);
  int ix, iy;
  for (i = 0; i < n_quad; ++i)
  {
    ix = (coor_quad[3 * i] - minx) / ds;
    iy = (coor_quad[3 * i + 1] - miny) / ds;
    if (ix < 0)
      Error("wrong here 11");
    if (iy < 0)
      Error("wrong here 12");
    if (ix > nx - 1)
      Error("wrong here 13");
    if (iy > ny - 1)
      Error("wrong here 14");
    surfelev[ix + iy * nx] = std::max(surfelev[ix + iy * nx], coor_quad[3 * i + 2]);
  }
  for (i = 0; i < n_quad; ++i)
  {
    ix = (coor_quad[3 * i] - minx) / ds;
    iy = (coor_quad[3 * i + 1] - miny) / ds;
    if (coor_quad[3 * i + 2] > surfelev[ix + iy * nx] - 3.0) // 3.0 m depth
      surfflag[i] = 1;
  }

  bool flag;
  std::vector<int> elementlist;
  for (i = 0; i < conn_tet10.size() / 11; ++i)
  {
    flag = false;
    for (j = 0; j < 10; ++j)
    {
      if (surfflag[conn_tet10[11 * i + j]] == 1)
        flag = true;
    }
    if (flag)
    {
      elementlist.push_back(i);
      for (j = 0; j < 10; ++j)
      {
        if (surfflag[conn_tet10[11 * i + j]] != 1)
          surfflag[conn_tet10[11 * i + j]] = 2;
      }
    }
  }

  surfnodes_quad.clear();
  std::vector<int> map_surf(n_quad, -1);
  for (j = 0; j < n_quad; ++j)
  {
    if (surfflag[j] != -1)
    {
      map_surf[j] = surfnodes_quad.size();
      surfnodes_quad.push_back(j);
    }
  }

  long localsize(surfnodes_quad.size());
  std::vector<float> coor(localsize * 3);
  std::vector<int> conn(11 * elementlist.size());
  std::vector<int> inttype(elementlist.size(), 24);
  std::vector<int> matid(elementlist.size());

  std::ofstream fout;
  check_open("./mdata/" + ToString(myrank, 6, '0') + ".2Doutputflag.dat", fout);
  fout << localsize << "\n";
  for (j = 0; j < localsize; ++j)
    fout << surfnodes_quad[j] + 1 << "\n"; // ToFortranStyle
  for (j = 0; j < localsize; ++j)
  {
    i = surfnodes_quad[j];
    coor[3 * j + 0] = coor_quad[3 * i + 0];
    coor[3 * j + 1] = coor_quad[3 * i + 1];
    coor[3 * j + 2] = coor_quad[3 * i + 2];
    fout << coor_quad[3 * i] << " " << coor_quad[3 * i + 1] << " " << coor_quad[3 * i + 2] << "\n";
  }
  fout << elementlist.size() << "\n"; // ToFortranStyle
  for (j = 0; j < elementlist.size(); ++j)
  {
    i = elementlist[j];
    conn[11 * j] = 10;
    for (k = 0; k < 10; ++k)
    {
      conn[11 * j + k + 1] = map_surf[conn_tet10[11 * i + k]];
      fout << map_surf[conn_tet10[11 * i + k]] + 1 << " "; // ToFortranStyle
    }
    fout << "\n";
  }
  fout << "material\n"; // ToFortranStyle
  for (j = 0; j < elementlist.size(); ++j)
  {
    i = elementlist[j];
    matid[j] = conn_tet10[11 * i + 10];
    fout << conn_tet10[11 * i + 10] + 1 << "\n"; // ToFortranStyle
  }
  fout.close();

#ifdef IAMINTEL
  vec_swap_bytes(coor);
  vec_swap_bytes(conn);
  vec_swap_bytes(inttype);
  vec_swap_bytes(matid);
#endif

  check_open("./mdata/output_part." + ToString(myrank, 6, '0') + ".vtk", fout);
  fout << "# vtk DataFile Version 2.0\n";
  fout << "file:output_part." + ToString(myrank, 6, '0') + ".vtk\n";
  fout << "BINARY\n";
  fout << "DATASET UNSTRUCTURED_GRID\n";
  fout << "POINTS " << coor.size() / 3 << " float\n";
  fout.write((char *)&(coor[0]), sizeof(float) * coor.size());
  fout << "CELLS " << conn.size() / 11 << " " << conn.size() << std::endl;
  fout.write((char *)&(conn[0]), sizeof(int) * conn.size());
  fout << "CELL_TYPES " << inttype.size() << std::endl;
  fout.write((char *)&(inttype[0]), sizeof(int) * inttype.size());
  fout << "CELL_DATA " << inttype.size() << std::endl;
  fout << "SCALARS material int 1\n";
  fout << "LOOKUP_TABLE default\n";
  fout.write((char *)&(matid[0]), sizeof(int) * inttype.size());
  fout.close();
#ifdef OUTPUT_VTK
#endif

  long glbsize;
#ifndef DEBUG_NO_MPI
  MPI_Reduce(&localsize, &glbsize, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  glbsize = localsize;
#endif
  if (myrank == 0)
  {
    check_open("./data/2Doutput.dat", fout);
    fout << glbsize << "\n";
    fout.close();
  }
}

#ifndef DEBUG_NO_HDF5
void write_hdf_numofelements(long fid, int ne, int nthread)
{
  std::vector<int> tmparr;
  for (int i = 0; i < nthread; ++i)
  {
    tmparr.push_back(i + 1);
    tmparr.push_back(ne / nthread);
  }
  tmparr[2 * nthread - 1] -= ((ne / nthread) * nthread - ne);
  long size(tmparr.size());
  hdf_write_int_array_(&fid, "/numofelements", (int *)&(tmparr[0]), &size);
}

void write_hdf_coor(long fid, const std::vector<double> &coor)
{
  long size(coor.size());
  hdf_write_double_array_(&fid, "/coor", (double *)&(coor[0]), &size);
}

void write_hdf_conn(long fid, const std::vector<int> &conn)
{
  std::vector<int> tmparr;
  long size(conn.size());
  ToFortranStyle(conn, tmparr);
  hdf_write_int_array_(&fid, "/conn", (int *)&(tmparr[0]), &size);
}

void write_hdf_MPInode(long fid, const std::vector<int> &neighborrankid,
                       const std::vector<int> &mpinode_pointer, const std::vector<int> &mpinode_index)
{
  std::vector<int> tmparr;
  tmparr.push_back(neighborrankid.size());
  for (int i = 0; i < neighborrankid.size(); ++i)
  {
    tmparr.push_back(neighborrankid[i]);
    tmparr.push_back(mpinode_pointer[i + 1] - mpinode_pointer[i]);
  }
  for (int i = 0; i < mpinode_index.size(); ++i)
    tmparr.push_back(mpinode_index[i] + 1); // ToFortranStyle
  long size(tmparr.size());
  hdf_write_int_array_(&fid, "/MPInode", (int *)&(tmparr[0]), &size);
}

void write_hdf_BC(long fid, const std::vector<int> &BC)
{
  std::vector<int> tmparr;
  tmparr.push_back(BC.size());
  for (int i = 0; i < BC.size(); ++i)
    tmparr.push_back(BC[i] + 1); // ToFortranStyle
  long size(tmparr.size());
  hdf_write_int_array_(&fid, "/BC", (int *)&(tmparr[0]), &size);
}

void write_hdf_duplicate(long fid, int n, const std::vector<int> &mpinode_index)
{
  std::vector<int> duplicate(n, 1);
  long size(n);
  for (int i = 0; i < mpinode_index.size(); ++i)
    duplicate[mpinode_index[i]]++;
  hdf_write_int_array_(&fid, "/duplicate", (int *)&(duplicate[0]), &size);
}

void write_hdf_global_node_id(long fid, const std::vector<int> &global_node_id)
{
  if (global_node_id.empty())
    return;
  long size(global_node_id.size());
  hdf_write_int_array_(&fid, "/global_node_id", (int *)&(global_node_id[0]), &size);
}

void write_hdf_global_node_id64(long fid, const std::vector<long> &global_node_id64)
{
  if (global_node_id64.empty())
    return;
  long size(global_node_id64.size());
  hdf_write_long_array_(&fid, "/global_node_id64", (long *)&(global_node_id64[0]), &size);
}

void write_hdf_mdata(const std::vector<double> &coor, const std::vector<int> &conn_alltet10,
                     const std::vector<int> &neighborrankid,
                     const std::vector<int> &mpinode_pointer, const std::vector<int> &mpinode_index,
                     const std::vector<int> &BC, int nmat, long fid,
                     const std::vector<int> &global_node_id,
                     const std::vector<long> &global_node_id64)
{
  int setting[3];
  setting[0] = coor.size() / 3;           // number of nodes
  setting[1] = conn_alltet10.size() / 11; // number of tet10 elements
  setting[2] = nmat;
  long size(3);
  hdf_write_int_array_(&fid, "/setting", setting, &size);

  write_hdf_coor(fid, coor);
  write_hdf_conn(fid, conn_alltet10);
  write_hdf_MPInode(fid, neighborrankid, mpinode_pointer, mpinode_index);
  write_hdf_BC(fid, BC);
  write_hdf_duplicate(fid, setting[0], mpinode_index);
  write_hdf_global_node_id(fid, global_node_id);
  write_hdf_global_node_id64(fid, global_node_id64);
}

void write_hdf_cdata(const std::vector<double> &coor, const std::vector<int> &conn_vox8,
                     const std::vector<int> &conn_tet4, const std::vector<int> &neighborrankid,
                     const std::vector<int> &mpinode_pointer, const std::vector<int> &mpinode_index,
                     const std::vector<int> &BC, int nmat, long fid,
                     const std::vector<int> &global_node_id,
                     const std::vector<long> &global_node_id64)
{
  std::vector<int> conn_alltet4;
  convert_vox8tet4_to_tet4(conn_vox8, conn_tet4, conn_alltet4);

  // setting
  int setting[3];
  setting[0] = coor.size() / 3;         // number of nodes
  setting[1] = conn_alltet4.size() / 5; // number of elements
  setting[2] = nmat;                    // nmat
  long size(3);
  hdf_write_int_array_(&fid, "/setting", setting, &size);

  write_hdf_coor(fid, coor);
  write_hdf_conn(fid, conn_alltet4);
  write_hdf_MPInode(fid, neighborrankid, mpinode_pointer, mpinode_index);
  write_hdf_BC(fid, BC);
  write_hdf_duplicate(fid, setting[0], mpinode_index);
  write_hdf_global_node_id(fid, global_node_id);
  write_hdf_global_node_id64(fid, global_node_id64);
}
#endif

void write_settings(const std::vector<double> &coor_quad, const std::vector<int> &conn_alltet10,
                    const std::vector<int> &neighborrankid, const std::vector<int> &mpinode_quad_pointer,
                    const std::vector<int> &mpinode_quad_index, int nmat, int myrank)
{
  long ne(conn_alltet10.size() / 11);
  long nall(coor_quad.size() / 3), n(0);
  std::vector<int> flag(nall, 1);
  int nbrs(neighborrankid.size());
  for (int i = 0; i < nbrs; ++i)
  {
    for (int j = mpinode_quad_pointer[i]; j < mpinode_quad_pointer[i + 1]; ++j)
    {
      if (myrank > neighborrankid[i])
        flag[mpinode_quad_index[j]] = 0;
    }
  }
  for (int i = 0; i < nall; ++i)
    n += flag[i];

  long n_g, ne_g;
#ifndef DEBUG_NO_MPI
  MPI_Reduce(&n, &n_g, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ne, &ne_g, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  n_g = n;
  ne_g = ne;
#endif
  if (myrank == 0)
  {
    std::ofstream fout;
    check_open("./data/setting.dat", fout);
    fout << "num of node\n";
    fout << n_g << "\n";
    fout << "num of tet elem\n";
    fout << ne_g << "\n";
    fout << "num of material properties\n";
    fout << nmat << "\n";
    fout.close();
  }
}

void set_rmat(int myrank, int nmat, std::vector<double> &rmat)
{
  rmat.resize(nmat * 10);
  if (myrank == 0)
  {
    std::string temp;
    std::ifstream fin;
    check_open("./data/material.dat", fin);
    for (int i = 0; i < nmat; ++i)
    {
      fin >> temp; // material
      fin >> temp; // "material-number"
      for (int j = 0; j < 10; ++j)
        fin >> rmat[10 * i + j];
    }
    fin.close();
  }
#ifndef DEBUG_NO_MPI
  MPI_Bcast(&(rmat[0]), rmat.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void convert_quad2linear(const std::vector<double> &coor_quad, const std::vector<int> &conn_alltet10,
                         int netet, int nevox, std::vector<double> &coor, std::vector<int> &conn_vox8, std::vector<int> &conn_tet4)
{
  int i, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
  int n11, n12, n13, n14, n15, n16, n17, n18, n19, n20;
  int n21, n22, n23, n24, n25, n26, matid, ind, offset;
  conn_vox8.clear();
  conn_tet4.clear();
  for (i = 0; i < nevox; ++i)
  {
    matid = conn_alltet10[11 * 6 * i + 10];
    ind = 11 * 6 * i;
    n0 = conn_alltet10[ind + 11 * 1 + 1];
    n1 = conn_alltet10[ind + 11 * 2 + 1];
    n2 = conn_alltet10[ind + 11 * 0 + 0];
    n3 = conn_alltet10[ind + 11 * 1 + 3];
    n4 = conn_alltet10[ind + 11 * 0 + 2];
    n5 = conn_alltet10[ind + 11 * 0 + 3];
    n6 = conn_alltet10[ind + 11 * 0 + 1];
    n7 = conn_alltet10[ind + 11 * 5 + 3];
    n8 = conn_alltet10[ind + 11 * 2 + 8];
    n9 = conn_alltet10[ind + 11 * 2 + 4];
    n10 = conn_alltet10[ind + 11 * 1 + 7];
    n11 = conn_alltet10[ind + 11 * 1 + 8];
    n12 = conn_alltet10[ind + 11 * 1 + 5];
    n13 = conn_alltet10[ind + 11 * 3 + 8];
    n14 = conn_alltet10[ind + 11 * 0 + 4];
    n15 = conn_alltet10[ind + 11 * 5 + 7];
    n16 = conn_alltet10[ind + 11 * 0 + 9];
    n17 = conn_alltet10[ind + 11 * 0 + 8];
    n18 = conn_alltet10[ind + 11 * 5 + 9];
    n19 = conn_alltet10[ind + 11 * 5 + 8];
    n20 = conn_alltet10[ind + 11 * 1 + 4];
    n21 = conn_alltet10[ind + 11 * 3 + 9];
    n22 = conn_alltet10[ind + 11 * 0 + 7];
    n23 = conn_alltet10[ind + 11 * 4 + 8];
    n24 = conn_alltet10[ind + 11 * 5 + 4];
    n25 = conn_alltet10[ind + 11 * 5 + 5];
    n26 = conn_alltet10[ind + 11 * 4 + 6];
    // HERE
    n0 = conn_alltet10[ind + 11 * 3 + 3]; //
    n1 = conn_alltet10[ind + 11 * 0 + 0]; //
    n2 = conn_alltet10[ind + 11 * 1 + 0]; //
    n3 = conn_alltet10[ind + 11 * 0 + 2]; //
    n4 = conn_alltet10[ind + 11 * 0 + 3]; //
    n5 = conn_alltet10[ind + 11 * 1 + 1]; //
    n6 = conn_alltet10[ind + 11 * 1 + 2]; //
    n7 = conn_alltet10[ind + 11 * 1 + 3]; //

    n8 = conn_alltet10[ind + 11 * 3 + 8];  //
    n9 = conn_alltet10[ind + 11 * 2 + 6];  //
    n10 = conn_alltet10[ind + 11 * 2 + 9]; //
    n11 = conn_alltet10[ind + 11 * 3 + 7]; //
    n12 = conn_alltet10[ind + 11 * 3 + 9]; //
    n13 = conn_alltet10[ind + 11 * 4 + 8]; //
    n14 = conn_alltet10[ind + 11 * 1 + 6]; //
    n15 = conn_alltet10[ind + 11 * 0 + 5]; //
    n16 = conn_alltet10[ind + 11 * 5 + 7]; //
    n17 = conn_alltet10[ind + 11 * 1 + 5]; //
    n18 = conn_alltet10[ind + 11 * 1 + 9]; //
    n19 = conn_alltet10[ind + 11 * 0 + 8]; //
    n20 = conn_alltet10[ind + 11 * 0 + 6]; //
    n21 = conn_alltet10[ind + 11 * 0 + 7]; //
    n22 = conn_alltet10[ind + 11 * 1 + 4]; //
    n23 = conn_alltet10[ind + 11 * 1 + 7]; //
    n24 = conn_alltet10[ind + 11 * 0 + 9]; //
    n25 = conn_alltet10[ind + 11 * 1 + 8]; //
    n26 = conn_alltet10[ind + 11 * 0 + 4]; //

    // element 0
    conn_vox8.push_back(n0);
    conn_vox8.push_back(n8);
    conn_vox8.push_back(n20);
    conn_vox8.push_back(n11);
    conn_vox8.push_back(n12);
    conn_vox8.push_back(n21);
    conn_vox8.push_back(n26);
    conn_vox8.push_back(n24);
    conn_vox8.push_back(matid);
    // element 1
    conn_vox8.push_back(n8);
    conn_vox8.push_back(n1);
    conn_vox8.push_back(n9);
    conn_vox8.push_back(n20);
    conn_vox8.push_back(n21);
    conn_vox8.push_back(n13);
    conn_vox8.push_back(n22);
    conn_vox8.push_back(n26);
    conn_vox8.push_back(matid);
    // element 2
    conn_vox8.push_back(n20);
    conn_vox8.push_back(n9);
    conn_vox8.push_back(n2);
    conn_vox8.push_back(n10);
    conn_vox8.push_back(n26);
    conn_vox8.push_back(n22);
    conn_vox8.push_back(n14);
    conn_vox8.push_back(n23);
    conn_vox8.push_back(matid);
    // element 3
    conn_vox8.push_back(n11);
    conn_vox8.push_back(n20);
    conn_vox8.push_back(n10);
    conn_vox8.push_back(n3);
    conn_vox8.push_back(n24);
    conn_vox8.push_back(n26);
    conn_vox8.push_back(n23);
    conn_vox8.push_back(n15);
    conn_vox8.push_back(matid);
    // element 4
    conn_vox8.push_back(n12);
    conn_vox8.push_back(n21);
    conn_vox8.push_back(n26);
    conn_vox8.push_back(n24);
    conn_vox8.push_back(n4);
    conn_vox8.push_back(n16);
    conn_vox8.push_back(n25);
    conn_vox8.push_back(n19);
    conn_vox8.push_back(matid);
    // element 5
    conn_vox8.push_back(n21);
    conn_vox8.push_back(n13);
    conn_vox8.push_back(n22);
    conn_vox8.push_back(n26);
    conn_vox8.push_back(n16);
    conn_vox8.push_back(n5);
    conn_vox8.push_back(n17);
    conn_vox8.push_back(n25);
    conn_vox8.push_back(matid);
    // element 6
    conn_vox8.push_back(n26);
    conn_vox8.push_back(n22);
    conn_vox8.push_back(n14);
    conn_vox8.push_back(n23);
    conn_vox8.push_back(n25);
    conn_vox8.push_back(n17);
    conn_vox8.push_back(n6);
    conn_vox8.push_back(n18);
    conn_vox8.push_back(matid);
    // element 7
    conn_vox8.push_back(n24);
    conn_vox8.push_back(n26);
    conn_vox8.push_back(n23);
    conn_vox8.push_back(n15);
    conn_vox8.push_back(n19);
    conn_vox8.push_back(n25);
    conn_vox8.push_back(n18);
    conn_vox8.push_back(n7);
    conn_vox8.push_back(matid);
  }
  offset = 11 * 6 * nevox;
  for (i = 0; i < netet; ++i)
  {
    matid = conn_alltet10[offset + 11 * i + 10];
    // element 0
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 0]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 4]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 6]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 7]);
    conn_tet4.push_back(matid);
    // element 1
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 4]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 1]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 5]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 8]);
    conn_tet4.push_back(matid);
    // element 2
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 6]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 5]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 2]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 9]);
    conn_tet4.push_back(matid);
    // element 3
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 7]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 8]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 9]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 3]);
    conn_tet4.push_back(matid);
    // element 4
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 6]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 8]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 9]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 7]);
    conn_tet4.push_back(matid);
    // element 5
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 6]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 4]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 8]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 7]);
    conn_tet4.push_back(matid);
    // element 6
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 4]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 8]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 5]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 6]);
    conn_tet4.push_back(matid);
    // element 7
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 5]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 8]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 9]);
    conn_tet4.push_back(conn_alltet10[offset + 11 * i + 6]);
    conn_tet4.push_back(matid);
  }
  coor = coor_quad;
}

void refine_linear_mesh(
    double dsorg, double gmaxx, double gmaxy,
    std::vector<double> &coor, std::vector<int> &conn_vox8, std::vector<int> &conn_tet4,
    double &ds, int myrank, int totalrank, int neighbors, const std::vector<int> &neighborrankid,
    std::vector<int> &mpinode_linear_pointer, std::vector<int> &mpinode_linear_index, std::vector<int> &surfnodes,
    std::vector<int> &nodetype)
{
  std::vector<double> coor_quad;
  std::vector<int> conn_alltet10; // voxel elements are stored in beginning part, tetra elements in latter part
  std::vector<int> nodemap_quad2linear;
  std::vector<int> nodemap_linear2quad;
  std::vector<int> surfnodes_quad;
  std::vector<int> mpinode_quad_pointer;
  std::vector<int> mpinode_quad_index;
  int netet(conn_tet4.size() / 5);
  int nevox(conn_vox8.size() / 9);
  bool merge4and10(false);

  convert_linear2quad(
      myrank, dsorg, gmaxx, gmaxy,
      merge4and10, coor, conn_vox8, conn_tet4, ds,
      coor_quad, conn_alltet10,
      nodemap_quad2linear, nodemap_linear2quad, nodetype);

  std::vector<int> nodetype_linear(nodetype);
#ifdef RENUMBER_NODES
  renumber_nodes(coor_quad, nodetype, conn_alltet10, nodemap_quad2linear, nodemap_linear2quad);
#endif

  MPI_set_comm_linear_settings(coor, ds, neighborrankid,
                               mpinode_linear_pointer, mpinode_linear_index,
                               myrank, totalrank, nodetype_linear);
  MPI_set_comm_quad_settings(coor, coor_quad, ds, conn_alltet10, neighborrankid,
                             mpinode_linear_pointer, mpinode_linear_index,
                             nodemap_quad2linear, nodemap_linear2quad,
                             mpinode_quad_pointer, mpinode_quad_index, myrank, totalrank, nodetype);
  MPI_extract_surface_nodes_quad(surfnodes, nodemap_linear2quad, ds,
                                 coor_quad, conn_alltet10, surfnodes_quad,
                                 myrank, neighborrankid, mpinode_quad_pointer, mpinode_quad_index);
  mpinode_linear_pointer = mpinode_quad_pointer;
  mpinode_linear_index = mpinode_quad_index;
  surfnodes = surfnodes_quad;
  convert_quad2linear(coor_quad, conn_alltet10, netet, nevox, coor, conn_vox8, conn_tet4);
  ds *= 0.5;
}

//******************************************************************************
// main function
//
#ifdef MERGE_SURFACE_VTK
int main(int argc, char **argv)
{
  int myrank, totalrank;
#ifndef DEBUG_NO_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &totalrank);
#else
  myrank = 0;
  totalrank = 1;
#endif
  if (totalrank != 1)
    Error("use only one mpi process for MERGE_SURFACE_VTK");
  int nproc, nthread;
  std::string temp;
  double ds, gminx, gminy, gminz, gmaxx, gmaxy;
  int nmat;
  read_global_settings(myrank, ds, gminx, gminy, gminz, gmaxx, gmaxy, nmat, nproc, nthread);
  std::cout << "nproc,ds: " << nproc << " " << ds << "\n";
  double dsh(ds * 0.5 * NROUGH);
  int nx(round(gmaxx - gminx) / dsh + 1);
  int ny(round(gmaxy - gminy) / dsh + 1);
  std::vector<double> zval(nx * ny, gminz);
  std::vector<int> zindex(nx * ny, -1);
  int localnodes, totnodes, itmp, j, k, ix, iy;
  totnodes = 0;
  double x, y, z;
  std::vector<float> coor;
  std::vector<int> rankid;
  for (int iproc = 0; iproc < nproc; ++iproc)
  {
    std::ifstream fin1;
    std::string temp1;
    std::cout << "reading rank " << iproc << "\n";
    check_open("./mdata/" + ToString(iproc, 6, '0') + ".2Doutputflag.dat", fin1);
    fin1 >> temp1;
    localnodes = atoi(temp1.c_str());
    for (j = 0; j < localnodes; ++j)
      fin1 >> temp1; // local node id
    for (j = 0; j < localnodes; ++j)
    {
      fin1 >> temp1;
      x = atof(temp1.c_str());
      fin1 >> temp1;
      y = atof(temp1.c_str());
      fin1 >> temp1;
      z = atof(temp1.c_str());
      ix = round((x - gminx) / dsh);
      iy = round((y - gminy) / dsh);
      if (ix < 0)
        ix = 0;
      if (iy < 0)
        iy = 0;
      if (ix > nx - 1)
        ix = nx - 1;
      if (iy > ny - 1)
        iy = ny - 1;
      if (fabs(ix * dsh - x + gminx) < 0.01 * dsh &&
          fabs(iy * dsh - y + gminy) < 0.01 * dsh && zval[nx * iy + ix] < z)
      {
        zval[nx * iy + ix] = z;
        zindex[nx * iy + ix] = totnodes;
      }
      totnodes++;
      coor.push_back(x);
      coor.push_back(y);
      coor.push_back(z);
      rankid.push_back(iproc);
    }
    fin1.close();
  }

  int n(totnodes), i;
  int ne((nx - 1) * (ny - 1));
  std::vector<int> conn(5 * ne);
  std::vector<int> inttype(ne, 9); // 9 for square elements

  for (iy = 0; iy < ny - 1; ++iy)
    for (ix = 0; ix < nx - 1; ++ix)
    {
      i = (nx - 1) * iy + ix;
      conn[5 * i] = 4;
      conn[5 * i + 1] = zindex[nx * iy + ix];
      conn[5 * i + 2] = zindex[nx * iy + ix + 1];
      conn[5 * i + 3] = zindex[nx * (iy + 1) + ix + 1];
      conn[5 * i + 4] = zindex[nx * (iy + 1) + ix];
    }
  vec_swap_bytes(coor);
  vec_swap_bytes(conn);
  vec_swap_bytes(inttype);
  vec_swap_bytes(rankid);

  std::ofstream fout;

  std::string name("./vtk_surface.geom");
  check_open(name, fout);
  fout << "# vtk DataFile Version 2.0\n";
  fout << "file: " << name << std::endl;
  fout << "BINARY\n";
  fout << "DATASET UNSTRUCTURED_GRID\n";
  fout << "POINTS " << n << " float\n";
  fout.write((char *)&(coor[0]), sizeof(float) * 3 * n);
  fout << "CELLS " << ne << " " << ne * 5 << std::endl;
  fout.write((char *)&(conn[0]), sizeof(int) * 5 * ne);
  fout << "CELL_TYPES " << ne << std::endl;
  fout.write((char *)&(inttype[0]), sizeof(int) * ne);
  fout.close();

  name = "./vtk_surface.rank";
  check_open(name, fout);
  fout << "POINT_DATA " << n << "\n";
  fout << "SCALARS rank int 1\n";
  fout << "LOOKUP_TABLE default\n";
  fout.write((char *)&(rankid[0]), sizeof(int) * n);
  fout.close();

  name = "./vtk_surface.displshead";
  check_open(name, fout);
  fout << "POINT_DATA " << n << "\n";
  fout << "VECTORS displs float\n";
  fout.close();

#ifndef DEBUG_NO_MPI
  MPI_Finalize();
#endif
  return 0;
}

#else
int main(int argc, char **argv)
{
  int myrank, totalrank;
#ifndef DEBUG_NO_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &totalrank);
#else
  myrank = 0;
  totalrank = 0;
#endif

  double ds, gminx, gminy, gminz, gmaxx, gmaxy;
  int nmat, nproc, nthread;
  read_global_settings(myrank, ds, gminx, gminy, gminz, gmaxx, gmaxy, nmat, nproc, nthread);

  int n, ne_tet4, ne_vox8, neighbors;
  std::vector<double> coor;
  std::vector<int> conn_tet4;
  std::vector<int> conn_vox8;
  std::vector<int> neighborrankid;         // list of neighboring processes
  std::vector<int> mpinode_linear_pointer; // pointer to list of local nodes to be communicated
  std::vector<int> mpinode_linear_index;   // list of local nodes to be communicated
  //  mpinode_index[mpinode_linear_pointer[i]]->start of nodes to be communicated in mpinode_index for neighborrank i
  std::vector<int> surfnodes; // list of local nodes on surface
  std::vector<int> global_node_id; // optional global node id (1-based)
  std::vector<long> global_node_id64; // optional global node id (64-bit)

  read_local_data("./hdata/" + ToString(myrank, 6, '0') + ".data.h5",
                  n, ne_tet4, ne_vox8, neighbors,
                  coor, conn_tet4, conn_vox8,
                  neighborrankid, mpinode_linear_pointer, mpinode_linear_index, surfnodes,
                  global_node_id, global_node_id64);

  double dsorg(ds);
#ifdef MAKE_INFINITE_BOUNDARY
  if (myrank == 0)
    std::cout << "using INFINITE_BOUNDARY\n"
              << std::flush;
  bool flagx, flagy;
  for (int i = 0; i < coor.size() / 3; ++i)
  {
    flagx = (fabs(coor[3 * i]) < 0.01 * ds);
    flagy = (fabs(coor[3 * i + 1]) < 0.01 * ds);
    //    if(flagx) coor[3*i] += 0.0222*ds;
    //    if(flagy) coor[3*i+1] += 0.0222*ds;
  }
#endif

  // refine mesh
  // wrong->  //0 for tet4 nodes, 1 for edge nodes, (n+1) for level n refinement edge nodes
  std::vector<int> nodetype(coor.size() / 3, 0);
  for (int i = 1; i < REFINE_FACTOR; ++i)
  {
#ifdef OUTPUT_VTK
    std::vector<int> node_intscalar(coor.size() / 3, 0);
    std::vector<double> node_doublevec(coor.size(), 0.0);
    for (int j = 0; j < surfnodes.size(); ++j)
      node_intscalar[surfnodes[j]] = 1;
    for (int j = 0; j < mpinode_linear_index.size(); ++j)
      node_doublevec[3 * mpinode_linear_index[j]] = 1.0;
    for (int j = 0; j < nodetype.size(); ++j)
      node_doublevec[3 * j + 1] = nodetype[j];
    for (int j = 0; j < nodetype.size(); ++j)
      node_doublevec[3 * j + 2] = j;
    write_vtk_hybrid_mesh(coor, conn_tet4, conn_vox8, node_intscalar, node_doublevec,
                          "./hdata/hybridmesh_" + ToString(i) + "_" + ToString(myrank) + ".vtk");
#endif
    if (myrank == 0)
      std::cout << "refining mesh..." << i << "/" << REFINE_FACTOR - 1 << "\n"
                << std::flush;
    refine_linear_mesh(
        dsorg, gmaxx, gmaxy, coor, conn_vox8, conn_tet4, ds, myrank, totalrank,
        neighbors, neighborrankid, mpinode_linear_pointer, mpinode_linear_index, surfnodes, nodetype);
    ne_tet4 *= 8;
    ne_vox8 *= 8;
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
      std::cout << "refining mesh..." << i << "/" << REFINE_FACTOR - 1 << " finished\n"
                << std::flush;
  }
#ifdef OUTPUT_VTK
  std::vector<int> node_intscalar1(coor.size() / 3, 0);
  std::vector<double> node_doublevec1(coor.size(), 0.0);
  for (int i = 0; i < surfnodes.size(); ++i)
    node_intscalar1[surfnodes[i]] = 1;
  for (int i = 0; i < mpinode_linear_index.size(); ++i)
    node_doublevec1[3 * mpinode_linear_index[i]] = 1.0;
  for (int j = 0; j < nodetype.size(); ++j)
    node_doublevec1[3 * j + 1] = nodetype[j];
  for (int j = 0; j < nodetype.size(); ++j)
    node_doublevec1[3 * j + 2] = j;
  write_vtk_hybrid_mesh(coor, conn_tet4, conn_vox8, node_intscalar1, node_doublevec1,
                        "./hdata/hybridmesh_" + ToString(REFINE_FACTOR) + "_" + ToString(myrank) + ".vtk");
#endif

  // convert linear to quadratic
  std::vector<double> coor_quad;
  std::vector<int> conn_alltet10; // voxel elements are stored in beginning part, tetra elements in latter part
  std::vector<int> nodemap_quad2linear;
  std::vector<int> nodemap_linear2quad;
  bool merge4and10(false);
  convert_linear2quad(
      myrank, dsorg, gmaxx, gmaxy,
      merge4and10, coor, conn_vox8, conn_tet4, ds,
      coor_quad, conn_alltet10,
      nodemap_quad2linear, nodemap_linear2quad, nodetype);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "convert_linear2quad finished\n"
              << std::flush;

  std::vector<int> nodetype_linear(nodetype);
#ifdef RENUMBER_ELEMENTS // renumber elements for overlapping communication and computation
  renumber_elements(mpinode_linear_index, nodemap_quad2linear, nodemap_linear2quad,
                    conn_vox8, conn_tet4, conn_alltet10, myrank);
#endif
#ifdef RENUMBER_NODES
  renumber_nodes(coor_quad, nodetype, conn_alltet10, nodemap_quad2linear, nodemap_linear2quad);
#endif

  // mpi_settings
  std::vector<int> mpinode_quad_pointer;
  std::vector<int> mpinode_quad_index;
#ifndef DEBUG_NO_COMM
  MPI_set_comm_linear_settings(coor, ds, neighborrankid,
                               mpinode_linear_pointer, mpinode_linear_index,
                               myrank, totalrank, nodetype_linear);

  // -- check mpi settings - linear
  MPI_check_comm_settings(coor, ds, neighborrankid,
                          mpinode_linear_pointer, mpinode_linear_index, myrank, totalrank, nodetype_linear);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "MPI_set_comm_linear_settings finished\n"
              << std::flush;

  MPI_set_comm_quad_settings(coor, coor_quad, ds, conn_alltet10, neighborrankid,
                             mpinode_linear_pointer, mpinode_linear_index,
                             nodemap_quad2linear, nodemap_linear2quad,
                             mpinode_quad_pointer, mpinode_quad_index, myrank, totalrank, nodetype);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "MPI_set_comm_quad_settings finished\n"
              << std::flush;

#else
  for (int i = 0; i < neighborrankid.size() + 1; ++i)
    ;
  mpinode_quad_pointer.push_back(0);
#endif

  // set boundary
  std::vector<int> boundarynode_linear;
  std::vector<int> boundarynode_quad;
  set_boundary(coor, ds, gminx, gminy, gminz, gmaxx, gmaxy, boundarynode_linear);
  set_boundary(coor_quad, ds, gminx, gminy, gminz, gmaxx, gmaxy, boundarynode_quad);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "set_boundary finished\n"
              << std::flush;

#ifdef STRUCTURED_MESH_SET_MATERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "structured_mesh_set_material\n"
              << std::flush;
  modify_mat(myrank, nmat, coor, conn_tet4, conn_vox8, conn_alltet10);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "structured_mesh_set_material finished\n"
              << std::flush;
#endif

#ifndef DEBUG_NO_HDF5
  long fidc(1), fidm(2);
  hdf_create_file_(("./cdata/" + ToString(myrank, 6, '0') + ".data.h5").c_str(), &fidc);
  hdf_create_file_(("./mdata/" + ToString(myrank, 6, '0') + ".data.h5").c_str(), &fidm);

#ifdef MAKE_INFINITE_BOUNDARY
  convert_infinite_boundary_coor(myrank, dsorg, gmaxx, gmaxy, conn_alltet10, nodemap_linear2quad, coor_quad, coor);
#endif

  // surface_nodes
  std::vector<int> surfnodes_quad;
  MPI_extract_surface_nodes_quad(surfnodes, nodemap_linear2quad, ds, coor_quad, conn_alltet10, surfnodes_quad,
                                 myrank, neighborrankid, mpinode_quad_pointer, mpinode_quad_index);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "MPI_extract_surface_nodes_quad finished\n"
              << std::flush;

#ifdef OUTPUT_3DRESULT
  double minelev(-30);
  MPI_extract_3D_nodes_quad(coor_quad, conn_alltet10, surfnodes_quad, myrank, minelev);
#endif

  std::vector<int> global_node_id_quad;
  std::vector<long> global_node_id64_quad;
  if (!global_node_id64.empty() || !global_node_id.empty()) {
    if (global_node_id64.empty() && !global_node_id.empty()) {
      global_node_id64.assign(global_node_id.begin(), global_node_id.end());
    }
    if (global_node_id.empty() && !global_node_id64.empty()) {
      global_node_id.resize(global_node_id64.size(), 0);
      for (int i = 0; i < static_cast<int>(global_node_id64.size()); ++i) {
        if (global_node_id64[i] > 0 && global_node_id64[i] <= INT_MAX)
          global_node_id[i] = static_cast<int>(global_node_id64[i]);
      }
    }

    global_node_id_quad.assign(coor_quad.size() / 3, 0);
    global_node_id64_quad.assign(coor_quad.size() / 3, 0);
    for (int inode = 0; inode < static_cast<int>(global_node_id64.size()); ++inode) {
      int qnode = nodemap_linear2quad[inode];
      if (qnode >= 0 && qnode < static_cast<int>(global_node_id_quad.size()))
        global_node_id_quad[qnode] = (inode < static_cast<int>(global_node_id.size())) ? global_node_id[inode] : 0;
      if (qnode >= 0 && qnode < static_cast<int>(global_node_id64_quad.size()))
        global_node_id64_quad[qnode] = global_node_id64[inode];
    }

    auto edge_gid64 = [](long a, long b) -> long {
      if (a > b) {
        long t = a;
        a = b;
        b = t;
      }
      unsigned long long ua = static_cast<unsigned long long>(a);
      unsigned long long ub = static_cast<unsigned long long>(b);
      unsigned long long key = (ua << 32) | ub;
      return static_cast<long>(key);
    };

    int ne = conn_alltet10.size() / 11;
    for (int ie = 0; ie < ne; ++ie) {
      int n0 = conn_alltet10[11 * ie + 0];
      int n1 = conn_alltet10[11 * ie + 1];
      int n2 = conn_alltet10[11 * ie + 2];
      int n3 = conn_alltet10[11 * ie + 3];
      int e01 = conn_alltet10[11 * ie + 4];
      int e12 = conn_alltet10[11 * ie + 5];
      int e02 = conn_alltet10[11 * ie + 6];
      int e03 = conn_alltet10[11 * ie + 7];
      int e13 = conn_alltet10[11 * ie + 8];
      int e23 = conn_alltet10[11 * ie + 9];

      long g0 = global_node_id64_quad[n0];
      long g1 = global_node_id64_quad[n1];
      long g2 = global_node_id64_quad[n2];
      long g3 = global_node_id64_quad[n3];

      if (e01 >= 0 && e01 < static_cast<int>(global_node_id64_quad.size()) && global_node_id64_quad[e01] == 0)
        global_node_id64_quad[e01] = edge_gid64(g0, g1);
      if (e12 >= 0 && e12 < static_cast<int>(global_node_id64_quad.size()) && global_node_id64_quad[e12] == 0)
        global_node_id64_quad[e12] = edge_gid64(g1, g2);
      if (e02 >= 0 && e02 < static_cast<int>(global_node_id64_quad.size()) && global_node_id64_quad[e02] == 0)
        global_node_id64_quad[e02] = edge_gid64(g0, g2);
      if (e03 >= 0 && e03 < static_cast<int>(global_node_id64_quad.size()) && global_node_id64_quad[e03] == 0)
        global_node_id64_quad[e03] = edge_gid64(g0, g3);
      if (e13 >= 0 && e13 < static_cast<int>(global_node_id64_quad.size()) && global_node_id64_quad[e13] == 0)
        global_node_id64_quad[e13] = edge_gid64(g1, g3);
      if (e23 >= 0 && e23 < static_cast<int>(global_node_id64_quad.size()) && global_node_id64_quad[e23] == 0)
        global_node_id64_quad[e23] = edge_gid64(g2, g3);
    }
  }

  write_hdf_mdata(coor_quad, conn_alltet10, neighborrankid,
                  mpinode_quad_pointer, mpinode_quad_index, boundarynode_quad, nmat, fidm,
                  global_node_id_quad, global_node_id64_quad);
  write_hdf_numofelements(fidm, conn_alltet10.size() / 11, nthread);

  write_hdf_cdata(coor, conn_vox8, conn_tet4, neighborrankid,
                  mpinode_linear_pointer, mpinode_linear_index, boundarynode_linear, nmat, fidc,
                  global_node_id, global_node_id64);

  write_settings(coor_quad, conn_alltet10, neighborrankid, mpinode_quad_pointer, mpinode_quad_index, nmat, myrank);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "write_hdf finished\n"
              << std::flush;

  // set mapping -- call after write_hdf_mdata
  int n_quad(coor_quad.size() / 3);
  int n_linear(coor.size() / 3);
  int ne(conn_alltet10.size() / 11);
  int ne_linear(ne);
  mapping_pairs_lib_(&n_quad, &ne, &n_linear, &ne_linear, &fidm, &fidc);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "here1\n"
              << std::flush;

  // set ABC linear -- call after write_hdf_cdata
  int ib_linear(boundarynode_linear.size());
  abc4_lib_(&n_linear, &ne, &nmat, &ib_linear, &fidc);
#ifdef USECRS4
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "USECRS4 1\n"
              << std::flush;
  crs_nodes4_lib_(&n_linear, &ne, &fidc);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "USECRS4 2\n"
              << std::flush;
  crs_connect_block4_lib_(&n_linear, &ne, &ib_linear, &fidc);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "here2\n"
              << std::flush;

  // set ABC quad -- call after write_hdf_mdata
  long n_quad_long(n_quad);
  long ne_long(ne);
  long nmat_long(nmat);
  long bufsl(20 * n_quad);
  abc2_lib3_(&n_quad_long, &ne_long, &nmat_long, &fidm, &bufsl);
#ifdef USECRS10
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "USECRS10 1\n"
              << std::flush;
  crs_nodes_lib_(&n_quad, &ne, &fidm);
  int ib_;
  long ibs_(0), ibc_(1);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "USECRS10 2\n"
              << std::flush;
  hdf_read_int_array_part_(&fidm, "/BC", &ibs_, &ibc_, &ib_);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "USECRS10 3\n"
              << std::flush;
  crs_connect_block_lib_(&n_quad, &ne, &ib_, &fidm);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "here3\n"
              << std::flush;

  hdf_close_file_(&fidc);
  hdf_close_file_(&fidm);
#endif

#ifdef OUTPUT_VTK
  std::vector<int> node_intscalar(coor_quad.size() / 3, 0);
  std::vector<double> node_doublevec(coor_quad.size(), 0.0);
  for (int i = 0; i < surfnodes_quad.size(); ++i)
    node_intscalar[surfnodes_quad[i]] = 1;
  for (int i = 0; i < mpinode_quad_index.size(); ++i)
    node_doublevec[3 * mpinode_quad_index[i]] = 1.0;
  for (int j = 0; j < nodetype.size(); ++j)
    node_doublevec[3 * j + 1] = nodetype[j];
  for (int j = 0; j < nodetype.size(); ++j)
    node_doublevec[3 * j + 2] = j;
  write_vtk_tet10_mesh(coor_quad, conn_alltet10, node_intscalar, node_doublevec,
                       "./hdata/" + ToString(myrank, 6, '0') + ".part_quad_nodetype.vtk");
  convert_quad2linear(coor_quad, conn_alltet10, ne_tet4, ne_vox8, coor, conn_vox8, conn_tet4);
  write_vtk_hybrid_mesh(coor, conn_tet4, conn_vox8, node_intscalar, node_doublevec,
                        "./hdata/hybridmesh_" + ToString(REFINE_FACTOR + 1) + "_" + ToString(myrank) + ".vtk");
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "Finished\n"
              << std::flush;

#ifndef DEBUG_NO_MPI
  MPI_Finalize();
#endif
  return 0;
}
#endif

void modify_mat(int rank, int nummat, const std::vector<double> &coor4, std::vector<int> &conn_tet4, std::vector<int> &conn_vox8, std::vector<int> &conn_tet10)
{
  std::vector<std::vector<double>> range(3, std::vector<double>(2)); // x,y,z, max,min

  long numnodes, numelems, numelemsvox;
  numnodes = coor4.size() / 3;
  numelems = conn_tet4.size() / 5;
  numelemsvox = conn_vox8.size() / 9;

  // set min max range
  for (long j = 0; j < 3; ++j)
  {
    range[j][0] = coor4[j];
    range[j][1] = coor4[j];
  }
  for (long i = 0; i < numnodes; ++i)
  {
    for (long j = 0; j < 3; ++j)
    {
      range[j][0] = std::min<double>(range[j][0], coor4[3 * i + j]);
      range[j][1] = std::max<double>(range[j][1], coor4[3 * i + j]);
    }
  }

  // set delta x at surface (ds)
  // set delta x at bottom (ds_bot)
  double ds(range[0][1]);
  double ds_sdata;
  double dtmp;
  for (long i = 0; i < numnodes; ++i)
  {
    if (range[0][0] + 0.00001 > coor4[3 * i])
      continue;
    ds = std::min<double>(ds, coor4[3 * i] - range[0][0]);
  }
  MPI_Allreduce(&ds, &dtmp, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  ds = dtmp;

  MPI_Allreduce(&range[0][0], &dtmp, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  range[0][0] = dtmp;
  MPI_Allreduce(&range[1][0], &dtmp, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  range[1][0] = dtmp;
  MPI_Allreduce(&range[2][0], &dtmp, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  range[2][0] = dtmp;
  MPI_Allreduce(&range[0][1], &dtmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  range[0][1] = dtmp;
  MPI_Allreduce(&range[1][1], &dtmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  range[1][1] = dtmp;
  MPI_Allreduce(&range[2][1], &dtmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  range[2][1] = dtmp;

  long ix, iy, iz, ltmp;
  std::string temp;
  int nx_sdata, ny_sdata, nummat_sdata;
  std::ifstream fin;
  if (rank == 0)
  {
    check_open("./sdata/modeling_setting.dat", fin);
    fin >> temp; // nx
    fin >> temp; // ny
    fin >> nx_sdata;
    fin >> ny_sdata;
    fin >> temp; // ds
    fin >> ds_sdata;
    fin >> temp; // num
    fin >> temp; // of
    fin >> temp; // layer
    fin >> nummat_sdata;
    fin.close();
    nx_sdata++;
    ny_sdata++;
    if (nummat_sdata != nummat)
      Error("number of layers does not match");

    std::cout << "setting material for structured mesh\n";
    std::cout << "ds: " << ds << "\n";
    std::cout << "ds_sdata: " << ds_sdata << "\n";
    std::cout << "nx_sdata: " << nx_sdata << "\n";
    std::cout << "ny_sdata: " << ny_sdata << "\n";
  }
  MPI_Bcast(&ds_sdata, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nx_sdata, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ny_sdata, 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<std::vector<double>> demelev(nummat, std::vector<double>(nx_sdata * ny_sdata, 0.0));
  double dsrat(ds / ds_sdata);
  int maxrat(ds_sdata / ds);
  if (rank == 0)
  {
    std::cout << "dsrat: " << dsrat << "\n";
    std::cout << "invdsrat: " << maxrat << "\n";
    if (dsrat > 1.01)
      Error("dsrat too big");
    if (fabs(1.0 / dsrat - maxrat) > 0.0001)
      Error("maxrat");

    for (long i = 0; i < nummat; ++i)
    {
      check_open("./sdata/sur" + ToString(i + 1, 4, '0') + ".dat", fin);
      fin >> temp; // nx,
      fin >> temp; // ny
      fin >> ltmp; // nx
      if (ltmp + 1 != nx_sdata)
        Error("something wrong with nx");
      fin >> ltmp; // ny
      if (ltmp + 1 != ny_sdata)
        Error("something wrong with ny");
      fin >> temp; // DEM
      fin >> temp; // data
      for (long j = 0; j < ny_sdata; ++j)
      {
        for (long k = 0; k < nx_sdata; ++k)
        {
          fin >> demelev[i][j * nx_sdata + k];
        }
      }
      fin.close();
    }
  }
  for (long i = 0; i < nummat; ++i)
    MPI_Bcast(&(demelev[i][0]), nx_sdata * ny_sdata, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  double x, y, z;
  long j, iix, iiy, iiz;

  for (long i = 0; i < numelems; ++i)
  {
    x = y = z = 0.0;
    for (j = 0; j < 4; ++j)
    {
      x += coor4[3 * conn_tet4[5 * i + j] + 0];
      y += coor4[3 * conn_tet4[5 * i + j] + 1];
      z += coor4[3 * conn_tet4[5 * i + j] + 2];
    }
    x *= 0.25;
    y *= 0.25;
    z *= 0.25;
    ix = (x - range[0][0]) / ds_sdata;
    iy = (y - range[1][0]) / ds_sdata;
    iz = (z - range[2][0]) / ds_sdata;
    iix = (x - range[0][0] - ix * ds_sdata) / ds;
    iiy = (y - range[1][0] - iy * ds_sdata) / ds;
    iiz = (z - range[2][0] - iz * ds_sdata) / ds;

    if (ix < 0 || ix >= nx_sdata - 1)
      Error("ix,nx,wrong " + ToString(ix) + " " + ToString(nx_sdata));
    if (iy < 0 || iy >= ny_sdata - 1)
      Error("iy,ny,wrong " + ToString(iy) + " " + ToString(ny_sdata));
    if (iix < 0 || iix >= maxrat)
      Error("iix,maxrat,wrong " + ToString(iix) + " " + ToString(maxrat));
    if (iiy < 0 || iiy >= maxrat)
      Error("iiy,maxrat,wrong " + ToString(iiy) + " " + ToString(maxrat));
    if (iiz < 0 || iiz >= maxrat)
      Error("iiz,maxrat,wrong " + ToString(iiz) + " " + ToString(maxrat));

    if (ix * ds_sdata + iix * ds > x)
      Error("x computation wrong 3? " + ToString(ix) + " " + ToString(iix) + " " + ToString(x));
    if (iy * ds_sdata + iiy * ds > y)
      Error("y computation wrong 3? " + ToString(iy) + " " + ToString(iiy) + " " + ToString(y));
    if (iz * ds_sdata + iiz * ds > z)
      Error("z computation wrong 3? " + ToString(iz) + " " + ToString(iiz) + " " + ToString(z));

    /*
        if( ix*ds_sdata+(iix+0.5)*ds > x )
          Error("x computation wrong 1? "+ToString(ix)+" "+ToString(iix)+" "+ToString(x));
        if( ix*ds_sdata+(iix+0.5)*ds < x-0.5*ds )
          Error("x computation wrong 2? "+ToString(ix)+" "+ToString(iix)+" "+ToString(x));
        if( iy*ds_sdata+(iiy+0.5)*ds > y )
          Error("y computation wrong 1? "+ToString(iy)+" "+ToString(iiy)+" "+ToString(y));
        if( iy*ds_sdata+(iiy+0.5)*ds < y-0.5*ds )
          Error("y computation wrong 2? "+ToString(iy)+" "+ToString(iiy)+" "+ToString(y));
        if( iz*ds_sdata+(iiz+0.5)*ds > z )
          Error("z computation wrong 1? "+ToString(iz)+" "+ToString(iiz)+" "+ToString(z));
        if( iz*ds_sdata+(iiz+0.5)*ds < z-0.5*ds )
          Error("z computation wrong 2? "+ToString(iz)+" "+ToString(iiz)+" "+ToString(z));
    */
    conn_tet4[5 * i + 4] = conn_tet10[11 * (6 * numelemsvox + i) + 10] = 0;
    for (j = nummat - 1; j >= 0; --j)
    {
      if (iz * ds_sdata + (iiz + 0.5) * ds < (demelev[j][iy * nx_sdata + ix] * (1.0 - (iix + 0.5) * dsrat) * (1.0 - (iiy + 0.5) * dsrat) +
                                              demelev[j][(iy + 1) * nx_sdata + ix] * (1.0 - (iix + 0.5) * dsrat) * ((iiy + 0.5) * dsrat) +
                                              demelev[j][(iy + 1) * nx_sdata + (ix + 1)] * ((iix + 0.5) * dsrat) * ((iiy + 0.5) * dsrat) +
                                              demelev[j][iy * nx_sdata + (ix + 1)] * ((iix + 0.5) * dsrat) * (1.0 - (iiy + 0.5) * dsrat)))
      {
        conn_tet4[5 * i + 4] = conn_tet10[11 * (6 * numelemsvox + i) + 10] = j;
        break;
      }
    }
    /*
        if( rank==0 && fabs(x-244500)<1000.0 && fabs(y-384500)<1000.0 && fabs(z - 400000)<1000.0 )
        {
          std::cout << "write elem number " << i << " " << conn_tet4[5*i+4] << "\n";
          for( int kk = 0; kk < 4; ++kk )
            std::cout << "node " << kk << " " << coor4[3*conn_tet4[5*i+kk]+0] << " " << coor4[3*conn_tet4[5*i+kk]+1] << " " << coor4[3*conn_tet4[5*i+kk]+2] << "\n";
          std::cout << "x " << ix << " " << iix << " " << x << "\n";
          std::cout << "y " << iy << " " << iiy << " " << y << "\n";
          std::cout << "z " << iz << " " << iiz << " " << z << "\n";
        }
    */
  }
  long kk;
  for (long i = 0; i < numelemsvox; ++i)
  {
    x = y = z = 0.0;
    for (j = 0; j < 8; ++j)
    {
      x += coor4[3 * conn_vox8[9 * i + j] + 0];
      y += coor4[3 * conn_vox8[9 * i + j] + 1];
      z += coor4[3 * conn_vox8[9 * i + j] + 2];
    }
    x *= 0.125;
    y *= 0.125;
    z *= 0.125;
    ix = (x - range[0][0]) / ds_sdata;
    iy = (y - range[1][0]) / ds_sdata;
    iz = (z - range[2][0]) / ds_sdata;
    iix = (x - range[0][0] - ix * ds_sdata) / ds;
    iiy = (y - range[1][0] - iy * ds_sdata) / ds;
    iiz = (z - range[2][0] - iz * ds_sdata) / ds;

    if (ix < 0 || ix >= nx_sdata - 1)
      Error("ix,nx,wrong " + ToString(ix) + " " + ToString(nx_sdata));
    if (iy < 0 || iy >= ny_sdata - 1)
      Error("iy,ny,wrong " + ToString(iy) + " " + ToString(ny_sdata));
    if (iix < 0 || iix >= maxrat)
      Error("iix,maxrat,wrong " + ToString(iix) + " " + ToString(maxrat));
    if (iiy < 0 || iiy >= maxrat)
      Error("iiy,maxrat,wrong " + ToString(iiy) + " " + ToString(maxrat));
    if (iiz < 0 || iiz >= maxrat)
      Error("iiz,maxrat,wrong " + ToString(iiz) + " " + ToString(maxrat));

    if (ix * ds_sdata + iix * ds > x)
      Error("x computation wrong 3? " + ToString(ix) + " " + ToString(iix) + " " + ToString(x));
    if (iy * ds_sdata + iiy * ds > y)
      Error("y computation wrong 3? " + ToString(iy) + " " + ToString(iiy) + " " + ToString(y));
    if (iz * ds_sdata + iiz * ds > z)
      Error("z computation wrong 3? " + ToString(iz) + " " + ToString(iiz) + " " + ToString(z));

    /*
        if( ix*ds_sdata+(iix+0.5)*ds > x )
          Error("x computation wrong 1? "+ToString(ix)+" "+ToString(iix)+" "+ToString(x));
        if( ix*ds_sdata+(iix+0.5)*ds < x-0.5*ds )
          Error("x computation wrong 2? "+ToString(ix)+" "+ToString(iix)+" "+ToString(x));
        if( iy*ds_sdata+(iiy+0.5)*ds > y )
          Error("y computation wrong 1? "+ToString(iy)+" "+ToString(iiy)+" "+ToString(y));
        if( iy*ds_sdata+(iiy+0.5)*ds < y-0.5*ds )
          Error("y computation wrong 2? "+ToString(iy)+" "+ToString(iiy)+" "+ToString(y));
        if( iz*ds_sdata+(iiz+0.5)*ds > z )
          Error("z computation wrong 1? "+ToString(iz)+" "+ToString(iiz)+" "+ToString(z));
        if( iz*ds_sdata+(iiz+0.5)*ds < z-0.5*ds )
          Error("z computation wrong 2? "+ToString(iz)+" "+ToString(iiz)+" "+ToString(z));
    */
    conn_vox8[9 * i + 8] = 0;
    for (j = 0; j < 6; ++j)
      conn_tet10[11 * (6 * i + j) + 10] = 0;
    for (j = nummat - 1; j >= 0; --j)
    {
      if (iz * ds_sdata + (iiz + 0.5) * ds < (demelev[j][iy * nx_sdata + ix] * (1.0 - (iix + 0.5) * dsrat) * (1.0 - (iiy + 0.5) * dsrat) +
                                              demelev[j][(iy + 1) * nx_sdata + ix] * (1.0 - (iix + 0.5) * dsrat) * ((iiy + 0.5) * dsrat) +
                                              demelev[j][(iy + 1) * nx_sdata + (ix + 1)] * ((iix + 0.5) * dsrat) * ((iiy + 0.5) * dsrat) +
                                              demelev[j][iy * nx_sdata + (ix + 1)] * ((iix + 0.5) * dsrat) * (1.0 - (iiy + 0.5) * dsrat)))
      {
        conn_vox8[9 * i + 8] = j;
        for (kk = 0; kk < 6; ++kk)
          conn_tet10[11 * (6 * i + kk) + 10] = j;
        break;
      }
    }
    //    conn_vox8[9*i+8] = conn_tet10[11*i+10] = -1;
    /*
        if( rank==0 && fabs(x-244500)<1000.0 && fabs(y-384500)<1000.0 && fabs(z - 400000)<1000.0 )
        {
          std::cout << "write elem number " << i << " " << conn_tet4[5*i+4] << "\n";
          for( int kk = 0; kk < 4; ++kk )
            std::cout << "node " << kk << " " << coor4[3*conn_tet4[5*i+kk]+0] << " " << coor4[3*conn_tet4[5*i+kk]+1] << " " << coor4[3*conn_tet4[5*i+kk]+2] << "\n";
          std::cout << "x " << ix << " " << iix << " " << x << "\n";
          std::cout << "y " << iy << " " << iiy << " " << y << "\n";
          std::cout << "z " << iz << " " << iiz << " " << z << "\n";
        }
    */
  }
}

void convert_infinite_boundary_coor(int myrank, double dsorg, double gmaxx, double gmaxy, const std::vector<int> &conn_alltet10, const std::vector<int> &nodemap_linear2quad, std::vector<double> &coor_quad, std::vector<double> &coor)
{
  if (myrank == 0)
    std::cout << "using INFINITE_BOUNDARY\n";
  double meanx(0.0), meany(0.0);
  for (int i = 0; i < coor_quad.size() / 3; ++i)
  {
    //    if( fabs(coor_quad[3*i  ]-0.0222*dsorg) < 0.01*dsorg ) coor_quad[3*i] = 0.0;
    //    if( fabs(coor_quad[3*i+1]-0.0222*dsorg) < 0.01*dsorg ) coor_quad[3*i+1] = 0.0;
    meanx += coor_quad[3 * i];
    meany += coor_quad[3 * i + 1];
  }
  meanx /= (coor_quad.size() / 3);
  meany /= (coor_quad.size() / 3);

  int iii, index;
  double elemrange[2][2];
  std::vector<bool> flagclampx(coor_quad.size() / 3, false);
  std::vector<bool> flagclampy(coor_quad.size() / 3, false);
  for (int i = 0; i < conn_alltet10.size() / 11; ++i)
  {
    index = conn_alltet10[11 * i];
    elemrange[0][0] = elemrange[1][0] = coor_quad[3 * index];
    elemrange[0][1] = elemrange[1][1] = coor_quad[3 * index + 1];
    for (iii = 1; iii < 4; ++iii)
    {
      index = conn_alltet10[11 * i + iii];
      elemrange[0][0] = std::min(elemrange[0][0], coor_quad[3 * index]);
      elemrange[1][0] = std::max(elemrange[1][0], coor_quad[3 * index]);
      elemrange[0][1] = std::min(elemrange[0][1], coor_quad[3 * index + 1]);
      elemrange[1][1] = std::max(elemrange[1][1], coor_quad[3 * index + 1]);
    }
    if (elemrange[1][0] - elemrange[0][0] < gmaxx * 0.8)
      for (iii = 0; iii < 4; ++iii)
        flagclampx[conn_alltet10[11 * i + iii]] = true;
    if (elemrange[1][1] - elemrange[0][1] < gmaxy * 0.8)
      for (iii = 0; iii < 4; ++iii)
        flagclampy[conn_alltet10[11 * i + iii]] = true;
  }
  int iiix(0), iiiy(0);
  meanx = meany = 0.0;
  for (int i = 0; i < coor_quad.size() / 3; ++i)
  {
    if (flagclampx[i])
    {
      meanx += coor_quad[3 * i];
      ++iiix;
    }
    if (flagclampy[i])
    {
      meany += coor_quad[3 * i + 1];
      ++iiiy;
    }
  }
  meanx /= iiix;
  meany /= iiiy;
  for (int i = 0; i < coor_quad.size() / 3; ++i)
  {
    if (!flagclampx[i])
    {
      if (fabs(coor_quad[3 * i]) < 0.01 * dsorg && meanx > 0.5 * gmaxx)
        coor_quad[3 * i] += gmaxx;
      else if (fabs(coor_quad[3 * i]) > 0.01 * dsorg && meanx < 0.5 * gmaxx)
        coor_quad[3 * i] -= gmaxx;
    }
    if (!flagclampy[i])
    {
      if (fabs(coor_quad[3 * i + 1]) < 0.01 * dsorg && meany > 0.5 * gmaxy)
        coor_quad[3 * i + 1] += gmaxy;
      else if (fabs(coor_quad[3 * i + 1]) > 0.01 * dsorg && meany < 0.5 * gmaxy)
        coor_quad[3 * i + 1] -= gmaxy;
    }
  }
  double maxrangex(0.0), maxrangey(0.0);
  for (int i = 0; i < conn_alltet10.size() / 11; ++i)
  {
    index = conn_alltet10[11 * i];
    elemrange[0][0] = elemrange[1][0] = coor_quad[3 * index];
    elemrange[0][1] = elemrange[1][1] = coor_quad[3 * index + 1];
    for (iii = 1; iii < 4; ++iii)
    {
      index = conn_alltet10[11 * i + iii];
      elemrange[0][0] = std::min(elemrange[0][0], coor_quad[3 * index]);
      elemrange[1][0] = std::max(elemrange[1][0], coor_quad[3 * index]);
      elemrange[0][1] = std::min(elemrange[0][1], coor_quad[3 * index + 1]);
      elemrange[1][1] = std::max(elemrange[1][1], coor_quad[3 * index + 1]);
    }
    maxrangex = std::max(elemrange[1][0] - elemrange[0][0], maxrangex);
    maxrangey = std::max(elemrange[1][1] - elemrange[0][1], maxrangey);
  }
  if (maxrangex > 0.5 * gmaxx)
    std::cout << "ERROR: maxrangex " << myrank << " " << maxrangex << "\n"
              << std::flush;
  if (maxrangey > 0.5 * gmaxy)
    std::cout << "ERROR: maxrangey " << myrank << " " << maxrangey << "\n"
              << std::flush;

  for (int i = 0; i < conn_alltet10.size() / 11; ++i)
  {
    for (iii = 0; iii < 3; ++iii)
    {
      coor_quad[3 * conn_alltet10[11 * i + 4] + iii] = 0.5 * (coor_quad[3 * conn_alltet10[11 * i + 0] + iii] +
                                                              coor_quad[3 * conn_alltet10[11 * i + 1] + iii]);

      coor_quad[3 * conn_alltet10[11 * i + 5] + iii] = 0.5 * (coor_quad[3 * conn_alltet10[11 * i + 1] + iii] +
                                                              coor_quad[3 * conn_alltet10[11 * i + 2] + iii]);

      coor_quad[3 * conn_alltet10[11 * i + 6] + iii] = 0.5 * (coor_quad[3 * conn_alltet10[11 * i + 0] + iii] +
                                                              coor_quad[3 * conn_alltet10[11 * i + 2] + iii]);

      coor_quad[3 * conn_alltet10[11 * i + 7] + iii] = 0.5 * (coor_quad[3 * conn_alltet10[11 * i + 0] + iii] +
                                                              coor_quad[3 * conn_alltet10[11 * i + 3] + iii]);

      coor_quad[3 * conn_alltet10[11 * i + 8] + iii] = 0.5 * (coor_quad[3 * conn_alltet10[11 * i + 1] + iii] +
                                                              coor_quad[3 * conn_alltet10[11 * i + 3] + iii]);

      coor_quad[3 * conn_alltet10[11 * i + 9] + iii] = 0.5 * (coor_quad[3 * conn_alltet10[11 * i + 2] + iii] +
                                                              coor_quad[3 * conn_alltet10[11 * i + 3] + iii]);
    }
  }
  for (int i = 0; i < coor.size() / 3; ++i)
  {
    coor[3 * i + 0] = coor_quad[3 * nodemap_linear2quad[i] + 0];
    coor[3 * i + 1] = coor_quad[3 * nodemap_linear2quad[i] + 1];
    coor[3 * i + 2] = coor_quad[3 * nodemap_linear2quad[i] + 2];
  }
}
