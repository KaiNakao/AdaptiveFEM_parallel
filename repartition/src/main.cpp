#include <mpi.h>
#include <parmetis.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <numeric>
#include <cmath>

#include "hdf_lib.h"

namespace {

enum class PartitionMethod {
  AdaptiveRepart,
  PartMeshKway
};

struct FaceKey {
  std::array<int, 3> idx;

  bool operator<(const FaceKey &other) const {
    return idx < other.idx;
  }
};

std::string RankFilename(const std::string &dir, int rank, const std::string &suffix) {
  std::string rank_id = std::to_string(rank);
  rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
  return dir + "/" + rank_id + suffix;
}

void ReadSetting(int myid, int &nnode, int &nelem, int &nmaterial) {
  long fid = 0;
  const std::string filename = RankFilename("rdata", myid, ".data.h5");
  hdf_open_file_(filename.c_str(), &fid);
  long size = 3;
  int settings[3];
  long count = 0;
  hdf_read_int_array_(&fid, "/setting", settings, &size, &count);
  nnode = settings[0];
  nelem = settings[1];
  nmaterial = settings[2];
  hdf_close_file_(&fid);
}

void ReadConn(int myid, int nelem,
              std::vector<std::array<int, 4>> &conn,
              std::vector<int> &matid) {
  long fid = 0;
  const std::string filename = RankFilename("rdata", myid, ".data.h5");
  hdf_open_file_(filename.c_str(), &fid);

  long size = static_cast<long>(5) * nelem;
  long count = 0;
  std::vector<int> conn_buf(size);
  hdf_read_int_array_(&fid, "/conn", conn_buf.data(), &size, &count);

  conn.resize(nelem);
  matid.resize(nelem);
  for (int e = 0; e < nelem; ++e) {
    for (int i = 0; i < 4; ++i) {
      conn[e][i] = conn_buf[5 * e + i] - 1;  // to 0-based
    }
    matid[e] = conn_buf[5 * e + 4] - 1;
  }

  hdf_close_file_(&fid);
}

void ReadMPInode(int myid, std::map<int, std::vector<int>> &neighbor_nodes) {
  long fid = 0;
  const std::string filename = RankFilename("rdata", myid, ".data.h5");
  hdf_open_file_(filename.c_str(), &fid);

  if (!hdf_has_dataset_(&fid, "/MPInode")) {
    hdf_close_file_(&fid);
    return;
  }

  long size = 0;
  hdf_get_array_size_(&fid, "/MPInode", &size);
  long count = 0;
  std::vector<int> buf(static_cast<size_t>(size));
  hdf_read_int_array_(&fid, "/MPInode", buf.data(), &size, &count);
  hdf_close_file_(&fid);

  int idx = 0;
  const int nneighbor = buf[idx++];
  std::vector<int> neighbor_rank(nneighbor);
  std::vector<int> neighbor_size(nneighbor);
  for (int i = 0; i < nneighbor; ++i) {
    neighbor_rank[i] = buf[idx++];
    neighbor_size[i] = buf[idx++];
  }

  for (int i = 0; i < nneighbor; ++i) {
    neighbor_nodes[neighbor_rank[i]].resize(neighbor_size[i]);
    for (int j = 0; j < neighbor_size[i]; ++j) {
      neighbor_nodes[neighbor_rank[i]][j] = buf[idx++] - 1;  // to 0-based
    }
  }
}

void ReadCoor(int myid, int nnode, std::vector<std::array<double, 3>> &coor) {
  long fid = 0;
  const std::string filename = RankFilename("rdata", myid, ".data.h5");
  hdf_open_file_(filename.c_str(), &fid);

  long size = static_cast<long>(3) * nnode;
  long count = 0;
  std::vector<double> buf(size);
  hdf_read_double_array_(&fid, "/coor", buf.data(), &size, &count);
  hdf_close_file_(&fid);

  coor.resize(nnode);
  for (int i = 0; i < nnode; ++i) {
    coor[i][0] = buf[3 * i + 0];
    coor[i][1] = buf[3 * i + 1];
    coor[i][2] = buf[3 * i + 2];
  }
}

void AddAdj(std::vector<std::unordered_set<idx_t>> &adj,
            idx_t local_elem, idx_t global_neighbor) {
  adj[local_elem].insert(global_neighbor);
}

void BuildLocalAdjacency(const std::vector<std::array<int, 4>> &conn,
                         idx_t base_gid,
                         std::vector<std::unordered_set<idx_t>> &adj) {
  std::map<FaceKey, std::vector<int>> face_to_elem;
  for (int e = 0; e < static_cast<int>(conn.size()); ++e) {
    const auto &c = conn[e];
    const int faces[4][3] = {
        {c[1], c[2], c[3]},
        {c[0], c[2], c[3]},
        {c[0], c[1], c[3]},
        {c[0], c[1], c[2]},
    };
    for (int f = 0; f < 4; ++f) {
      std::array<int, 3> nodes = {faces[f][0], faces[f][1], faces[f][2]};
      std::sort(nodes.begin(), nodes.end());
      face_to_elem[FaceKey{nodes}].push_back(e);
    }
  }

  for (const auto &kv : face_to_elem) {
    if (kv.second.size() != 2) {
      continue;
    }
    int e0 = kv.second[0];
    int e1 = kv.second[1];
    AddAdj(adj, e0, base_gid + e1);
    AddAdj(adj, e1, base_gid + e0);
  }
}

void BuildRemoteAdjacency(const std::vector<std::array<int, 4>> &conn,
                          const std::map<int, std::vector<int>> &neighbor_nodes,
                          idx_t base_gid,
                          std::vector<std::unordered_set<idx_t>> &adj,
                          MPI_Comm comm) {
  for (const auto &p : neighbor_nodes) {
    int neighbor_rank = p.first;
    const std::vector<int> &nodes = p.second;

    std::unordered_map<int, int> node_to_pos;
    node_to_pos.reserve(nodes.size() * 2);
    for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
      node_to_pos[nodes[i]] = i;
    }

    std::map<FaceKey, int> face_to_elem;
    for (int e = 0; e < static_cast<int>(conn.size()); ++e) {
      const auto &c = conn[e];
      const int faces[4][3] = {
          {c[1], c[2], c[3]},
          {c[0], c[2], c[3]},
          {c[0], c[1], c[3]},
          {c[0], c[1], c[2]},
      };
      for (int f = 0; f < 4; ++f) {
        int pos[3];
        bool all_shared = true;
        for (int i = 0; i < 3; ++i) {
          auto it = node_to_pos.find(faces[f][i]);
          if (it == node_to_pos.end()) {
            all_shared = false;
            break;
          }
          pos[i] = it->second;
        }
        if (!all_shared) {
          continue;
        }
        std::array<int, 3> key = {pos[0], pos[1], pos[2]};
        std::sort(key.begin(), key.end());
        FaceKey fk{key};
        if (!face_to_elem.count(fk)) {
          face_to_elem[fk] = e;
        }
      }
    }

    int send_faces = static_cast<int>(face_to_elem.size());
    int recv_faces = 0;
    MPI_Sendrecv(&send_faces, 1, MPI_INT, neighbor_rank, 100,
                 &recv_faces, 1, MPI_INT, neighbor_rank, 100, comm, MPI_STATUS_IGNORE);

    std::vector<idx_t> sendbuf;
    sendbuf.reserve(static_cast<size_t>(send_faces) * 4);
    for (const auto &kv : face_to_elem) {
      sendbuf.push_back(static_cast<idx_t>(kv.first.idx[0]));
      sendbuf.push_back(static_cast<idx_t>(kv.first.idx[1]));
      sendbuf.push_back(static_cast<idx_t>(kv.first.idx[2]));
      sendbuf.push_back(base_gid + kv.second);
    }

    std::vector<idx_t> recvbuf(static_cast<size_t>(recv_faces) * 4);
    MPI_Sendrecv(sendbuf.data(), static_cast<int>(sendbuf.size()), IDX_T, neighbor_rank, 101,
                 recvbuf.data(), static_cast<int>(recvbuf.size()), IDX_T, neighbor_rank, 101,
                 comm, MPI_STATUS_IGNORE);

    for (int i = 0; i < recv_faces; ++i) {
      std::array<int, 3> key = {static_cast<int>(recvbuf[4 * i + 0]),
                                static_cast<int>(recvbuf[4 * i + 1]),
                                static_cast<int>(recvbuf[4 * i + 2])};
      std::sort(key.begin(), key.end());
      FaceKey fk{key};
      auto it = face_to_elem.find(fk);
      if (it == face_to_elem.end()) {
        continue;
      }
      idx_t local_elem = it->second;
      idx_t remote_gid = recvbuf[4 * i + 3];
      AddAdj(adj, local_elem, remote_gid);
    }
  }
}

void BuildParMETISGraph(const std::vector<std::unordered_set<idx_t>> &adj,
                        std::vector<idx_t> &xadj,
                        std::vector<idx_t> &adjncy) {
  idx_t nvtxs = static_cast<idx_t>(adj.size());
  xadj.resize(static_cast<size_t>(nvtxs) + 1);
  adjncy.clear();
  adjncy.reserve(adj.size() * 8);

  xadj[0] = 0;
  for (idx_t i = 0; i < nvtxs; ++i) {
    for (idx_t neigh : adj[i]) {
      adjncy.push_back(neigh);
    }
    xadj[i + 1] = static_cast<idx_t>(adjncy.size());
  }
}

void BuildParMETISMesh(const std::vector<std::array<int, 4>> &conn,
                       const std::vector<idx_t> &global_node_id,
                       std::vector<idx_t> &eptr,
                       std::vector<idx_t> &eind) {
  const idx_t nelem = static_cast<idx_t>(conn.size());
  eptr.resize(static_cast<size_t>(nelem) + 1);
  eind.resize(static_cast<size_t>(nelem) * 4);
  eptr[0] = 0;
  for (idx_t e = 0; e < nelem; ++e) {
    eptr[e + 1] = eptr[e] + 4;
    for (int i = 0; i < 4; ++i) {
      eind[4 * e + i] = global_node_id[conn[static_cast<size_t>(e)][i]];
    }
  }
}

bool RunAdaptiveRepart(const std::vector<idx_t> &vtxdist,
                       const std::vector<idx_t> &xadj,
                       const std::vector<idx_t> &adjncy,
                       real_t ipc2redist,
                       std::vector<idx_t> &part,
                       MPI_Comm comm) {
  int nprocs = 1;
  MPI_Comm_size(comm, &nprocs);

  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon = 1;
  idx_t nparts = static_cast<idx_t>(nprocs);
  std::vector<real_t> tpwgts(static_cast<size_t>(nparts) * ncon, 1.0 / static_cast<real_t>(nparts));
  std::vector<real_t> ubvec(static_cast<size_t>(ncon), 1.05);
  std::vector<idx_t> options(3, 0);
  idx_t edgecut = 0;

  int status = ParMETIS_V3_AdaptiveRepart(
      const_cast<idx_t *>(vtxdist.data()), const_cast<idx_t *>(xadj.data()),
      const_cast<idx_t *>(adjncy.data()),
      NULL, NULL, NULL,
      &wgtflag, &numflag, &ncon, &nparts,
      tpwgts.data(), ubvec.data(), &ipc2redist,
      options.data(), &edgecut, part.data(), &comm);

  return status == METIS_OK;
}

bool RunPartMeshKway(const std::vector<idx_t> &elmdist,
                     const std::vector<idx_t> &eptr,
                     const std::vector<idx_t> &eind,
                     idx_t ncommonnodes,
                     std::vector<idx_t> &part,
                     MPI_Comm comm) {
  int nprocs = 1;
  MPI_Comm_size(comm, &nprocs);

  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon = 1;
  idx_t nparts = static_cast<idx_t>(nprocs);
  std::vector<real_t> tpwgts(static_cast<size_t>(nparts) * ncon, 1.0 / static_cast<real_t>(nparts));
  std::vector<real_t> ubvec(static_cast<size_t>(ncon), 1.05);
  std::vector<idx_t> options(3, 0);
  idx_t edgecut = 0;

  int status = ParMETIS_V3_PartMeshKway(
      const_cast<idx_t *>(elmdist.data()), const_cast<idx_t *>(eptr.data()),
      const_cast<idx_t *>(eind.data()), NULL,
      &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts,
      tpwgts.data(), ubvec.data(), options.data(), &edgecut,
      part.data(), &comm);

  return status == METIS_OK;
}

void WritePartition(int myid, const std::vector<idx_t> &part) {
  long fid = 0;
  const std::string filename = RankFilename("rdata", myid, ".repart.h5");
  hdf_create_file_(filename.c_str(), &fid);
  long size = static_cast<long>(part.size());
  std::vector<int> part_i(part.size());
  for (size_t i = 0; i < part.size(); ++i) {
    part_i[i] = static_cast<int>(part[i]);
  }
  hdf_write_int_array_(&fid, "/part", part_i.data(), &size);
  hdf_close_file_(&fid);
}

void WriteRepartMesh(int myid,
                     int nmaterial,
                     const std::vector<std::array<int, 4>> &conn,
                     const std::vector<int> &matid,
                     const std::vector<std::array<double, 3>> &coor,
                     const std::map<int, std::vector<int>> &neighbor_nodes,
                     const std::vector<idx_t> &global_node_id) {
  long fid = 0;
  const std::string filename = RankFilename("rdata", myid, ".repart.h5");
  hdf_create_file_(filename.c_str(), &fid);

  const int ne_tet4 = static_cast<int>(conn.size());
  const int ne_vox8 = 0;

  std::vector<int> neighborrankid;
  neighborrankid.reserve(neighbor_nodes.size());
  std::vector<int> mpinode_pointer;
  std::vector<int> mpinode_index;
  mpinode_pointer.push_back(0);
  for (const auto &p : neighbor_nodes) {
    neighborrankid.push_back(p.first);
    for (int node_id : p.second) {
      mpinode_index.push_back(node_id);
    }
    mpinode_pointer.push_back(static_cast<int>(mpinode_index.size()));
  }

  const int nnode = static_cast<int>(coor.size());
  std::vector<int> settings(7, 0);
  settings[0] = nnode;
  settings[1] = ne_vox8;
  settings[2] = ne_tet4;
  settings[3] = 0;  // surf nodes
  settings[4] = static_cast<int>(neighborrankid.size());
  settings[5] = static_cast<int>(mpinode_pointer.size());
  settings[6] = static_cast<int>(mpinode_index.size());
  long size = static_cast<long>(settings.size());
  hdf_write_int_array_(&fid, "/settings", settings.data(), &size);

  if (!coor.empty()) {
    std::vector<double> coor_buf(coor.size() * 3);
    for (size_t i = 0; i < coor.size(); ++i) {
      coor_buf[3 * i + 0] = coor[i][0];
      coor_buf[3 * i + 1] = coor[i][1];
      coor_buf[3 * i + 2] = coor[i][2];
    }
    size = static_cast<long>(coor_buf.size());
    hdf_write_double_array_(&fid, "/coor", coor_buf.data(), &size);
  }

  if (!conn.empty()) {
    std::vector<int> conn_tet4(conn.size() * 5);
    for (size_t e = 0; e < conn.size(); ++e) {
      for (int i = 0; i < 4; ++i) {
        conn_tet4[5 * e + i] = conn[e][i];
      }
      conn_tet4[5 * e + 4] = matid[e];
    }
    size = static_cast<long>(conn_tet4.size());
    hdf_write_int_array_(&fid, "/conn_tet4", conn_tet4.data(), &size);
  }

  // empty vox8 and surf nodes
  std::vector<int> empty_int;
  size = 0;
  hdf_write_int_array_(&fid, "/conn_vox8", empty_int.data(), &size);
  hdf_write_int_array_(&fid, "/surfnodeid", empty_int.data(), &size);

  if (!neighborrankid.empty()) {
    size = static_cast<long>(neighborrankid.size());
    hdf_write_int_array_(&fid, "/neighborrankid", neighborrankid.data(), &size);
  } else {
    hdf_write_int_array_(&fid, "/neighborrankid", empty_int.data(), &size);
  }

  if (!mpinode_pointer.empty()) {
    size = static_cast<long>(mpinode_pointer.size());
    hdf_write_int_array_(&fid, "/mpinode_pointer", mpinode_pointer.data(), &size);
  } else {
    hdf_write_int_array_(&fid, "/mpinode_pointer", empty_int.data(), &size);
  }

  if (!mpinode_index.empty()) {
    size = static_cast<long>(mpinode_index.size());
    hdf_write_int_array_(&fid, "/mpinode_index", mpinode_index.data(), &size);
  } else {
    hdf_write_int_array_(&fid, "/mpinode_index", empty_int.data(), &size);
  }

  if (!global_node_id.empty()) {
    std::vector<int> global_id_i(global_node_id.size());
    for (size_t i = 0; i < global_node_id.size(); ++i) {
      global_id_i[i] = static_cast<int>(global_node_id[i]) + 1;  // 1-based
    }
    size = static_cast<long>(global_id_i.size());
    hdf_write_int_array_(&fid, "/global_node_id", global_id_i.data(), &size);
  }

  hdf_close_file_(&fid);
}

std::vector<int> ComputeNodeOwners(const std::map<int, std::vector<int>> &neighbor_nodes,
                                   int nnode, int myid) {
  std::vector<int> owner(nnode, myid);
  for (const auto &p : neighbor_nodes) {
    int neighbor_rank = p.first;
    for (int node_id : p.second) {
      owner[node_id] = std::min(owner[node_id], neighbor_rank);
    }
  }
  return owner;
}

void AssignGlobalNodeIds(const std::vector<int> &owner,
                         std::vector<idx_t> &global_id,
                         std::vector<idx_t> &node_dist,
                         MPI_Comm comm) {
  int myid = 0;
  int nprocs = 1;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &nprocs);

  int nnode = static_cast<int>(owner.size());
  std::vector<int> owned_nodes;
  owned_nodes.reserve(nnode);
  for (int i = 0; i < nnode; ++i) {
    if (owner[i] == myid) {
      owned_nodes.push_back(i);
    }
  }

  idx_t local_owned = static_cast<idx_t>(owned_nodes.size());
  std::vector<idx_t> owned_counts(nprocs, 0);
  MPI_Allgather(&local_owned, 1, IDX_T, owned_counts.data(), 1, IDX_T, comm);

  node_dist.assign(nprocs + 1, 0);
  for (int p = 0; p < nprocs; ++p) {
    node_dist[p + 1] = node_dist[p] + owned_counts[p];
  }

  global_id.assign(nnode, -1);
  for (idx_t i = 0; i < local_owned; ++i) {
    int node = owned_nodes[i];
    global_id[node] = node_dist[myid] + i;
  }
}

void ExchangeOwnedGlobalIds(const std::map<int, std::vector<int>> &neighbor_nodes,
                            const std::vector<int> &owner,
                            std::vector<idx_t> &global_id,
                            MPI_Comm comm) {
  int myid = 0;
  MPI_Comm_rank(comm, &myid);
  for (const auto &p : neighbor_nodes) {
    int neighbor_rank = p.first;
    const std::vector<int> &nodes = p.second;

    std::vector<idx_t> sendbuf(nodes.size(), -1);
    for (size_t i = 0; i < nodes.size(); ++i) {
      int node = nodes[i];
      if (owner[node] == myid) {
        sendbuf[i] = global_id[node];
      }
    }

    std::vector<idx_t> recvbuf(nodes.size(), -1);
    MPI_Sendrecv(sendbuf.data(), static_cast<int>(sendbuf.size()), IDX_T, neighbor_rank, 200,
                 recvbuf.data(), static_cast<int>(recvbuf.size()), IDX_T, neighbor_rank, 200,
                 comm, MPI_STATUS_IGNORE);

    for (size_t i = 0; i < nodes.size(); ++i) {
      int node = nodes[i];
      if (owner[node] != myid && recvbuf[i] != -1) {
        global_id[node] = recvbuf[i];
      }
    }
  }
}

int OwnerOfGid(const std::vector<idx_t> &node_dist, idx_t gid) {
  auto it = std::upper_bound(node_dist.begin(), node_dist.end(), gid);
  int owner = static_cast<int>(std::distance(node_dist.begin(), it)) - 1;
  return owner;
}

void ExchangeNodeCoordinates(const std::vector<idx_t> &node_dist,
                             const std::vector<idx_t> &local_gids,
                             const std::unordered_map<idx_t, std::array<double, 3>> &owned_coord,
                             std::unordered_map<idx_t, std::array<double, 3>> &coord_map,
                             MPI_Comm comm) {
  int myid = 0;
  int nprocs = 1;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &nprocs);

  std::vector<std::vector<idx_t>> req_to_rank(nprocs);
  for (idx_t gid : local_gids) {
    if (coord_map.find(gid) != coord_map.end()) {
      continue;
    }
    int owner = OwnerOfGid(node_dist, gid);
    if (owner == myid) {
      continue;
    }
    req_to_rank[owner].push_back(gid);
  }

  std::vector<int> send_counts(nprocs, 0), recv_counts(nprocs, 0);
  for (int r = 0; r < nprocs; ++r) {
    send_counts[r] = static_cast<int>(req_to_rank[r].size());
  }
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);

  std::vector<int> send_displs(nprocs, 0), recv_displs(nprocs, 0);
  int send_total = 0;
  int recv_total = 0;
  for (int r = 0; r < nprocs; ++r) {
    send_displs[r] = send_total;
    recv_displs[r] = recv_total;
    send_total += send_counts[r];
    recv_total += recv_counts[r];
  }

  std::vector<idx_t> sendbuf(send_total);
  for (int r = 0; r < nprocs; ++r) {
    std::copy(req_to_rank[r].begin(), req_to_rank[r].end(), sendbuf.begin() + send_displs[r]);
  }

  std::vector<idx_t> recvbuf(recv_total);
  MPI_Alltoallv(sendbuf.data(), send_counts.data(), send_displs.data(), IDX_T,
                recvbuf.data(), recv_counts.data(), recv_displs.data(), IDX_T, comm);

  std::vector<double> send_coord(recv_total * 3, 0.0);
  for (int i = 0; i < recv_total; ++i) {
    idx_t gid = recvbuf[i];
    auto it = owned_coord.find(gid);
    if (it != owned_coord.end()) {
      send_coord[3 * i + 0] = it->second[0];
      send_coord[3 * i + 1] = it->second[1];
      send_coord[3 * i + 2] = it->second[2];
    }
  }

  std::vector<int> send_counts3(nprocs, 0), recv_counts3(nprocs, 0);
  std::vector<int> send_displs3(nprocs, 0), recv_displs3(nprocs, 0);
  for (int r = 0; r < nprocs; ++r) {
    send_counts3[r] = recv_counts[r] * 3;
    recv_counts3[r] = send_counts[r] * 3;
    send_displs3[r] = recv_displs[r] * 3;
    recv_displs3[r] = send_displs[r] * 3;
  }

  std::vector<double> recv_coord(send_total * 3, 0.0);
  MPI_Alltoallv(send_coord.data(), send_counts3.data(), send_displs3.data(), MPI_DOUBLE,
                recv_coord.data(), recv_counts3.data(), recv_displs3.data(), MPI_DOUBLE, comm);

  for (int r = 0; r < nprocs; ++r) {
    int count = send_counts[r];
    int disp = send_displs[r];
    for (int i = 0; i < count; ++i) {
      idx_t gid = sendbuf[disp + i];
      std::array<double, 3> c = {recv_coord[3 * (disp + i) + 0],
                                 recv_coord[3 * (disp + i) + 1],
                                 recv_coord[3 * (disp + i) + 2]};
      coord_map[gid] = c;
    }
  }

}

void BuildNeighborNodesFromGlobalIds(const std::vector<idx_t> &local_gids,
                                     const std::unordered_map<idx_t, int> &gid_to_local,
                                     std::map<int, std::vector<int>> &neighbor_nodes,
                                     MPI_Comm comm) {
  int myid = 0;
  int nprocs = 1;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &nprocs);

  std::vector<int> counts(nprocs, 0);
  int local_count = static_cast<int>(local_gids.size());
  MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

  std::vector<int> displs(nprocs, 0);
  int total = 0;
  for (int r = 0; r < nprocs; ++r) {
    displs[r] = total;
    total += counts[r];
  }

  std::vector<idx_t> all_gids(total);
  MPI_Allgatherv(local_gids.data(), local_count, IDX_T,
                 all_gids.data(), counts.data(), displs.data(), IDX_T, comm);

  std::unordered_set<idx_t> local_set(local_gids.begin(), local_gids.end());
  for (int r = 0; r < nprocs; ++r) {
    if (r == myid) {
      continue;
    }
    int offset = displs[r];
    int count = counts[r];
    std::vector<idx_t> shared;
    shared.reserve(count);
    for (int i = 0; i < count; ++i) {
      idx_t gid = all_gids[offset + i];
      if (local_set.find(gid) != local_set.end()) {
        shared.push_back(gid);
      }
    }
    if (shared.empty()) {
      continue;
    }
    std::sort(shared.begin(), shared.end());
    std::vector<int> nodes;
    nodes.reserve(shared.size());
    for (idx_t gid : shared) {
      nodes.push_back(gid_to_local.at(gid));
    }
    neighbor_nodes[r] = std::move(nodes);
  }
}

void print_n(int nnode, int nelem, int myid, int numprocs) {
  std::vector<int> all_nnode(numprocs, 0);
  MPI_Gather(&nnode, 1, MPI_INT, all_nnode.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<int> all_nelem(numprocs, 0);
  MPI_Gather(&nelem, 1, MPI_INT, all_nelem.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid == 0) {
    std::cout << "nnode: " << std::endl;
    int nmin, nmax;
    long nave, nstd;
    nmin = *std::min_element(all_nnode.begin(), all_nnode.end());
    nmax = *std::max_element(all_nnode.begin(), all_nnode.end());
    nave = std::accumulate(all_nnode.begin(), all_nnode.end(), 0) / static_cast<double>(numprocs);
    nstd = 0;
    for (int n : all_nnode) {
      nstd += (n - nave) * (n - nave);
    }
    nstd = std::sqrt(nstd / numprocs);
    std::cout << "  min: " << nmin << ", max: " << nmax << ", ave: " << nave << ", std: " << nstd << std::endl;

    std::cout << "nelem: " << std::endl;
    nmin = *std::min_element(all_nelem.begin(), all_nelem.end());
    nmax = *std::max_element(all_nelem.begin(), all_nelem.end());
    nave = std::accumulate(all_nelem.begin(), all_nelem.end(), 0) / static_cast<double>(numprocs);
    nstd = 0;
    for (int n : all_nelem) {
      nstd += (n - nave) * (n - nave);
    }
    nstd = std::sqrt(nstd / numprocs);
    std::cout << "  min: " << nmin << ", max: " << nmax << ", ave: " << nave << ", std: " << nstd << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

PartitionMethod ParsePartitionMethod(int argc, char **argv, real_t &ipc2redist) {
  ipc2redist = 1.0;
  PartitionMethod method = PartitionMethod::AdaptiveRepart;

  if (argc >= 2) {
    std::string arg1(argv[1]);
    if (arg1 == "kway" || arg1 == "meshkway" || arg1 == "partmeshkway") {
      method = PartitionMethod::PartMeshKway;
    } else if (arg1 == "adaptive") {
      method = PartitionMethod::AdaptiveRepart;
    } else {
      ipc2redist = static_cast<real_t>(std::atof(argv[1]));
    }
  }

  if (argc >= 3) {
    std::string arg2(argv[2]);
    if (arg2 == "kway" || arg2 == "meshkway" || arg2 == "partmeshkway") {
      method = PartitionMethod::PartMeshKway;
    } else if (arg2 == "adaptive") {
      method = PartitionMethod::AdaptiveRepart;
    }
  }

  return method;
}

}  // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int myid = 0;
  int nprocs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  real_t ipc2redist = 1.0;
  PartitionMethod method = ParsePartitionMethod(argc, argv, ipc2redist);

  int nnode = 0;
  int nelem = 0;
  int nmaterial = 0;
  ReadSetting(myid, nnode, nelem, nmaterial);

  std::vector<std::array<int, 4>> conn;
  std::vector<int> matid;
  ReadConn(myid, nelem, conn, matid);

  std::vector<std::array<double, 3>> coor;
  ReadCoor(myid, nnode, coor);

  print_n(nnode, nelem, myid, nprocs);

  std::map<int, std::vector<int>> neighbor_nodes;
  ReadMPInode(myid, neighbor_nodes);

  // global node id construction (needed for mesh-based partitioning and redistribution)
  std::vector<int> owner = ComputeNodeOwners(neighbor_nodes, nnode, myid);
  std::vector<idx_t> global_id;
  std::vector<idx_t> node_dist;
  AssignGlobalNodeIds(owner, global_id, node_dist, MPI_COMM_WORLD);
  ExchangeOwnedGlobalIds(neighbor_nodes, owner, global_id, MPI_COMM_WORLD);

  for (int i = 0; i < nnode; ++i) {
    if (global_id[i] < 0) {
      std::cerr << "ERROR: global node id unresolved for node " << i
                << " on rank " << myid << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }

  idx_t local_nelem = static_cast<idx_t>(nelem);
  std::vector<idx_t> vtxdist(static_cast<size_t>(nprocs) + 1, 0);
  std::vector<idx_t> nelems_all(static_cast<size_t>(nprocs), 0);
  MPI_Allgather(&local_nelem, 1, IDX_T, nelems_all.data(), 1, IDX_T, MPI_COMM_WORLD);
  for (int p = 0; p < nprocs; ++p) {
    vtxdist[p + 1] = vtxdist[p] + nelems_all[p];
  }
  idx_t base_gid = vtxdist[myid];

  std::vector<idx_t> part(static_cast<size_t>(nelem), 0);
  MPI_Comm comm = MPI_COMM_WORLD;
  bool ok = true;
  if (method == PartitionMethod::AdaptiveRepart) {
    std::vector<std::unordered_set<idx_t>> adj(static_cast<size_t>(nelem));
    BuildLocalAdjacency(conn, base_gid, adj);
    BuildRemoteAdjacency(conn, neighbor_nodes, base_gid, adj, MPI_COMM_WORLD);

    std::vector<idx_t> xadj;
    std::vector<idx_t> adjncy;
    BuildParMETISGraph(adj, xadj, adjncy);

    ok = RunAdaptiveRepart(vtxdist, xadj, adjncy, ipc2redist, part, comm);
  } else {
    const idx_t ncommonnodes = 3;  // tet4: face sharing
    std::vector<idx_t> eptr;
    std::vector<idx_t> eind;
    BuildParMETISMesh(conn, global_id, eptr, eind);
    ok = RunPartMeshKway(vtxdist, eptr, eind, ncommonnodes, part, comm);
  }

  if (!ok) {
    if (myid == 0) {
      if (method == PartitionMethod::AdaptiveRepart) {
        std::cerr << "ParMETIS_V3_AdaptiveRepart failed." << std::endl;
      } else {
        std::cerr << "ParMETIS_V3_PartMeshKway failed." << std::endl;
      }
    }
    MPI_Finalize();
    return 1;
  }

  // WritePartition(myid, part);

  // redistribute elements based on part
  MPI_Barrier(MPI_COMM_WORLD);
  std::vector<std::vector<int>> send_elems(nprocs);
  std::vector<std::vector<int>> send_mat(nprocs);
  for (int e = 0; e < nelem; ++e) {
    int dest = static_cast<int>(part[e]);
    for (int i = 0; i < 4; ++i) {
      send_elems[dest].push_back(static_cast<int>(global_id[conn[e][i]]));
    }
    send_mat[dest].push_back(matid[e]);
  }

  std::vector<int> send_counts(nprocs, 0), recv_counts(nprocs, 0);
  for (int r = 0; r < nprocs; ++r) {
    send_counts[r] = static_cast<int>(send_elems[r].size());
    // if (r != myid) {
    //   std::cout << "rank " << myid << " -> " << r << ": " << send_counts[r] / 4 << " elements" << std::endl;
    // }
  }
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> send_displs(nprocs, 0), recv_displs(nprocs, 0);
  int send_total = 0;
  int recv_total = 0;
  for (int r = 0; r < nprocs; ++r) {
    send_displs[r] = send_total;
    recv_displs[r] = recv_total;
    send_total += send_counts[r];
    recv_total += recv_counts[r];
  }

  std::vector<int> sendbuf(send_total);
  for (int r = 0; r < nprocs; ++r) {
    std::copy(send_elems[r].begin(), send_elems[r].end(), sendbuf.begin() + send_displs[r]);
  }
  std::vector<int> recvbuf(recv_total);
  MPI_Alltoallv(sendbuf.data(), send_counts.data(), send_displs.data(), MPI_INT,
                recvbuf.data(), recv_counts.data(), recv_displs.data(), MPI_INT, MPI_COMM_WORLD);

  std::vector<int> send_mat_counts(nprocs, 0), recv_mat_counts(nprocs, 0);
  for (int r = 0; r < nprocs; ++r) {
    send_mat_counts[r] = static_cast<int>(send_mat[r].size());
  }
  MPI_Alltoall(send_mat_counts.data(), 1, MPI_INT,
               recv_mat_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  std::vector<int> send_mat_displs(nprocs, 0), recv_mat_displs(nprocs, 0);
  int send_mat_total = 0;
  int recv_mat_total = 0;
  for (int r = 0; r < nprocs; ++r) {
    send_mat_displs[r] = send_mat_total;
    recv_mat_displs[r] = recv_mat_total;
    send_mat_total += send_mat_counts[r];
    recv_mat_total += recv_mat_counts[r];
  }
  std::vector<int> sendbuf_mat(send_mat_total);
  for (int r = 0; r < nprocs; ++r) {
    std::copy(send_mat[r].begin(), send_mat[r].end(), sendbuf_mat.begin() + send_mat_displs[r]);
  }
  std::vector<int> recvbuf_mat(recv_mat_total);
  MPI_Alltoallv(sendbuf_mat.data(), send_mat_counts.data(), send_mat_displs.data(), MPI_INT,
                recvbuf_mat.data(), recv_mat_counts.data(), recv_mat_displs.data(), MPI_INT, MPI_COMM_WORLD);

  int new_nelem = recv_mat_total;
  std::vector<std::array<idx_t, 4>> new_elem_gids(new_nelem);
  for (int e = 0; e < new_nelem; ++e) {
    for (int i = 0; i < 4; ++i) {
      new_elem_gids[e][i] = static_cast<idx_t>(recvbuf[4 * e + i]);
    }
  }

  std::unordered_set<idx_t> gid_set;
  gid_set.reserve(new_nelem * 4);
  for (const auto &elem : new_elem_gids) {
    for (int i = 0; i < 4; ++i) {
      gid_set.insert(elem[i]);
    }
  }

  std::vector<idx_t> local_gids(gid_set.begin(), gid_set.end());
  std::sort(local_gids.begin(), local_gids.end());
  std::unordered_map<idx_t, int> gid_to_local;
  gid_to_local.reserve(local_gids.size() * 2);
  for (size_t i = 0; i < local_gids.size(); ++i) {
    gid_to_local[local_gids[i]] = static_cast<int>(i);
  }

  std::unordered_map<idx_t, std::array<double, 3>> coord_map;
  coord_map.reserve(local_gids.size() * 2);
  std::unordered_map<idx_t, std::array<double, 3>> owned_coord;
  owned_coord.reserve(nnode);
  for (int i = 0; i < nnode; ++i) {
    if (owner[i] == myid) {
      owned_coord[global_id[i]] = coor[i];
      if (gid_to_local.find(global_id[i]) != gid_to_local.end()) {
        coord_map[global_id[i]] = coor[i];
      }
    }
  }
  ExchangeNodeCoordinates(node_dist, local_gids, owned_coord, coord_map, MPI_COMM_WORLD);

  std::vector<std::array<double, 3>> new_coor(local_gids.size());
  for (size_t i = 0; i < local_gids.size(); ++i) {
    auto it = coord_map.find(local_gids[i]);
    if (it == coord_map.end()) {
      std::cerr << "ERROR: missing coordinate for gid " << local_gids[i]
                << " on rank " << myid << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 3);
    }
    new_coor[i] = it->second;
  }

  std::vector<std::array<int, 4>> new_conn(new_nelem);
  for (int e = 0; e < new_nelem; ++e) {
    for (int i = 0; i < 4; ++i) {
      new_conn[e][i] = gid_to_local[new_elem_gids[e][i]];
    }
  }

  print_n(new_coor.size(), new_conn.size(), myid, nprocs);

  std::vector<int> new_matid = recvbuf_mat;

  std::map<int, std::vector<int>> new_neighbor_nodes;
  BuildNeighborNodesFromGlobalIds(local_gids, gid_to_local, new_neighbor_nodes, MPI_COMM_WORLD);

  WriteRepartMesh(myid, nmaterial, new_conn, new_matid, new_coor, new_neighbor_nodes, local_gids);

  MPI_Finalize();
  return 0;
}
