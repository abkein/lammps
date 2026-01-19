/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// TODO: NUCC FILE

#include "compute_cluster_volume.h"
#include "compute_cluster_size_ext.h"
#include "nucc_cspan.hpp"
#include <cmath>
#include <cstddef>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <unordered_map>

using namespace LAMMPS_NS;
using NUCC::cspan;
using NUCC::cluster_data;

/* ---------------------------------------------------------------------- */

ComputeClusterVolume::ComputeClusterVolume(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), nloc_grid(0), n_cells(0), precompute(false), noff(0), nloc_recv(0)
{
  vector_flag = 1;
  size_vector = 0;
  extvector = 0;
  local_flag = 1;
  size_local_rows = 0;
  size_local_cols = 0;

  if (narg < 8) { utils::missing_cmd_args(FLERR, "compute cluster/volume", error); }

  // Parse arguments //

  // Get cluster/size compute
  compute_cluster_size =
      dynamic_cast<ComputeClusterSizeExt *>(lmp->modify->get_compute_by_id(arg[3]));
  if (compute_cluster_size == nullptr) {
    error->all(FLERR, "compute {}: Cannot find compute with style 'cluster/size' with id: {}",
               style, arg[3]);
  }

  // Get the critical size
  size_cutoff = compute_cluster_size->get_size_cutoff();
  if (::strcmp(arg[4], "inherit") != 0) {
    const int t_size_cutoff = utils::inumeric(FLERR, arg[4], true, lmp);
    if (t_size_cutoff < 1) {
      error->all(FLERR, "size_cutoff for compute {} must be greater than 0", style);
    }
    if (t_size_cutoff > size_cutoff) {
      error->all(FLERR,
                 "size_cutoff for compute {} cannot be greater than it of compute cluster/size",
                 style);
    }
  }

  if (::strcmp(arg[5], "rect") == 0) {
    mode = VOLUMEMODE::RECTANGLE;
    // } else if (::strcmp(arg[5], "sphere") == 0) {
    //   mode = VOLUMEMODE::SPHERE;
  } else if (::strcmp(arg[5], "calc") == 0) {
    mode = VOLUMEMODE::CALC;
  } else {
    error->all(FLERR, "Unknown mode for compute {}: {}", style, arg[5]);
  }

  overlap = utils::numeric(FLERR, arg[6], true, lmp);
  overlap_sq = overlap * overlap;
  if (overlap < 0) { error->all(FLERR, "Minimum distance for {} must be non-negative", style); }

  voxel_size = utils::numeric(FLERR, arg[7], true, lmp);
  if (voxel_size < 0) { error->all(FLERR, "voxel size for {} must be non-negative", style); }
  if (voxel_size > overlap) {
    error->all(FLERR, "voxel size for {} must be less than overlap", style);
  }

  int iarg = 8;

  while (iarg < narg) {
    if (::strcmp(arg[iarg], "precompute") == 0) {
      precompute = true;
      iarg += 1;
    } else {
      error->all(FLERR, "Illegal {} option {}", style, arg[iarg]);
    }
  }

  n_cells = static_cast<int>(overlap / voxel_size);

  size_local_rows = size_cutoff + 1;
  size_vector = size_cutoff + 1;

  dist_local.create(memory, size_local_rows, "compute:cluster/ke:local_kes");
  vector_local = dist_local.data();
  dist.create(memory, size_vector, "compute:cluster/ke:kes");
  vector = dist.data();
  dist_global.create(memory, size_vector * comm->nprocs, "dist_global");
}

/* ---------------------------------------------------------------------- */

ComputeClusterVolume::~ComputeClusterVolume() noexcept(true)
{
  volumes.destroy(memory);
  dist.destroy(memory);
  dist_local.destroy(memory);
  occupancy_grid.destroy(memory);
  offsets.destroy(memory);
  bboxes.destroy(memory);
  recv_buf.destroy(memory);
  dist_global.destroy(memory);
}

/* ---------------------------------------------------------------------- */

void ComputeClusterVolume::init()
{
  // ::MPI_Comm_set_errhandler(world, MPI_ERRORS_RETURN);
  if ((modify->get_compute_by_style(style).size() > 1) && (comm->me == 0)) {
    error->warning(FLERR, "More than one compute {}", style);
  }

  subbonds[0 + 0] = domain->sublo[0];
  subbonds[0 + 1] = domain->subhi[0];
  subbonds[2 + 0] = domain->sublo[1];
  subbonds[2 + 1] = domain->subhi[1];
  subbonds[4 + 0] = domain->sublo[2];
  subbonds[4 + 1] = domain->subhi[2];

  if ((nloc < atom->nlocal) || volumes.empty()) {
    nloc = static_cast<int>(atom->nlocal * LMP_NUCC_ALLOC_COEFF);
    volumes.grow(memory, nloc, "ComputeClusterVolume:volumes");
    if (mode == VOLUMEMODE::RECTANGLE) {
      bboxes.create(memory, 6 * nloc, "ComputeClusterVolume:bboxes");
    }
  }


  // memory->create(recv_comm_matrix_local, comm->nprocs, "recv_comm_matrix_local");
  // memory->create(recv_comm_matrix_global, comm->nprocs * comm->nprocs, "recv_comm_matrix_global");
  // memory->create(send_comm_matrix_local, comm->nprocs, "send_comm_matrix_local");
  // memory->create(send_comm_matrix_global, comm->nprocs * comm->nprocs, "send_comm_matrix_global");

  if (precompute) {
    const auto radius_voxels_int = static_cast<int>(std::ceil(overlap / voxel_size));
    noff =
        3 * 2 * (radius_voxels_int + 1) * 2 * (radius_voxels_int + 1) * 2 * (radius_voxels_int + 1);
    if (offsets.empty()) { offsets.create(memory, noff, "ComputeClusterVolume:offsets"); }

    double const voxel_size_sq = voxel_size * voxel_size;
    noffsets = 0;
    for (int dx = -radius_voxels_int; dx <= radius_voxels_int; ++dx) {
      for (int dy = -radius_voxels_int; dy <= radius_voxels_int; ++dy) {
        for (int dz = -radius_voxels_int; dz <= radius_voxels_int; ++dz) {
          double const dist_sq = dx * dx + dy * dy + dz * dz;
          if (dist_sq * voxel_size_sq <= overlap_sq) {
            offsets[noffsets + 0] = dx;
            offsets[noffsets + 1] = dy;
            offsets[noffsets + 2] = dz;
            noffsets += 3;
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeClusterVolume::compute_vector()
{
  invoked_vector = update->ntimestep;

  compute_local();
  dist.reset();
  ::MPI_Allgather(dist_local.data(), size_vector, MPI_DOUBLE, dist_global.data(), size_vector,
                  MPI_DOUBLE, world);
  for (int i = 0; i < size_vector; ++i) {
    int n = 0;
    for (int j = 0; j < comm->nprocs; ++j) {
      const int k = j * comm->nprocs + i;
      if (dist_global[k] > 0) {
        dist[i] += dist_global[k];
        ++n;
      }
    }
    if (dist[i] > 0) { dist[i] /= n; }
  }
}

/* ---------------------------------------------------------------------- */

inline void maximize_bbox(const double *const src, double *const dest) noexcept
{
  if (src[0] < dest[0]) { dest[0] = src[0]; }
  if (src[1] < dest[1]) { dest[1] = src[1]; }
  if (src[2] < dest[2]) { dest[2] = src[2]; }
  if (src[3] > dest[3]) { dest[3] = src[3]; }
  if (src[4] > dest[4]) { dest[4] = src[4]; }
  if (src[5] > dest[5]) { dest[5] = src[5]; }
}

/* ---------------------------------------------------------------------- */

void ComputeClusterVolume::compute_local()
{
  invoked_local = update->ntimestep;

  if (compute_cluster_size->invoked_vector != update->ntimestep) {
    compute_cluster_size->compute_vector();
  }

  if (nloc < atom->nlocal) {
    nloc = static_cast<int>(atom->nlocal * LMP_NUCC_ALLOC_COEFF);
    volumes.grow(memory, nloc, "ComputeClusterVolume:volumes");

    if (mode == VOLUMEMODE::RECTANGLE) {
      bboxes.grow(memory, 6 * nloc, "ComputeClusterVolume:bboxes");
    }
  }

  if (nloc_recv < compute_cluster_size->get_nonexclusive()) {
    nloc_recv =
        static_cast<bigint>(compute_cluster_size->get_nonexclusive() * LMP_NUCC_ALLOC_COEFF);
    grow_ptr_array(in_reqs, comm->nprocs * nloc_recv, "ComputeClusterVolume:in_reqs");
    grow_ptr_array(out_reqs, comm->nprocs * nloc_recv, "ComputeClusterVolume:out_reqs");

    if (mode == VOLUMEMODE::CALC) {
      recv_buf.grow(memory, comm->nprocs * nloc_recv, "ComputeClusterVolume:recv_buf");
    } else if (mode == VOLUMEMODE::RECTANGLE) {
      recv_buf.grow(memory, static_cast<bigint>(6) * comm->nprocs * nloc_recv,
                    "ComputeClusterVolume:recv_buf");
    } else {
      // compliant
    }
  }

  // ::memset(recv_comm_matrix_local, 0, comm->nprocs * sizeof(int));
  // ::memset(recv_comm_matrix_global, 0, comm->nprocs * comm->nprocs * sizeof(int));
  // ::memset(send_comm_matrix_local, 0, comm->nprocs * sizeof(int));
  // ::memset(send_comm_matrix_global, 0, comm->nprocs * comm->nprocs * sizeof(int));

  int to_send = 0;
  int to_receive = 0;
  const auto &cmap = compute_cluster_size->get_cIDs_by_size();
  const auto &clusters = compute_cluster_size->get_clusters();
  if (mode == VOLUMEMODE::CALC) {
    for (const auto &[size, clidxs] : cmap) {
      for (const auto clidx : clidxs) {
        const cluster_data &clstr = clusters[clidx];
        bool const nonexclusive = clstr.nowners > 0;
        // volumes[clidx] = occupied_volume_precomputed(clstr.atoms_all(), clstr.l_size, clstr.nghost,
        //                                              nonexclusive);
        volumes[clidx] = occupied_volume_precomputed(clstr.atoms(), clstr.l_size, 0,
        nonexclusive);
        if (nonexclusive) {
          if (clstr.host != comm->me) {
            // ++send_comm_matrix_local[clstr.host];
            // if ((clstr.host < 0) || (clstr.host >= comm->nprocs)) { utils::logmesg(lmp, "Invalid rank: {}", clstr.host); }
            int err = ::MPI_Send_init(volumes.offset(clidx), 1, MPI_DOUBLE, clstr.host, clstr.clid,
                                      world, out_reqs + to_send);
            if (err > 0) {
              utils::logmesg(lmp, "MPI returned non-zero code: {}, rank: {}", err, clstr.host);
            }
            ++to_send;
          } else {
            for (int i = 0; i < clstr.nowners; ++i) {
              // if (clstr.owners[i] == comm->me) { utils::logmesg(lmp, "Me in owners"); }
              // ++recv_comm_matrix_local[clstr.owners[i]];
              ::MPI_Recv_init(recv_buf.offset(to_receive), 1, MPI_DOUBLE, clstr.owners()[i],
                              clstr.clid, world, in_reqs + to_receive);
              ++to_receive;
            }
          }
        }
      }
    }
  } else if (mode == VOLUMEMODE::RECTANGLE) {
    for (const auto &[size, clidxs] : cmap) {
      for (const auto clidx : clidxs) {
        const cluster_data &clstr = clusters[clidx];
        bigint const bbox_start = static_cast<bigint>(6) * clidx;
        cluster_bbox(clstr.atoms(), clstr.l_size /*+ clstr.nghost*/, bboxes.offset(bbox_start));
        if (clstr.nowners > 0) {
          if (clstr.host != comm->me) {
            ::MPI_Send_init(bboxes.offset(bbox_start), 6, MPI_DOUBLE, clstr.host, clstr.clid, world,
                            out_reqs + to_send);
            ++to_send;
          } else {
            for (int i = 0; i < clstr.nowners; ++i) {
              ::MPI_Recv_init(recv_buf.offset(6 * to_receive), 6, MPI_DOUBLE, clstr.owners()[i],
                              clstr.clid, world, in_reqs + to_receive);
              ++to_receive;
            }
          }
        }
      }
    }
  } else {
    // compliant
  }

  // ::MPI_Allgather(send_comm_matrix_local, comm->nprocs, MPI_INT, send_comm_matrix_global, comm->nprocs, MPI_INT, world);
  // ::MPI_Allgather(recv_comm_matrix_local, comm->nprocs, MPI_INT, recv_comm_matrix_global, comm->nprocs, MPI_INT, world);
  // if (comm->me == 0) {
  //   int _counts[2*comm->nprocs];
  //   ::memset(_counts, 0, 2 * comm->nprocs * sizeof(int));
  //   utils::logmesg(lmp, "Send matrix:\n");
  //   for (int i = 0; i < comm->nprocs; ++i) {
  //     for (int j = 0; j < comm->nprocs; ++j) {
  //       int s = send_comm_matrix_global[i*comm->nprocs + j];
  //       _counts[i] += s;
  //       _counts[comm->nprocs + j] += s;
  //       utils::logmesg(lmp, "{:3} ", s);
  //     }
  //     utils::logmesg(lmp, "|{:3}\n", _counts[i]);
  //   }
  //   for (int i = 0; i < comm->nprocs; ++i) {
  //     utils::logmesg(lmp, "----");
  //   }
  //   utils::logmesg(lmp, "\n");
  //   for (int i = 0; i < comm->nprocs; ++i) {
  //     utils::logmesg(lmp, "{:3} ", _counts[comm->nprocs + i]);
  //   }
  //   ::memset(_counts, 0, 2 * comm->nprocs * sizeof(int));
  //   utils::logmesg(lmp, "\nRecv matrix:\n");
  //   for (int i = 0; i < comm->nprocs; ++i) {
  //     for (int j = 0; j < comm->nprocs; ++j) {
  //       int s = recv_comm_matrix_global[i*comm->nprocs + j];
  //       _counts[i] += s;
  //       _counts[comm->nprocs + j] += s;
  //       utils::logmesg(lmp, "{:3} ", s);
  //     }
  //     utils::logmesg(lmp, "|{:3}\n", _counts[i]);
  //   }
  //   for (int i = 0; i < comm->nprocs; ++i) {
  //     utils::logmesg(lmp, "----");
  //   }
  //   utils::logmesg(lmp, "\n");
  //   for (int i = 0; i < comm->nprocs; ++i) {
  //     utils::logmesg(lmp, "{:3} ", _counts[comm->nprocs + i]);
  //   }
  //   utils::logmesg(lmp, "\n");
  // }

  if (to_send > 0) { ::MPI_Startall(to_send, out_reqs); }
  if (to_receive > 0) { ::MPI_Startall(to_receive, in_reqs); }
  if (to_send > 0) { ::MPI_Waitall(to_send, out_reqs, MPI_STATUSES_IGNORE /*out_stats*/); }
  if (to_receive > 0) { ::MPI_Waitall(to_receive, in_reqs, MPI_STATUSES_IGNORE /*in_stats*/); }

  int received = 0;
  if (mode == VOLUMEMODE::CALC) {
    for (const auto &[size, clidxs] : cmap) {
      for (const auto clidx : clidxs) {
        const cluster_data &clstr = clusters[clidx];
        if ((clstr.nowners > 0) && (clstr.host == comm->me)) {
          for (int i = 0; i < clstr.nowners; ++i) { volumes[clidx] += recv_buf[received++]; }
        }
      }
    }
  } else if (mode == VOLUMEMODE::RECTANGLE) {
    for (const auto &[size, clidxs] : cmap) {
      for (const auto clidx : clidxs) {
        const cluster_data &clstr = clusters[clidx];
        bigint const bbox_start = static_cast<bigint>(6) * clidx;
        if (clstr.host == comm->me) {
          if ((clstr.nowners > 0)) {
            for (int i = 0; i < clstr.nowners; ++i) {
              ::maximize_bbox(recv_buf.offset(6 * received), bboxes.offset(bbox_start));
              ++received;
            }
          }
          volumes[clidx] = (bboxes[bbox_start + 3] - bboxes[bbox_start + 0]) *
              (bboxes[bbox_start + 4] - bboxes[bbox_start + 1]) *
              (bboxes[bbox_start + 5] - bboxes[bbox_start + 2]);
        }
      }
    }
  } else {
    // compliant
  }

  dist_local.reset();
  for (const auto &[size, clidxs] : cmap) {
    int num_clst = 0;
    for (const auto clidx : clidxs) {
      const cluster_data &clstr = clusters[clidx];
      if (clstr.host == comm->me) {
        dist_local[size] += volumes[clidx];
        ++num_clst;
      }
    }
    if (num_clst > 0) { dist_local[size] /= num_clst; }
  }

  for (int i = 0; i < to_send; ++i) { ::MPI_Request_free(out_reqs + i); }
  for (int i = 0; i < to_receive; ++i) { ::MPI_Request_free(in_reqs + i); }
}

/* ---------------------------------------------------------------------- */

void ComputeClusterVolume::cluster_bbox(const cspan<const int> &list, const int n,
                                        double *bbox) const noexcept
{
  double **const centers = atom->x;
  // double bbox[6] =
  //   {
  //     centers[list[0]][0],
  //     centers[list[0]][1],
  //     centers[list[0]][2],
  //     centers[list[0]][0],
  //     centers[list[0]][1],
  //     centers[list[0]][2]
  //   };

  bbox[0] = centers[list[0]][0];
  bbox[1] = centers[list[0]][1];
  bbox[2] = centers[list[0]][2];
  bbox[3] = centers[list[0]][0];
  bbox[4] = centers[list[0]][1];
  bbox[5] = centers[list[0]][2];

  // find min and max center coordinates, i.e. bounding box
  for (int i = 1; i < n; ++i) {
    if (centers[list[i]][0] < bbox[0]) { bbox[0] = centers[list[i]][0]; }
    if (centers[list[i]][1] < bbox[1]) { bbox[1] = centers[list[i]][1]; }
    if (centers[list[i]][2] < bbox[2]) { bbox[2] = centers[list[i]][2]; }
    if (centers[list[i]][0] > bbox[3]) { bbox[3] = centers[list[i]][0]; }
    if (centers[list[i]][1] > bbox[4]) { bbox[4] = centers[list[i]][1]; }
    if (centers[list[i]][2] > bbox[5]) { bbox[5] = centers[list[i]][2]; }
  }

  bbox[0] -= overlap;
  bbox[1] -= overlap;
  bbox[2] -= overlap;
  bbox[3] += overlap;
  bbox[4] += overlap;
  bbox[5] += overlap;
}

/* ---------------------------------------------------------------------- */

[[nodiscard]] inline static bigint idx(const bigint i, const bigint j, const bigint k,
                                       const bigint ny, const bigint nz) noexcept
{
  return i * (nz * ny) + j * nz + k;
}

/* ---------------------------------------------------------------------- */

double ComputeClusterVolume::occupied_volume_precomputed(const cspan<const int> &list, const int n,
                                                         const int nghost,
                                                         const bool nonexclusive) noexcept
{
  double **const centers = atom->x;
  double bbox[6]{};
  cluster_bbox(list, n, bbox);

  if (nonexclusive) {
    bbox[0] = MAX(bbox[0], subbonds[0]);    // cut should be on only needed sides
    bbox[1] = MAX(bbox[1], subbonds[2]);
    bbox[2] = MAX(bbox[2], subbonds[4]);
    bbox[3] = MIN(bbox[3], subbonds[1]);
    bbox[4] = MIN(bbox[4], subbonds[3]);
    bbox[5] = MIN(bbox[5], subbonds[5]);
  }

  const bigint steps_x = static_cast<bigint>((bbox[3] - bbox[0]) / voxel_size) + 1;
  const bigint steps_y = static_cast<bigint>((bbox[4] - bbox[1]) / voxel_size) + 1;
  const bigint steps_z = static_cast<bigint>((bbox[5] - bbox[2]) / voxel_size) + 1;

  if (nloc_grid < steps_x * steps_y * steps_z) {
    nloc_grid = steps_x * steps_y * steps_z;
    occupancy_grid.grow(memory, nloc_grid, "compute_cluster_volume:grid");
  }
  occupancy_grid.reset();

  for (int l = 0; l < n + nghost; ++l) {
    const double *const center = centers[list[l]];
    const auto ni = static_cast<bigint>((center[0] - bbox[0]) / voxel_size);
    const auto nj = static_cast<bigint>((center[1] - bbox[1]) / voxel_size);
    const auto nk = static_cast<bigint>((center[2] - bbox[2]) / voxel_size);

    for (int noffset = 0; noffset < noffsets; noffset += 3) {
      const bigint i = ni + offsets[noffset + 0];
      const bigint j = nj + offsets[noffset + 1];
      const bigint k = nk + offsets[noffset + 2];

      if ((i >= 0) && (i < steps_x) && (j >= 0) && (j < steps_y) && (k >= 0) && (k < steps_z)) {
        occupancy_grid[::idx(i, j, k, steps_y, steps_z)] = true;
      }
    }
  }

  // count occupied cells
  bigint occupied = 0;
  for (bigint i = 0; i < steps_x; ++i) {
    for (bigint j = 0; j < steps_y; ++j) {
      for (bigint k = 0; k < steps_z; ++k) {
        if (occupancy_grid[::idx(i, j, k, steps_y, steps_z)]) { ++occupied; }
      }
    }
  }

  return occupied * voxel_size * voxel_size * voxel_size;
}

/* ---------------------------------------------------------------------- */

double ComputeClusterVolume::occupied_volume_grid(const int *const list, const int n,
                                                  const int nghost,
                                                  const bool nonexclusive) noexcept
{
  double **const centers = atom->x;
  double bbox[6]{};
  cluster_bbox(cspan(list, n), n, bbox);

  if (nonexclusive) {
    bbox[0] = MAX(bbox[0], subbonds[0]);    // cut should be on only needed sides
    bbox[1] = MAX(bbox[1], subbonds[2]);
    bbox[2] = MAX(bbox[2], subbonds[4]);
    bbox[3] = MIN(bbox[3], subbonds[1]);
    bbox[4] = MIN(bbox[4], subbonds[3]);
    bbox[5] = MIN(bbox[5], subbonds[5]);
  }

  const bigint steps_x = static_cast<bigint>((bbox[3] - bbox[0]) / voxel_size) + 1;
  const bigint steps_y = static_cast<bigint>((bbox[4] - bbox[1]) / voxel_size) + 1;
  const bigint steps_z = static_cast<bigint>((bbox[5] - bbox[2]) / voxel_size) + 1;

  if (nloc_grid < steps_x * steps_y * steps_z) {
    nloc_grid = steps_x * steps_y * steps_z;
    occupancy_grid.grow(memory, nloc_grid, "compute_cluster_volume:grid");
  }
  occupancy_grid.reset();

  for (int l = 0; l < n + nghost; ++l) {
    const double *const center = centers[list[l]];
    const auto ni = static_cast<bigint>((center[0] - bbox[0]) / voxel_size);
    const auto nj = static_cast<bigint>((center[1] - bbox[1]) / voxel_size);
    const auto nk = static_cast<bigint>((center[2] - bbox[2]) / voxel_size);
    const bigint min_i = MAX(ni - n_cells, 0);
    const bigint max_i = MIN(ni + n_cells, steps_x - 1);
    const bigint min_j = MAX(nj - n_cells, 0);
    const bigint max_j = MIN(nj + n_cells, steps_y - 1);
    const bigint min_k = MAX(nk - n_cells, 0);
    const bigint max_k = MIN(nk + n_cells, steps_z - 1);

    for (bigint i = min_i; i <= max_i; ++i) {
      for (bigint j = min_j; j <= max_j; ++j) {
        for (bigint k = min_k; k <= max_k; ++k) {
          // Compute the grid cell center coordinates
          // double x = bbox[0] + i * voxel_size + voxel_size / 2;
          // double y = bbox[1] + j * voxel_size + voxel_size / 2;
          // double z = bbox[2] + k * voxel_size + voxel_size / 2;
          // Check if the voxel center is within the sphere of radius r
          const double dx = bbox[0] + i * voxel_size /*x*/ - center[0];
          const double dy = bbox[1] + j * voxel_size /*y*/ - center[1];
          const double dz = bbox[2] + k * voxel_size /*z*/ - center[2];
          if (dx * dx + dy * dy + dz * dz <= overlap_sq) {
            occupancy_grid[::idx(i, j, k, steps_y, steps_z)] = true;
          }
        }
      }
    }
  }

  // count occupied cells
  bigint occupied = 0;
  for (bigint i = 0; i < steps_x; ++i) {
    for (bigint j = 0; j < steps_y; ++j) {
      for (bigint k = 0; k < steps_z; ++k) {
        if (occupancy_grid[::idx(i, j, k, steps_y, steps_z)]) { ++occupied; }
      }
    }
  }

  return occupied * voxel_size * voxel_size * voxel_size;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeClusterVolume::memory_usage()
{
  std::size_t num = volumes.memory_usage() + dist_global.memory_usage() + dist_local.memory_usage() + dist.memory_usage();
  num += recv_buf.memory_usage();
  if (mode == VOLUMEMODE::CALC) {
    num += occupancy_grid.memory_usage();
    if (precompute) { num += offsets.memory_usage(); }
  } else if (mode == VOLUMEMODE::RECTANGLE)
  {
    num += bboxes.memory_usage();
  }
  return static_cast<double>(num);
}
