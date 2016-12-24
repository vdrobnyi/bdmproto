#ifndef NEIGHBOR_OP_H_
#define NEIGHBOR_OP_H_

#include <iostream>
#include <nanoflann.h>

//#include "spatial/kd_tree.h"
#include "spatial/octree.hpp"
#include "backend.h"
#include "array.h"

namespace bdm {

    using nanoflann::KDTreeSingleIndexAdaptorParams;
    using nanoflann::L2_Simple_Adaptor;
    using nanoflann::KDTreeSingleIndexAdaptor;

// https://github.com/jlblancoc/nanoflann/blob/master/examples/pointcloud_adaptor_example.cpp
// And this is the "dataset to kd-tree" adaptor class:
    template<typename Derived>
    struct NanoFlannDaosoaAdapter {
        using coord_t = VcBackend::real_t;

        const Derived &obj;  //!< A const ref to the data set origin

        /// The constructor that sets the data set source
        NanoFlannDaosoaAdapter(const Derived &obj_) : obj(obj_) {}

        /// CRTP helper method
        inline const Derived &derived() const { return obj; }

        /// Must return the number of data points
        inline size_t kdtree_get_point_count() const { return derived().elements(); }

        /// Returns the distance between the vector "p1[0:size-1]" and the data point
        /// with index "idx_p2" stored in the class:
        inline coord_t kdtree_distance(const coord_t *p1, const size_t idx_p2,
                                       size_t /*size*/) const {
            // fixme makes a lot of copies (GetScalar)
            const coord_t d0 = p1[0] - derived().GetScalar(idx_p2).GetPosition()[0][0];
            const coord_t d1 = p1[1] - derived().GetScalar(idx_p2).GetPosition()[1][0];
            const coord_t d2 = p1[2] - derived().GetScalar(idx_p2).GetPosition()[2][0];
            return d0 * d0 + d1 * d1 + d2 * d2;
        }

        /// Returns the dim'th component of the idx'th point in the class:
        /// Since this is inlined and the "dim" argument is typically an immediate
        /// value, the "if/else's" are actually solved at compile time.
        inline coord_t kdtree_get_pt(const size_t idx, int dim) const {
            // fixme makes a lot of copies (GetScalar)
            return derived().GetScalar(idx).GetPosition()[dim][0];
        }

        /// Optional bounding-box computation: return false to default to a standard
        /// bbox computation loop.
        ///   Return true if the BBOX was already computed by the class and returned
        ///   in
        ///   "bb" so it can be avoided to redo it again.
        ///   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3
        ///   for point clouds)
        template<class BBOX>
        bool kdtree_get_bbox(BBOX & /*bb*/) const {
            return false;
        }
    };

    class NeighborOp {
    public:
        NeighborOp() {}

        NeighborOp(double distance) : distance_(distance) {}

        ~NeighborOp() {}

        NeighborOp(const NeighborOp &) = delete;

        NeighborOp &operator=(const NeighborOp &) = delete;

        template<typename daosoa>
        Vc_ALWAYS_INLINE void Compute(daosoa *cells) const {

            std::vector<std::array<bdm::array<int, 8>, VcBackend::kVecLen> > neighbors(
                    cells->vectors());


            // Tree search
            spatial_tree_node<size_t> *tree = new octree_node<size_t>(bound(-1000, -1000, -1000, 1000, 1000, 1000), 100,
                                                                      100);
            size_t *neighbor_counter = new size_t[cells->elements()];

            for (size_t i = 0; i < cells->elements(); i++) {
                neighbor_counter[i] = 0;
                auto cell = cells->GetScalar(i);
                const auto &position = cell.GetPosition();
                const VcBackend::real_t query_pt[3] = {};
                point pos(position[0][0], position[1][0], position[2][0]);
                tree->put(pos, i);
            }
            const VcBackend::real_t search_radius =
                    static_cast<VcBackend::real_t>(distance_);

            bdm::array<int, 8> * i_neighbors = new bdm::array<int, 8>[cells->elements()];

            auto tree_neighbors = tree->get_neighbors_with_points(25);
            int amount_of_pairs = tree_neighbors->size();

            for (int i = 0; i < amount_of_pairs; i++) {
                size_t neighbor_a = (*tree_neighbors)[i].first.second;
                size_t neighbor_b = (*tree_neighbors)[i].second.second;
                if (neighbor_counter[neighbor_a] < 8)
                    i_neighbors[neighbor_a][neighbor_counter[neighbor_a]++] = neighbor_b;
                if  (neighbor_counter[neighbor_b] < 8)
                    i_neighbors[neighbor_b][neighbor_counter[neighbor_b]++] = neighbor_a;
            }
            for (size_t i = 0; i < cells->elements(); i++) {
                const auto vector_idx = i / VcBackend::kVecLen;
                const auto scalar_idx = i % VcBackend::kVecLen;
                neighbors[vector_idx][scalar_idx] = std::move(i_neighbors[i]);
                neighbors[vector_idx][scalar_idx].SetSize(neighbor_counter[i]);
            }
            delete tree_neighbors;
            delete tree;
            delete[] neighbor_counter;
            delete[] i_neighbors;
            // End of tree search

            // update neighbors
#pragma omp parallel for
            for (size_t i = 0; i < cells->vectors(); i++) {
                (*cells)[i].SetNeighbors(neighbors[i]);
            }
        }

    private:
        double distance_ = 3000;
    };

}  // namespace bdm

#endif  // NEIGHBOR_OP_H_
