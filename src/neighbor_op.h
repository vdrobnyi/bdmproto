#ifndef NEIGHBOR_OP_H_
#define NEIGHBOR_OP_H_

#include <iostream>
#include <math.h>

#include "spatial/octree.hpp"
#include "backend.h"
#include "array.h"

namespace bdm {

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
            spatial_tree_node<size_t> *tree = new octree_node<size_t>(100, 100);
            size_t *neighbor_counter = new size_t[cells->elements()];

            for (size_t i = 0; i < cells->elements(); i++) {
                neighbor_counter[i] = 0;
                auto cell = cells->GetScalar(i);
                const auto &position = cell.GetPosition();
                point pos(position[0][0], position[1][0], position[2][0]);
                tree->put(pos, i);
            }
            const VcBackend::real_t search_radius =
                    sqrt(static_cast<VcBackend::real_t>(distance_));

            bdm::array<int, 8> * i_neighbors = new bdm::array<int, 8>[cells->elements()];

            std::cout << "Neighbor search. Distance " << search_radius << std::endl;
            auto tree_neighbors = tree->get_neighbors(search_radius);
            int amount_of_pairs = tree_neighbors->size();

            for (int i = 0; i < amount_of_pairs; i++) {
                size_t neighbor_a = (*tree_neighbors)[i].first;
                size_t neighbor_b = (*tree_neighbors)[i].second;
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

            std::cout << "End of search." << std::endl;
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
