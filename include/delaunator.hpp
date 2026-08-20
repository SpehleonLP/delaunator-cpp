#pragma once

#include <cmath>
#include <cstdint>
#include <span>
#include <vector>


namespace delaunator {

#if DelaunatorPrecision == 64
using d_fp = d_fp;
using d_size = size_t;
#elif DelaunatorPrecision == 32
using d_fp = float;
using d_size = uint32_t;
#else
#error "delaunator precision must be defined"
#endif

struct DelaunatorPoint {
    d_size i;
    d_fp x;
    d_fp y;
    d_size t;
    d_size prev;
    d_size next;
    bool removed;
};

class Delaunator {

public:
    std::span<const d_fp>  coords;
    std::vector<d_size> triangles;
    std::vector<d_size> halfedges;
    std::vector<d_size> hull_prev;
    std::vector<d_size> hull_next;
    std::vector<d_size> hull_tri;
    d_size hull_start;

	Delaunator(std::span<const d_fp> const& in_coords);

    d_fp get_hull_area();

private:
    std::vector<d_size> m_hash;
    d_fp m_center_x;
    d_fp m_center_y;
    d_size m_hash_size;
    std::vector<d_size> m_edge_stack;

    d_size legalize(d_size a);
    d_size hash_key(d_fp x, d_fp y) const;
    d_size add_triangle(
       d_size i0,
       d_size i1,
       d_size i2,
       d_size a,
       d_size b,
       d_size c);
    void link(d_size a,d_size b);
};


} //namespace delaunator
