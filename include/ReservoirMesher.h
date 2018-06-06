#ifndef RESERVOIRMESHER_H
#define RESERVOIRMESHER_H

#include <mpi.h>

#include <p4est_to_p8est.h>
#include <p8est_vtk.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_iterate.h>
#include <p8est_bits.h>

#include <vector>
#include <map>
#include <set>
#include <unordered_set>

class ReservoirMesher
{
public:
    ReservoirMesher(int *argc, char ***argv);
    ~ReservoirMesher();

    void build_facets(std::vector<int> &facets, std::vector<int> &facet_ids);
    void set_vel_file(std::string filename);

    void finalise();

    void init_domain();

    bool is_rank_master();

    void refine_p4est(int lb_its=3);

    void enable_debugging();
    void enable_verbose();
    void set_feature_resolution(float);

    void triangulate();
    void construct_adjancies();

    void dump_stats(std::string tag);

    double total_volume() const;
    double total_area() const;

    void write_vtu(std::string basename);
    void write_pvtu(std::string basename);
    void write_gmsh(std::string basename);

    void sanity_mesh(const char *file=NULL, int line=0) const;
private:
    // MPI communicator
    sc_MPI_Comm mpicomm;

    // Initialise coordinate hash.
    int init_coord_hash();

    void get_element_facet(int eid, int fid, int facet[3]) const;

    // Create the initial p4est.
    void create_p4est();

    void sanity_facet(const int *facet) const;

    int get_number_points();
    int get_number_tets();

    void get_lnn2unn(std::vector<int> &_lnn2unn);

    std::string filename_vtu(const char *_filename);

    // p4est data.
    p4est_t *p4est;
    p4est_connectivity_t *conn;
    p4est_ghost_t *ghost;

    std::vector<int> boundary;
    std::map<std::string, int> label_boundary_map;

    bool debugging;
    bool verbose;
    int nranks, rank;
};

#endif
