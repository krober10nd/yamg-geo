#ifndef YAMG_H
#define YAMG_H

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
#include <cmath> 
#include <string>

class Yamg
{
public:
    // default constructor 
    Yamg(int *argc, char ***argv);

    // default destructor 
    ~Yamg();
    
    static int refine_fn(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant);

    int parse_arguments(int *argc, char ***argv, Yamg &mesher);

    void read_velocity_file(Yamg &mesher, std::string filename);

    void usage(char *cmd);

    static void append_cell(int n0, int n1, int n2, int n3,
                            int b0, int b1, int b2, int b3);

    void build_facets(std::vector<int> &facets, std::vector<int> &facet_ids);

    void dump_stats(std::string tag);

    void enable_debugging();
    void enable_verbose();

    void finalise();

    static void generate_quadcoords(const p4est_t *p4est, p4est_topidx_t treeid, const p4est_quadrant_t *quad, double *x, int *boundary_id, int *index_range);

    const double *get_point(int nid) const;

    float get_scalar(double x, double y, double z) const;
    static float get_scalar(int i, int j, int k);
    float get_scalar_p0(int eid) const;
    static const int *get_tet(int eid);

    void init_domain();

    void refine_p4est(p4est_refine_t refine_fn, int recursive);

    void sanity_mesh(const char *file=NULL, int line=0) const;
    void set_data(const std::vector<float> &data_in, const double *origin, const double *spacing, const int *dims);

    static long double tet_volume(int e);

    double total_volume() const;
    double total_area() const;

    void triangulate();

    void write_vtu(std::string basename);
    void write_vti(std::string basename) const; 
    void write_gmsh(std::string basename);

    std::string ofname="YamgMesh"; // output mesh filename

private:

    // MPI communicator
    sc_MPI_Comm mpicomm;

    // Create the initial p4est.
    void create_p4est();

    void get_element_facet(int eid, int fid, int facet[3]) const;

    void get_lnn2unn(std::vector<int> &_lnn2unn);

    int get_number_points();
    int get_number_tets();

    // Initialise coordinate hash.
    int init_coord_hash();

    void sanity_facet(const int *facet) const;

    std::string filename_vtu(const char *_filename);

    // p4est data.
    p4est_t *p4est;
    p4est_connectivity_t *conn;
    p4est_ghost_t *ghost;

    bool debugging;
    bool verbose;
    int nranks, rank;
};

#endif
