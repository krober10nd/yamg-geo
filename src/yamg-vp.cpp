#include "Yamg.h"

int main(int argc, char **argv)
{
    // Initialise MPI and P4EST.
    Yamg mesher(&argc, &argv);

    // and read in vp from segy file (onto all cores)
    mesher.parse_arguments(&argc, &argv, mesher);

    mesher.init_domain();

    // Refine p4est to geometry
    int lb_its = 3; // number of times to load balance
    mesher.refine_p4est(Yamg::refine_fn,lb_its);

    // Generate tetrahedra from quads 
    mesher.triangulate();

    mesher.sanity_mesh(__FILE__, __LINE__);

    mesher.write_vtu(mesher.ofname);

    mesher.write_gmsh(mesher.ofname);

    mesher.dump_stats("this is the end");

    mesher.finalise();

    return 0;
}
