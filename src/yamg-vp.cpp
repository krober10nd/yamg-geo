#include "Yamg.h"
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <sys/stat.h>

void usage(char *cmd)
{
    std::cout<<"Usage: "<<cmd<<" <vel-file> [options]\n"
             <<"\nOptions:\n"
             <<" -h, --help\n\tHelp! Prints this message.\n"
             <<" -v, --verbose\n\tVerbose output.\n"
             <<" -d, --debug\n\tDebug mode (slow!).\n"
             <<" -r <value>, --resolution <value>\n\tSet minimum feature resolution.\n"
             <<" -s <value>, --scale <value>\n\tSet scale length.\n";
    return;
}

double __yamg_feature_resolution=-1;
double __yamg_scale=200.0;

double vp_sw=1484.0; // Velocity of sound in sea water

// marked for deletion, will use python to accomplish this
void read_velocity_file(Yamg &mesher, std::string filename)
{
    ///*
    //    scaling info:
    //    H0.vel contains velocities vp in m/s
    //    a length scale is typically  scale*vp/1500
    //    with scale, for instance, 100, 50, 25, 10 or 5 m.

    //    unformatted file, little-endian, with
    //     3 doubles for ox (origin in x,y,z),
    //     3 doubles for dx (spacing in x,y,z),
    //     3 4-byte ints for nx (gridpoints in x,y,z),
    //     nx(1)*nx(2)*nx(3) 4-byte floats for velocity vp,
    //       x is fastest index, z is slowest
    //*/
    //
    //

    std::ifstream input(filename.c_str(), std::ios::binary);
    std::vector<float> data;
    double origin[3], spacing[3];
    int dims[3];
    input.read((char *)origin, 3*sizeof(double));
    input.read((char *)spacing, 3*sizeof(double));
    input.read((char *)dims, 3*sizeof(int));

    data.resize(dims[0]*dims[1]*dims[2]);
    input.read((char *)data.data(), data.size()*sizeof(float));

    input.close();

    for(auto &val:data) {
        if(!std::isfinite(val))
        val = vp_sw;
    }

    mesher.set_data(data, origin, spacing, dims);

    if(__yamg_feature_resolution==-1) {
        __yamg_feature_resolution = std::min(spacing[0], std::min(spacing[1], spacing[2]));
    }
}

int parse_arguments(int *argc, char ***argv, Yamg &mesher)
{

    // Set defaults
    if(*argc==1) {
        usage(*argv[0]);
        exit(-1);
    }

    struct option longOptions[] = {
        {"help", 0, 0, 'h'},
        {"verbose", 0, 0, 'v'},
        {"debug", 0, 0, 'd'},
        {"resolution", optional_argument, 0, 'r'},
        {"scale", optional_argument, 0, 's'},
        {0, 0, 0, 0}
    };

    int optionIndex = 0;
    int c;
    const char *shortopts = "hvdr:s:";

    // Set opterr to nonzero to make getopt print error messages
    opterr=1;
    while (true) {
        c = getopt_long(*argc, *argv, shortopts, longOptions, &optionIndex);

        if (c == -1) break;

        switch (c) {
        case 'h':
            usage(*argv[0]);
            break;
        case 'v':
            mesher.enable_verbose();
            break;
        case 'd':
            mesher.enable_debugging();
            break;
        case 'r':
            __yamg_feature_resolution = atof(optarg);
            break;
        case 's':
            __yamg_scale = atof(optarg);
            break;
        case '?':
            // missing argument only returns ':' if the option string starts with ':'
            // but this seems to stop the printing of error messages by getopt?
            std::cerr<<"ERROR: unknown option or missing argument\n";
            usage(*argv[0]);
            exit(-1);
        case ':':
            std::cerr<<"ERROR: missing argument\n";
            usage(*argv[0]);
            exit(-1);
        default:
            // unexpected:
            std::cerr<<"ERROR: getopt returned unrecognized character code\n";
            exit(-1);
        }
    }

    std::string filename((*argv)[*argc-1]);
    struct stat buffer;   
    if(stat (filename.c_str(), &buffer) == 0) {
      read_velocity_file(mesher, filename);
    }else{
        std::cerr<<"ERROR: missing velocity file"<<std::endl;
        usage(*argv[0]);
        exit(-1);
    }

    return 0;
}

extern "C" {
    int refine_fn(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
    {
        double quad[3*8];
        int index_range[6];
        Yamg::generate_quadcoords(p4est, which_tree, quadrant, quad, NULL, index_range);

        // Check if we have already reached the smallest feature size.
        double maxl=pow(quad[0]-quad[3], 2);
        maxl = std::max(maxl, pow(quad[1]-quad[2*3+1], 2));
        maxl = std::max(maxl, pow(quad[2]-quad[4*3+2], 2));
        maxl = sqrt(maxl);

        if(maxl<=__yamg_feature_resolution) {
            return 0;
        }

        for(int i=index_range[0]; i<=index_range[1]; i++) {
            for(int j=index_range[2]; j<=index_range[3]; j++) {
                for(int k=index_range[4]; k<=index_range[5]; k++) {
                    float vp = Yamg::get_scalar(i, j, k);

                    if(std::isfinite(vp)) {
                        float l = __yamg_scale*vp/1500.;
                        if(l<maxl) {
                            return 1;
                        }
                    }
                }
            }
        }

        return 0;
    }
}

int main(int argc, char **argv)
{
    // Initialise MPI and P4EST.
    Yamg mesher(&argc, &argv);

    // and read in vp from segy file
    parse_arguments(&argc, &argv, mesher);

    std::string basename("basename");

    //mesher.write_vti(basename);

    mesher.init_domain();

    // Refine p4est to geometry
    mesher.refine_p4est(refine_fn);

    // Generate tetrahedra from quads 
    mesher.triangulate();

    mesher.sanity_mesh(__FILE__, __LINE__);

#ifdef HAVE_VTK
    mesher.write_vtu(basename);
<<<<<<< HEAD
#endif
=======

>>>>>>> 2d73ac2... Implementation of CI building with Travis.
    mesher.write_gmsh(basename);

    mesher.dump_stats("this is the end");

    mesher.finalise();

    return 0;
}
