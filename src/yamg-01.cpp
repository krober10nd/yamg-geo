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
             <<" -r <value>, --resolution <value>\n\tSet minimum feature resolution.\n";
    return;
}

double __yamg_feature_resolution=-1;
double __yamg_scale=200.0;

void read_segmented_file(Yamg &mesher, std::string filename)
{
    std::ifstream file(filename);
    std::vector<float> data;
    std::copy(std::istream_iterator<int>(file), {}, std::back_inserter(data));
    file.close();

    double spacing[] = {20, 10, 2};
    int dims[] = {60, 220, 35};
    double origin[] = {0, 0, 0};

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
        {0, 0, 0, 0}
    };

    int optionIndex = 0;
    int c;
    const char *shortopts = "hvdr:";

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
        read_segmented_file(mesher, filename);
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

        float sum=0;
        int tcnt=0;
        for(int i=index_range[0]; i<=index_range[1]; i++) {
            for(int j=index_range[2]; j<=index_range[3]; j++) {
                for(int k=index_range[4]; k<=index_range[5]; k++) {
                    sum += Yamg::get_scalar(i, j, k);
                    tcnt++;
                }
            }
        }

        if (tcnt==0) {
            return 0;
        }

        float vf = sum/tcnt, var=0;
        for(int i=index_range[0]; i<index_range[1]; i++) {
            for(int j=index_range[2]; j<index_range[3]; j++) {
                for(int k=index_range[4]; k<index_range[5]; k++) {
                    var+=std::pow(Yamg::get_scalar(i, j, k)-vf, 2);
                }
            }
        }
        var/=tcnt;

        // Refine if the variance is "high".
        if(var>0.05) {
            return 1;
        }

        return 0;
    }
}

int main(int argc, char **argv)
{
    // Initialise MPI and P4EST.
    Yamg mesher(&argc, &argv);

    parse_arguments(&argc, &argv, mesher);

    std::string basename("basename");

    mesher.write_vti(basename);

    mesher.init_domain();

    // Refine p4est to geometry
    mesher.refine_p4est(refine_fn);

    // Generate tetrahedra
    mesher.triangulate();

    mesher.sanity_mesh(__FILE__, __LINE__);

    mesher.write_vtu(basename);
    mesher.write_gmsh(basename);

    mesher.dump_stats("this is the end");

    mesher.finalise();

    return 0;
}
