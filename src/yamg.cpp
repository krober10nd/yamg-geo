#include "ReservoirMesher.h"
#include <getopt.h>

#include <iostream>
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

int parse_arguments(int *argc, char ***argv, ReservoirMesher &mesher)
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
            mesher.set_feature_resolution(atof(optarg));
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

    std::string vel_filename((*argv)[*argc-1]);
    struct stat buffer;   
    if(stat (vel_filename.c_str(), &buffer) == 0) {
        mesher.set_vel_file(std::string((*argv)[*argc-1]));
    }else{
        std::cerr<<"ERROR: missing velocity file"<<std::endl;
        usage(*argv[0]);
        exit(-1);
    }

    return 0;
}

int main(int argc, char **argv)
{
    // Initialise MPI and P4EST.
    ReservoirMesher mesher(&argc, &argv);

    parse_arguments(&argc, &argv, mesher);

    std::string basename("basename");

    mesher.init_domain();

    // Refine p4est to geometry
    mesher.refine_p4est();

    // Generate tetrahedra
    mesher.triangulate();

    mesher.sanity_mesh(__FILE__, __LINE__);

    mesher.write_vtu(basename);
    mesher.write_gmsh(basename);

    mesher.dump_stats("this is the end");

    mesher.finalise();

    return 0;
}
