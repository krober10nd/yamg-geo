/*
Test basic functionality of yamg-geo 
 ??does the application run and produce a valid mesh??
*/
#include <Yamg.h> 

#include <vtkTetra.h>
#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <string>
#include <cassert> 
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector> 

double yamg_feature_resolution=-1; // default is minimum of nx[3]
double yamg_scale=200.0; // speed of sound in sea water gets this resolution in m 
double vp_sw=1484.0; // Velocity of sound in sea water

// call back function for p4est to refine octrees
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

        if(maxl<=yamg_feature_resolution) {
            return 0;
        }

        for(int i=index_range[0]; i<=index_range[1]; i++) {
            for(int j=index_range[2]; j<=index_range[3]; j++) {
                for(int k=index_range[4]; k<=index_range[5]; k++) {
                    float vp = Yamg::get_scalar(i, j, k);

                    if(std::isfinite(vp)) {
                        float l = yamg_scale*vp/vp_sw;
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


int main(int argc, char *argv[])
{

  // call the mesher 
  Yamg mesher(&argc, &argv);

  // and read in vp from segy file
  mesher.parse_arguments(&argc, &argv, mesher);

  // TODO: enable user to pass name of desired mesh file
  std::string basename("UnitCube");

  mesher.init_domain();

  // Refine p4est to geometry
  mesher.refine_p4est(refine_fn);

  // Generate tetrahedra from quads 
  mesher.triangulate();

  mesher.write_vtu(basename);

  mesher.finalise();

  mesher.dump_stats("this is the end");

  return 0 ; 
  

}


