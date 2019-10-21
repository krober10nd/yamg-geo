/*
Test basic functionality of yamg-geo 
*/
#include <Yamg.h> 

int main(int argc, char *argv[])
{

  Yamg mesher(&argc, &argv);

  mesher.parse_arguments(&argc, &argv, mesher);

  mesher.init_domain();

  int lb_its = 3; // number of times to load balance
  mesher.refine_p4est(Yamg::refine_fn,lb_its);

  mesher.triangulate();

  mesher.write_vtu(mesher.ofname);

  mesher.testMeshGeometry();

  mesher.finalise();

  return 0 ; 

}


