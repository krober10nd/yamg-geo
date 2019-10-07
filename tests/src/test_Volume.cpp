/*
 calculate the total volume of a tetrahedral mesh that's 10 by 10 x 10 elements in each dimension 
*/
#include <Yamg.h> 

#include <vtkTetra.h>
#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <cassert> 
#include <iostream>
#include <cmath>
#include <vector> 

// subtract two 3-dimensional vectors
template <typename T> 
std::vector<T> subtract(const std::vector<T> &a, const std::vector<T> &b) {
   	std::vector<T> c; 
        c.reserve(a.size()); 
        for(size_t i = 0; i < a.size(); ++i) {
           c.push_back(a[i]-b[i]);
	}
        return c; 
}


// calculate the determinant of a 3x3 matrix
double det_3x3(std::vector<double> e1, std::vector<double> e2, std::vector<double> e3) {
    return(e1[0] * (e2[1] * e3[2] - e3[1]*e2[2]) -
           e1[1] * (e2[0] * e3[2] - e3[0]*e2[2]) +
           e1[2] * (e2[0] * e3[1] - e3[0]*e2[1]));
}


// calculate the volume of a tetrahedron 
double tet_volume(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d) {
    // volume = 1/6 * abs( det( a-d, b-d, c-d) ) 
    auto e1 = subtract(a,d);
    auto e2 = subtract(b,d);
    auto e3 = subtract(c,d);
    double m=det_3x3(e1,e2,e3);
    return (std::abs(m/6.0));
}

// load in the mesh and calculate the volume of all the tetraheadron
int main()
{
  
  std::string basename("box10x10x10.vtu");
  
  // call the mesher 
  system("./bin/yamg-vp --help");

  return 0; 

  vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

  reader->SetFileName(basename.c_str());
  reader->Update();
  ug->DeepCopy(reader->GetOutput());

  int NNodes, NElements;
  int nloc, ndims;
  int cell_type;
  std::vector<double> x, y, z;
  std::vector<int> ENList;

  NNodes = ug->GetNumberOfPoints();
  NElements = ug->GetNumberOfCells();
  
  cell_type = ug->GetCell(0)->GetCellType();
   
  if(cell_type==VTK_TETRA) {
      nloc = 4;
      ndims = 3;
  } else {
      std::cerr<<"ERROR("<<__FILE__<<"): unsupported element type\n";
      exit(-1);
  }

  x.reserve(NNodes);
  y.reserve(NNodes);
  z.reserve(NNodes);
  ENList.reserve(nloc*NElements);
  
  for(size_t i=0; i<NElements; i++) {
     vtkCell *cell = ug->GetCell(i);
     assert(cell->GetCellType()==cell_type);
     for(size_t j=0; j<nloc; j++) {
        ENList.push_back(cell->GetPointId(j));
        }
   }
    
   for(size_t i=0; i<NNodes; i++) {
      double r[3];
      ug->GetPoints()->GetPoint(i, r);
      x.push_back(r[0]);
      y.push_back(r[1]);
      z.push_back(r[2]);
   }

  cout << "INFO: this mesh has " << NNodes << " nodes"<< endl; 
  cout << "INFO: this mesh has " << NElements << " elements"<< endl; 
    
  double volume;
  double total_volume; 
  total_volume = 0.0;

  for(size_t i=0; i<NElements; i++) {
     int ai=ENList[i*4];
     int bi=ENList[i*4+1];
     int ci=ENList[i*4+2];
     int di=ENList[i*4+3];
     
     std::vector<double> a = {x[ai],y[ai],z[ai]}; 
     std::vector<double> b = {x[bi],y[bi],z[bi]}; 
     std::vector<double> c = {x[ci],y[ci],z[ci]}; 
     std::vector<double> d = {x[di],y[di],z[di]}; 

     volume = tet_volume(a,b,c,d);  

     total_volume += volume ;
  }

  cout << "INFO: this mesh has " << total_volume << endl; 

  return 0;

}


