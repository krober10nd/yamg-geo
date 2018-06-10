#include "Yamg.h"

#include <vtkTetra.h>
#include <vtkIdList.h>
#include <vtkTriangle.h>
#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkLongLongArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>

long long _hash_stride[2];
bool hash_trusty;
double hash_resolution;

double geom_bounds[6];
double ox[3], dx[3];
int nx[3];
std::vector<float> data;

std::unordered_map<long long, int> unn2lnn;
std::vector<long long> lnn2unn;

std::vector<double> gcoords;

std::vector<int> gcells;
std::vector<int> boundary;

const double *Yamg::get_point(int nid) const
{
    assert(nid>=0);
    assert(nid<gcoords.size()/3);
    return gcoords.data()+nid*3;
}

const int *Yamg::get_tet(int eid)
{
    assert(eid>=0);
    assert(eid<gcells.size()/4);
    return gcells.data()+eid*4;
}

// Calculate volume of tetrahedron.
long double Yamg::tet_volume(int e)
{
    assert(e>=0);
    assert(e<=gcells.size()/4);

    const int *n=Yamg::get_tet(e);

    long double x0[3], x1[3], x2[3];
    for(int i=0; i<3; i++) {
        x0[i] = gcoords[n[1]*3+i] - gcoords[n[0]*3+i];
        x1[i] = gcoords[n[2]*3+i] - gcoords[n[0]*3+i];
        x2[i] = gcoords[n[3]*3+i] - gcoords[n[0]*3+i];
    }

    long double det = x0[0]*x1[1]*x2[2] - x0[0]*x1[2]*x2[1] - x0[1]*x1[0]*x2[2] + x0[1]*x1[2]*x2[0] + x0[2]*x1[0]*x2[1] - x0[2]*x1[1]*x2[0];

    // Sign based on choice of orientation.
    long double v = -det/6;

    return v;
}

void Yamg::append_cell(int n0, int n1, int n2, int n3,
                       int b0, int b1, int b2, int b3) {
     gcells.push_back(n0);
     gcells.push_back(n1);
     gcells.push_back(n2);
     gcells.push_back(n3);

     boundary.push_back(b0);
     boundary.push_back(b1);
     boundary.push_back(b2);
     boundary.push_back(b3);

     assert(gcells.size()==boundary.size());
     assert(Yamg::tet_volume(gcells.size()/4-1)>0);
}

float Yamg::get_scalar(int i, int j, int k)
{
    i = std::min(i, nx[0]-1);
    j = std::min(j, nx[1]-1);
    k = std::min(k, nx[2]-1);

    return data[k*nx[0]*nx[1]+j*nx[0]+i];
}

float Yamg::get_scalar(double x, double y, double z) const
{
    int i = (x-ox[0])/dx[0];
    int j = (y-ox[1])/dx[1];
    int k = (z-ox[2])/dx[2];

    return get_scalar(i, j, k);
}

// This is really shit - should be projected value from image to tetrahedron
float Yamg::get_scalar_p0(int eid) const{
    double value = 0;
    int cnt=0;
    const int *tet = Yamg::get_tet(eid);
    for(int i=0; i<4; i++) {
        const double *x = get_point(tet[i]);
        value += get_scalar(x[0], x[1], x[2]);
        cnt++;
    }
    return value/cnt;
}

void Yamg::set_data(const std::vector<float> &data_in, const double *origin, const double *spacing, const int *dims)
{
    data = data_in;
    for(int i=0;i<3;i++) {
        ox[i] = origin[i];
        dx[i] = spacing[i];
        nx[i] = dims[i];
    }

    for(int i=0; i<3; i++) {
        geom_bounds[2*i] = ox[i]+dx[i];
        geom_bounds[2*i+1] = ox[i]+dx[i]*(nx[i]-1);
    }

    for(int i=0; i<2; i++)
        _hash_stride[i] = 0;

    hash_resolution = 0.01*std::min(std::min(dx[0], dx[1]), dx[2]);
    hash_trusty=false;
}

double edge_length(const double *x0, const double *x1)
{
    return sqrt((x0[0] - x1[0])*(x0[0] - x1[0])+
                (x0[1] - x1[1])*(x0[1] - x1[1])+
                (x0[2] - x1[2])*(x0[2] - x1[2]));
}

// Calculate area using Heron's Formula
long double tri_area(const double *x1, const double *x2, const double *x3) {
    long double a;
    {
        long double dx = ((long double)x1[0]-(long double)x2[0]);
        long double dy = ((long double)x1[1]-(long double)x2[1]);
        long double dz = ((long double)x1[2]-(long double)x2[2]);
        a = std::sqrt(dx*dx+dy*dy+dz*dz);
    }
    long double b;
    {
        long double dx = ((long double)x1[0]-(long double)x3[0]);
        long double dy = ((long double)x1[1]-(long double)x3[1]);
        long double dz = ((long double)x1[2]-(long double)x3[2]);
        b = std::sqrt(dx*dx+dy*dy+dz*dz);
    }
    long double c;
    {
        long double dx = ((long double)x2[0]-(long double)x3[0]);
        long double dy = ((long double)x2[1]-(long double)x3[1]);
        long double dz = ((long double)x2[2]-(long double)x3[2]);
        c = std::sqrt(dx*dx+dy*dy+dz*dz);
    }
    long double s = (a+b+c)/2;

    return std::sqrt(s*(s-a)*(s-b)*(s-c));
}

void Yamg::generate_quadcoords(const p4est_t *p4est, p4est_topidx_t treeid, const p4est_quadrant_t *quad, double *x, int *boundary_id, int *index_range)
{
    assert(P4EST_CHILDREN==8); // not pretending this is general

    p4est_quadrant_t node;
    for (int i=0; i<P4EST_CHILDREN; ++i) {

        p4est_quadrant_corner_node (quad, i, &node);
        p4est_qcoord_to_vertex (p4est->connectivity, treeid, node.x, node.y, node.z, x+i*3);
    }

    if(boundary_id!=NULL) {
        // Label boundary
        for (int i=0; i<6; ++i) {
            boundary_id[i] = -1;
        }

        for(int face=0;face<6;face++) {
            p4est_quadrant_t r;
            p4est_topidx_t neighbor = p4est_quadrant_face_neighbor_extra(quad, treeid, face, &r, NULL, p4est->connectivity);

            if(neighbor==-1) {
                boundary_id[face] = face;
            }
        }
    }

    if(index_range!=NULL) {
        // Index ranges for data
        index_range[0] = (x[0]-ox[0])/dx[0];
        index_range[1] = (x[7*3+0]-ox[0])/dx[0];

        index_range[2] = (x[1]-ox[1])/dx[1];
        index_range[3] = (x[7*3+1]-ox[1])/dx[1];

        index_range[4] = (x[2]-ox[2])/dx[2];
        index_range[5] = (x[7*3+2]-ox[2])/dx[2];
    }
}

// Give a unique identifier for a coordinate
long long coord_hash(double x, double y, double z)
{
    assert(hash_trusty);

    long long ix = (long long)std::round((x-geom_bounds[0])/hash_resolution);
    assert(ix>=0);

    long long iy = (long long)std::round((y-geom_bounds[2])/hash_resolution);
    assert(iy>=0);

    long long iz = (long long)std::round((z-geom_bounds[4])/hash_resolution);
    assert(iz>=0);

    long long unn = iz + iy*_hash_stride[1] + ix*_hash_stride[0];

    assert(unn>=0);

    return unn;
}

/* Take a brick, look at its neighbours and figure out if there are split edges, triangulate.
 */
void mesh_quad (p4est_iter_volume_info_t * info, void *user_data)
{
    double quad[3*8];
    int _boundary[6];
    Yamg::generate_quadcoords(info->p4est, info->treeid, info->quad, quad, _boundary, NULL);

    // Create a local table from local node numbers to unique (global) node number.
    long long unn[8];
    int lnn[8];
    std::unordered_map<int, long long> _lnn2unn;
    for(int i=0; i<8; i++) {
        unn[i] = coord_hash(quad[i*3], quad[i*3+1], quad[i*3+2]);
        assert(unn2lnn.find(unn[i])!=unn2lnn.end());

        lnn[i] = unn2lnn[ unn[i] ];
        _lnn2unn.insert({lnn[i], unn[i]});
    }

    std::vector<int> cavity(12*6, -1), cavity_boundary(12);

    // Two triangular facets on the -x side of brick
    cavity[0] = lnn[0];
    cavity[1] = lnn[2];
    cavity[2] = lnn[6];
    assert(gcoords[lnn[6]*3+1]-gcoords[lnn[0]*3+1] + gcoords[lnn[6]*3+2]-gcoords[lnn[0]*3+2] > 0);

    cavity[6] = lnn[0];
    cavity[7] = lnn[6];
    cavity[8] = lnn[4];

    cavity_boundary[0] = _boundary[0];
    cavity_boundary[1] = _boundary[0];

    // Two triangular facets on +x
    cavity[12] = lnn[1];
    cavity[13] = lnn[7];
    cavity[14] = lnn[3];
    assert(gcoords[lnn[7]*3+1]-gcoords[lnn[1]*3+1] + gcoords[lnn[7]*3+2]-gcoords[lnn[1]*3+2] > 0);

    cavity[18] = lnn[1];
    cavity[19] = lnn[5];
    cavity[20] = lnn[7];

    cavity_boundary[2] = _boundary[1];
    cavity_boundary[3] = _boundary[1];

    // -y
    cavity[24] = lnn[0];
    cavity[25] = lnn[5];
    cavity[26] = lnn[1];
    assert(gcoords[lnn[5]*3+0]-gcoords[lnn[0]*3+0] + gcoords[lnn[5]*3+2]-gcoords[lnn[0]*3+2] > 0);

    cavity[30] = lnn[0];
    cavity[31] = lnn[4];
    cavity[32] = lnn[5];

    cavity_boundary[4] = _boundary[2];
    cavity_boundary[5] = _boundary[2];

    // +y
    cavity[36] = lnn[2];
    cavity[37] = lnn[3];
    cavity[38] = lnn[7];
    assert(gcoords[lnn[7]*3+0]-gcoords[lnn[2]*3+0] + gcoords[lnn[7]*3+2]-gcoords[lnn[2]*3+2] > 0);

    cavity[42] = lnn[2];
    cavity[43] = lnn[7];
    cavity[44] = lnn[6];

    cavity_boundary[6] = _boundary[3];
    cavity_boundary[7] = _boundary[3];

    // -z
    cavity[48] = lnn[0];
    cavity[49] = lnn[1];
    cavity[50] = lnn[3];

    cavity[54] = lnn[0];
    cavity[55] = lnn[3];
    cavity[56] = lnn[2];

    cavity_boundary[8] = _boundary[4];
    cavity_boundary[9] = _boundary[4];

    // +z
    cavity[60] = lnn[4];
    cavity[61] = lnn[7];
    cavity[62] = lnn[5];

    cavity[66] = lnn[4];
    cavity[67] = lnn[6];
    cavity[68] = lnn[7];

    cavity_boundary[10] = _boundary[5];
    cavity_boundary[11] = _boundary[5];

    // Check for split edges.
    bool insert_center_point=false;
    for(int i=0; i<12; i++) {
        int *triangle=cavity.data()+i*6;
        for(int j=0; j<3; j++) {
            double x[3];
            for(int k=0; k<3; k++)
                x[k] = (gcoords[triangle[j]*3+k]+gcoords[triangle[(j+1)%3]*3+k])/2;
            long long _unn = coord_hash(x[0], x[1], x[2]);

            auto ilnn = unn2lnn.find(_unn);
            if(ilnn!=unn2lnn.end()) {
                // Edge is actually split.
                triangle[3+j] = ilnn->second;
                _lnn2unn.insert({ilnn->second, _unn});
                insert_center_point=true;
            }
        }
    }

    if(insert_center_point) {
        // Create and append centre point
        double x[]= {(gcoords[lnn[0]*3  ]+gcoords[lnn[7]*3  ])*0.5,
                     (gcoords[lnn[0]*3+1]+gcoords[lnn[7]*3+1])*0.5,
                     (gcoords[lnn[0]*3+2]+gcoords[lnn[7]*3+2])*0.5
                    };

        long long _unn = coord_hash(x[0], x[1], x[2]);
        assert(unn2lnn.find(_unn)==unn2lnn.end());

        int cid = gcoords.size()/3;
        unn2lnn.insert({_unn, cid});
        lnn2unn.push_back(_unn);

        for(int i=0; i<3; i++)
            gcoords.push_back(x[i]);

        // - loop over all facets in the cavity
        // - if one or more edges are split then refine in the right way
        // - conect facet to centre point to form a tet
        for(int i=0; i<12; i++) {
            // Count number of edges cut in this reference triangle.
            int cutcnt=0;
            for(int j=3; j<6; j++) {
                if(cavity[i*6+j]!=-1)
                    cutcnt++;
            }

            if(cutcnt==0) {
                Yamg::append_cell(cavity[i*6+1], cavity[i*6+0], cavity[i*6+2], cid,
                                             -1, -1, -1, cavity_boundary[i]);
            } else if(cutcnt==1) {
                int ref_copy[]= {-1, -1, -1, -1, -1, -1};
                for(int j=0; j<3; j++) {
                    if(cavity[i*6+3+j]!=-1) {
                        for(int k=0; k<3; k++) {
                            ref_copy[k  ] = cavity[i*6+(j+k)%3];
                            ref_copy[k+3] = cavity[i*6+3+(j+k)%3];
                        }
                        break;
                    }
                }

                Yamg::append_cell(ref_copy[3], ref_copy[0], ref_copy[2], cid,
                                             -1, -1, -1, cavity_boundary[i]);

                Yamg::append_cell(ref_copy[1], ref_copy[3], ref_copy[2], cid,
                                             -1, -1, -1, cavity_boundary[i]);
            } else if(cutcnt==2) {
                int ref_copy[]= {-1, -1, -1, -1, -1, -1};
                for(int j=0; j<3; j++) {
                    if(cavity[i*6+3+j]==-1) {
                        j = (j+1)%3;
                        for(int k=0; k<3; k++) {
                            ref_copy[k  ] = cavity[i*6+(j+k)%3];
                            ref_copy[k+3] = cavity[i*6+3+(j+k)%3];
                        }
                        break;
                    }
                }

                Yamg::append_cell(ref_copy[1], ref_copy[3], ref_copy[4], cid,
                                             -1, -1, -1, cavity_boundary[i]);
                // Need to determine unique orientation to split trapazoid.
                if(_lnn2unn[ref_copy[2]]>_lnn2unn[ref_copy[0]]) {
                    Yamg::append_cell(ref_copy[3], ref_copy[0], ref_copy[2], cid,
                                                 -1, -1, -1, cavity_boundary[i]);

                    Yamg::append_cell(ref_copy[4], ref_copy[3], ref_copy[2], cid,
                                                 -1, -1, -1, cavity_boundary[i]);
                } else {
                    Yamg::append_cell(ref_copy[3], ref_copy[0], ref_copy[4], cid,
                                                 -1, -1, -1, cavity_boundary[i]);

                    Yamg::append_cell(ref_copy[4], ref_copy[0], ref_copy[2], cid,
                                                 -1, -1, -1, cavity_boundary[i]);
                }
            } else if(cutcnt==3) {
                Yamg::append_cell(cavity[i*6+3], cavity[i*6+0], cavity[i*6+5], cid,
                                             -1, -1, -1, cavity_boundary[i]);

                Yamg::append_cell(cavity[i*6+1], cavity[i*6+3], cavity[i*6+4], cid,
                                             -1, -1, -1, cavity_boundary[i]);

                Yamg::append_cell(cavity[i*6+4], cavity[i*6+5], cavity[i*6+2], cid,
                                             -1, -1, -1, cavity_boundary[i]);

                Yamg::append_cell(cavity[i*6+4], cavity[i*6+3], cavity[i*6+5], cid,
                                             -1, -1, -1, cavity_boundary[i]);
            }
        }
    } else {
        // Mesh simple case.
        Yamg::append_cell(lnn[1], lnn[0], lnn[3], lnn[7],
                                     -1, _boundary[1], -1, _boundary[4]);
        Yamg::append_cell(lnn[7], lnn[4], lnn[5], lnn[0],
                                     _boundary[2], -1, -1, _boundary[5]);
        Yamg::append_cell(lnn[1], lnn[0], lnn[7], lnn[5],
                                     -1, _boundary[1], _boundary[2], -1);
        Yamg::append_cell(lnn[3], lnn[0], lnn[2], lnn[7],
                                     -1, _boundary[3], -1, _boundary[4]);
        Yamg::append_cell(lnn[6], lnn[4], lnn[7], lnn[0],
                                     -1, -1, _boundary[0], _boundary[5]);
        Yamg::append_cell(lnn[7], lnn[0], lnn[2], lnn[6],
                                     _boundary[0], _boundary[3], -1, -1);
    }
}

extern "C" {
    void create_verts (p4est_iter_volume_info_t * info, void *user_data)
    {
        double quad[3*8];
        Yamg::generate_quadcoords(info->p4est, info->treeid, info->quad, quad, NULL, NULL);

        // Insert these points
        for(int i=0; i<8; i++) {
            long long unn = coord_hash(quad[i*3], quad[i*3+1], quad[i*3+2]);
            if(unn2lnn.count(unn)==0) {
                int lnn = gcoords.size()/3;
                unn2lnn[unn] = lnn;
                lnn2unn.push_back(unn);

                for(int l=0; l<3; l++) {
                    gcoords.push_back(quad[i*3+l]);
                }
            }
        }
    }
}

Yamg::Yamg(int *argc, char ***argv)
{
    // Initialise MPI
    int mpiret = sc_MPI_Init (argc, argv);
    SC_CHECK_MPI (mpiret);

    mpicomm = sc_MPI_COMM_WORLD;

    // Initialise p4est
    sc_init(mpicomm, 1, 1, NULL, SC_LP_SILENT);
    p4est_init(NULL, SC_LP_SILENT);

    p4est = NULL;
    conn = NULL;
    ghost = NULL;

    debugging = false;
    verbose = false;

    MPI_Comm_size(mpicomm, &nranks);
    assert(nranks==1);

    MPI_Comm_rank(mpicomm, &rank);
    assert(rank==0);
}

Yamg::~Yamg() {}

void Yamg::enable_debugging()
{
    debugging = true;
}

void Yamg::enable_verbose()
{
    verbose = true;
}

void Yamg::finalise()
{
    if(verbose && rank==0)
        std::cout<<"INFO: void Yamg::finalise()\n";

    MPI_Barrier(mpicomm);

    // Destroy the p4est.
    if(p4est!=NULL)
        p4est_destroy (p4est);

    // Distroy connectivity.
    if(conn!=NULL)
        p4est_connectivity_destroy (conn);

    if(ghost!=NULL)
        p4est_ghost_destroy (ghost);

    /* Verify that allocations internal to p4est and sc do not leak memory.
     * This should be called if sc_init () has been called earlier. */
    sc_finalize ();

    /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
    int mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
}

void Yamg::init_domain()
{
    init_coord_hash();
    create_p4est();
}

// Init geometric hashing.
int Yamg::init_coord_hash()
{
    // long long lx = (long long)ceil((long double)1.1*(geom_bounds[1]-geom_bounds[0])/resolution);
    long long ly = (long long)ceil(1.1*(geom_bounds[3]-geom_bounds[2])/hash_resolution);
    long long lz = (long long)ceil(1.1*(geom_bounds[5]-geom_bounds[4])/hash_resolution);

    _hash_stride[0]=lz*ly;
    _hash_stride[1]=lz;

    hash_trusty = true;

    return 0;
}

void Yamg::create_p4est()
{
    // Determine a sensible coarse element size.
    double lx = geom_bounds[1] - geom_bounds[0];
    double ly = geom_bounds[3] - geom_bounds[2];
    double lz = geom_bounds[5] - geom_bounds[4];

    assert(lx>0 && ly>0 && lz>0);
    double base_size=std::min(lx, std::min(ly, lz));

    // Get dimensions.
    int ni=round(lx/base_size);
    int nj=round(ly/base_size);
    int nk=round(lz/base_size);
    double de[]= {lx/ni, ly/nj, lz/nk};

    // Create points.
    std::vector<double> pts;
    pts.reserve((ni+1)*(nj+1)*(nk+1));
    double x[3];
    for(int k=0; k<=nk; k++) {
        x[2] = geom_bounds[4] + de[2]*k;
        for(int j=0; j<=nj; j++) {
            x[1] = geom_bounds[2] + de[1]*j;
            for(int i=0; i<=ni; i++) {
                x[0] = geom_bounds[0] + de[0]*i;

                pts.insert(pts.end(), x, x+3);
            }
        }
    }

    // Create hexs.
    std::vector<int> hexs;
    int nij=(ni+1)*(nj+1);
    for(int k=0; k<nk; k++) {
        int k0=k*nij;
        int k1=(k+1)*nij;
        for(int j=0; j<nj; j++) {
            int j0=j*(ni+1);
            int j1=(j+1)*(ni+1);
            for(int i=0; i<ni; i++) {
                hexs.push_back(k0 + j0 + i);
                hexs.push_back(k0 + j0 + i+1);
                hexs.push_back(k0 + j1 + i);
                hexs.push_back(k0 + j1 + i+1);

                hexs.push_back(k1 + j0 + i);
                hexs.push_back(k1 + j0 + i+1);
                hexs.push_back(k1 + j1 + i);
                hexs.push_back(k1 + j1 + i+1);
            }
        }
    }

    // Create a forest.
    int num_vertices = pts.size()/3;
    int num_trees = hexs.size()/8;

    conn = p4est_connectivity_new(num_vertices, num_trees,
                                  0, 0, 0, 0);

    // Fill in verticies.
    for(int i=0; i<num_vertices*3; i++)
        conn->vertices[i] = pts[i];

    // Fill in elements.
    for(int i=0; i<num_trees; i++) {
        for(int j=0; j<8; j++) {
            conn->tree_to_vertex[i*8+j]   = hexs[i*8+j];
        }
    }

    // Create a forest that is not refined.
    p4est = p4est_new(mpicomm, conn, 0, NULL, NULL);

    // Fill tree_to_tree and tree_to_face to make sure we have a valid connectivity.
    for(int tree=0; tree<conn->num_trees; ++tree) {
        for(int face=0; face<P4EST_FACES; ++face) {
            conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
            conn->tree_to_face[P4EST_FACES * tree + face] = face;
        }
    }

    P4EST_ASSERT (p4est_connectivity_is_valid (conn));

    // Compute real tree_to_* fields and complete (edge and) corner fields.
    p4est_connectivity_complete (conn);
}

void Yamg::refine_p4est(p4est_refine_t refine_fn)
{
    double tic = MPI_Wtime();
    if(verbose && rank==0)
        std::cout<<"INFO: void Yamg::refine_p4est(...) :: ";

    p4est_refine (p4est, 1, refine_fn, NULL);
    p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

    if(verbose && rank==0)
        std::cout<<MPI_Wtime()-tic<<" seconds\n";
}

void Yamg::triangulate()
{
    double tic = MPI_Wtime();
    if(verbose && rank==0)
        std::cout<<"INFO: void Yamg::triangulate() :: ";

    p4est_iterate (p4est, /* the forest */
                   ghost,
                   NULL, /* the synchronized ghost data */
                   create_verts,
                   NULL,
                   NULL,
                   NULL);

    p4est_iterate (p4est, /* the forest */
                   ghost,
                   NULL, /* the synchronized ghost data */
                   mesh_quad,
                   NULL,
                   NULL,
                   NULL);

    unn2lnn.clear();
    int ntets = gcells.size()/4;
    for(int i=0; i<ntets; i++) {
        const int *n = Yamg::get_tet(i);
        for(int j=0; j<4; j++) {
            int lnn = n[j];
            long long unn = lnn2unn[lnn];
            unn2lnn[unn] = lnn;
        }
    }

    // Hereafter need to use the lnn2unn and unn2lnn luts.
    hash_trusty = false;

    if(verbose && rank==0)
        std::cout<<MPI_Wtime()-tic<<" seconds\n";
}

void Yamg::get_element_facet(int eid, int fid, int facet[3]) const
{
    assert(eid>=0);
    assert(eid<gcells.size()/4);

    const int *n = Yamg::get_tet(eid);

    switch(fid) {
    case 0:
        facet[0] = n[1];
        facet[1] = n[3];
        facet[2] = n[2];
        break;
    case 1:
        facet[0] = n[2];
        facet[1] = n[3];
        facet[2] = n[0];
        break;
    case 2:
        facet[0] = n[0];
        facet[1] = n[3];
        facet[2] = n[1];
        break;
    case 3:
        facet[0] = n[0];
        facet[1] = n[1];
        facet[2] = n[2];
        break;
    default:
        assert(false);
    }
}

#include <sstream>

int Yamg::get_number_points()
{
    return gcoords.size()/3;
}

int Yamg::get_number_tets()
{
    return gcells.size()/4;
}

void Yamg::write_vtu(std::string basename)
{
    double tic = MPI_Wtime();
    if(verbose)
        std::cout<<"INFO: void Yamg::write_vtu(std::string basename) :: ";

    assert(rank==0);

    int NNodes=gcoords.size()/3;

    std::string filename(basename);
    filename += ".vtu";

    NNodes=gcoords.size()/3;

    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> surface = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Set points
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetNumberOfPoints(NNodes);

    for(int i=0; i<NNodes; i++) {
        pts->SetPoint(i, Yamg::get_point(i));
    }
    ug->SetPoints(pts);
    surface->SetPoints(pts);

    int Ntets=gcells.size()/4;

    vtkSmartPointer<vtkDoubleArray> vtk_data = vtkSmartPointer<vtkDoubleArray>::New();
    vtk_data->SetNumberOfComponents(1);
    vtk_data->SetName("data");
    vtk_data->SetNumberOfTuples(Ntets);

    ug->Allocate(Ntets);
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    ids->SetNumberOfIds(4);
    for(int i=0; i<Ntets; i++) {
        assert(gcells[i*4]>=0);

        ids->SetId(0, gcells[i*4]);
        ids->SetId(1, gcells[i*4+1]);
        ids->SetId(2, gcells[i*4+3]);
        ids->SetId(3, gcells[i*4+2]);

        ug->InsertNextCell(VTK_TETRA, ids);

        vtk_data->SetTuple1(i, get_scalar_p0(i));
    }
    ug->GetCellData()->AddArray(vtk_data);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetCompressorTypeToZLib();
    writer->SetDataModeToBinary();

    writer->SetFileName(filename.c_str());

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(ug);
#else
    writer->SetInputData(ug);
#endif
    writer->Write();

    vtkSmartPointer<vtkIntArray> vtk_boundary = vtkSmartPointer<vtkIntArray>::New();
    vtk_boundary->SetNumberOfComponents(1);
    vtk_boundary->SetName("boundary");
    vtk_boundary->SetNumberOfTuples(Ntets*4);

    surface->Allocate(Ntets*4);
    vtkSmartPointer<vtkIdList> ids_tri = vtkSmartPointer<vtkIdList>::New();
    ids_tri->SetNumberOfIds(3);
    for(int i=0; i<Ntets; i++) {
        for(int j=0;j<4;j++) {
            for(int k=0;k<3;k++) {
                ids_tri->SetId(k, gcells[i*4+(j+1+k)%4]);
            }
            vtk_boundary->SetTuple1(i*4+j, boundary[i*4+j]);
            surface->InsertNextCell(VTK_TRIANGLE, ids_tri);
        }
    }
    surface->GetCellData()->AddArray(vtk_boundary);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> swriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    swriter->SetCompressorTypeToZLib();
    swriter->SetDataModeToBinary();

    swriter->SetFileName("surface.vtu");

#if VTK_MAJOR_VERSION <= 5
    swriter->SetInput(surface);
#else
    swriter->SetInputData(surface);
#endif
    swriter->Write();

    if(verbose)
        std::cout<<MPI_Wtime()-tic<<" seconds"<<std::endl;
}

void Yamg::write_vti(std::string basename) const
{
    std::string filename(basename);
    filename += ".vti";

    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
    image->SetDimensions(nx);
    image->SetOrigin(ox);
    image->SetSpacing(dx);
    image->AllocateScalars(VTK_FLOAT, 1);

    int ipos=0;
    for(int k=0;k<nx[2];k++) {
        for(int j=0;j<nx[1];j++) {
            for(int i=0;i<nx[0];i++) {
                image->SetScalarComponentFromFloat(i, j, k, 0, data[ipos]);
            }
        }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> image_writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    image_writer->SetFileName(filename.c_str());
    image_writer->SetInputData(image);
    image_writer->Write();
}

void Yamg::build_facets(std::vector<int> &facets, std::vector<int> &facet_ids)
{
    int NTetra = get_number_tets();

    for(int i=0; i<NTetra; i++) {
        for(int j=0; j<4; j++) {
            if(boundary[i*4+j]>0) {
                int facet[3];
                get_element_facet(i, j, facet);
                facets.insert(facets.end(), facet, facet+3);
                facet_ids.push_back(boundary[i*4+j]);
            }
        }
    }
}

void Yamg::write_gmsh(std::string basename)
{
    double tic = MPI_Wtime();
    if(verbose)
        std::cout<<"INFO: void Yamg::write_gmsh(std::string basename="<<basename<<") :: ";

    int NNodes = gcoords.size()/3;
    int NTetra = gcells.size()/4;

    std::vector<int> facets;
    std::vector<int> facet_ids;
    build_facets(facets, facet_ids);

    int NFacets = facet_ids.size();
    assert(NFacets==(int)facets.size()/3);

    ofstream file;
    file.open(std::string(basename+".msh").c_str());
    file<<"$MeshFormat"<<std::endl
        <<"2.2 0 8"<<std::endl
        <<"$EndMeshFormat"<<std::endl
        <<"$Nodes"<<std::endl
        <<NNodes<<std::endl;
    file<<std::setprecision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<NNodes; i++) {
        file<<i+1<<" "<<gcoords[i*3]<<" "<<gcoords[i*3+1]<<" "<<gcoords[i*3+2]<<std::endl;
    }
    file<<"$EndNodes"<<std::endl
        <<"$Elements"<<std::endl
        <<NTetra+NFacets<<std::endl;
    for(int i=0; i<NTetra; i++) {
        file<<i+1<<" 4 1 1 "<<gcells[i*4+1]+1<<" "<<gcells[i*4]+1<<" "<<gcells[i*4+2]+1<<" "<<gcells[i*4+3]+1<<std::endl;
    }
    for(int i=0; i<NFacets; i++) {
        file<<i+NTetra+1<<" 2 1 "<<facet_ids[i]<<" "<<facets[i*3]+1<<" "<<facets[i*3+1]+1<<" "<<facets[i*3+2]+1<<std::endl;
    }
    file<<"$EndElements"<<std::endl;

    file<<"$ElementData"<<std::endl
        <<"1"<<std::endl // number-of-string-tags
        <<"Vp"<<std::endl
        <<"0"<<std::endl
        <<"0"<<std::endl;
    for(int i=0; i<NTetra; i++) {
        file<<i+1<<" "<<Yamg::get_scalar_p0(i)<<std::endl;
    }
    file<<"$EndElementData"<<std::endl;

    file.close();

    if(verbose)
        std::cout<<MPI_Wtime()-tic<<" seconds\n";
}

void Yamg::sanity_facet(const int *facet) const
{
    if(!debugging)
        return;

    for(int i=0; i<3; i++) {
        assert(facet[i]>=0);
        assert(facet[i]<gcoords.size()/3);

        assert(facet[i]!=facet[(i+1)%3]);
        assert(edge_length(Yamg::get_point(facet[i]), Yamg::get_point(facet[(i+1)%3]))>std::numeric_limits<double>::epsilon());
    }
}

void Yamg::sanity_mesh(const char *file, int line) const
{
    if(!debugging)
        return;

    std::cout<<"Checking mesh sanity: ";
    if(file!=NULL) {
        std::cout<<file<<" "<<line;
    }
    std::cout<<std::endl;

    int NNodes = gcoords.size()/3;
    int NTetra = gcells.size()/4;

    for(int i=0; i<NTetra; i++) {
        const int *n = Yamg::get_tet(i);

        for(int j=0; j<4; j++) {
            assert(n[j]>=0);
            assert(n[j]<NNodes);
            for(int k=j+1; k<4; k++) {
                assert(n[j]!=n[k]);
            }
        }

        if(tet_volume(i)<0) {
            std::cerr<<"ERROR: found negative volume: "<<i<<", "<<tet_volume(i)<<std::endl;
            exit(-1);
        }
        assert(tet_volume(i)>-1.0e-6);
    }
}

double Yamg::total_volume() const{
    long double vol=0;
    int NTetra = gcells.size()/4;
    for(int i=0; i<NTetra; i++) {
        vol += tet_volume(i);
    }
    
    return vol;
}

double Yamg::total_area() const {
    long double area=0;
    int NTetra = gcells.size()/4;
    int facet[3];
    for(int i=0; i<NTetra; i++) {
        for(int j=0;j<4;j++) {
            get_element_facet(i, j, facet);

            if(boundary[i*4+j]>=0) {
                const double *x0 = Yamg::get_point(facet[0]);
                const double *x1 = Yamg::get_point(facet[1]);
                const double *x2 = Yamg::get_point(facet[2]);
                area += tri_area(x0, x1, x2);
            }
        }
    }

    return area;
}

void Yamg::dump_stats(std::string tag)
{
    int NTetra = gcells.size()/4;
    double vol = total_volume();
    double area = total_area();

    double dx = geom_bounds[1]-geom_bounds[0];
    double dy = geom_bounds[3]-geom_bounds[2];
    double dz = geom_bounds[5]-geom_bounds[4];

    double exact_volume = dx*dy*dz;
    double exact_area = 2*(dx*dy + dx*dz + dy*dz);

    std::cout<<"Total volume: "<<vol<<" (error="<<exact_volume-vol<<")"<<std::endl;
    std::cout<<"Total area: "<<area<<" (error="<<exact_area-area<<")"<<std::endl;
}
