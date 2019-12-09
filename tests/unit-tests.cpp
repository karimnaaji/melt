#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#define MELT_DEBUG
#define MELT_ASSERT(stmt) assert(stmt)
#define MELT_IMPLEMENTATION
#include "melt.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include <math.h>

#define FABS(x) ((float)fabs(x)) 
#define USE_EPSILON_TEST TRUE  
#define EPSILON 0.000001
#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; dest[1]=v1[1]-v2[1]; dest[2]=v1[2]-v2[2]; 
#define ADD(dest,v1,v2) dest[0]=v1[0]+v2[0]; dest[1]=v1[1]+v2[1]; dest[2]=v1[2]+v2[2];
#define MULT(dest,v,factor) dest[0]=factor*v[0]; dest[1]=factor*v[1]; dest[2]=factor*v[2];
#define SET(dest,src) dest[0]=src[0]; dest[1]=src[1]; dest[2]=src[2]; 
#define SORT(a,b) \
    if(a>b)       \
    {             \
        float c;  \
        c=a;      \
        a=b;      \
        b=c;      \
    }
#define ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
              isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
              isect1=VV0+(VV2-VV0)*D0/(D0-D2);
#define COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
    if(D0D1>0.0f)                                         \
    {                                                     \
    /* here we know that D0D2<=0.0 */                     \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
        ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);        \
    }                                                     \
    else if(D0D2>0.0f)                                    \
    {                                                     \
        /* here we know that d0d1<=0.0 */                 \
        ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);        \
    }                                                     \
    else if(D1*D2>0.0f || D0!=0.0f)                       \
    {                                                     \
        /* here we know that d0d1<=0.0 or that D0!=0.0 */ \
        ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);        \
    }                                                     \
    else if(D1!=0.0f)                                     \
    {                                                     \
        ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);        \
    }                                                     \
    else if(D2!=0.0f)                                     \
    {                                                     \
        ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);        \
    }                                                     \
    else                                                  \
    {                                                     \
        /* triangles are coplanar */                      \
        return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);    \
    }

// this edge to edge test is based on Franlin Antonio's gem:
// "Faster Line Segment Intersection", in Graphics Gems III,
// pp. 199-202
#define EDGE_EDGE_TEST(V0,U0,U1)                        \
    Bx=U0[i0]-U1[i0];                                   \
    By=U0[i1]-U1[i1];                                   \
    Cx=V0[i0]-U0[i0];                                   \
    Cy=V0[i1]-U0[i1];                                   \
    f=Ay*Bx-Ax*By;                                      \
    d=By*Cx-Bx*Cy;                                      \
    if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
    {                                                   \
        e=Ax*Cy-Ay*Cx;                                  \
        if(f>0)                                         \
        {                                               \
            if(e>=0 && e<=f) return 1;                  \
        }                                               \
        else                                            \
        {                                               \
            if(e<=0 && e>=f) return 1;                  \
        }                                               \
    }                                

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2)   \
{                                                \
    float Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
    Ax=V1[i0]-V0[i0];                            \
    Ay=V1[i1]-V0[i1];                            \
    /* test edge U0,U1 against V0,V1 */          \
    EDGE_EDGE_TEST(V0,U0,U1);                    \
    /* test edge U1,U2 against V0,V1 */          \
    EDGE_EDGE_TEST(V0,U1,U2);                    \
    /* test edge U2,U1 against V0,V1 */          \
    EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)             \
{                                             \
    float a,b,c,d0,d1,d2;                     \
    /* is T1 completly inside T2? */          \
    /* check if V0 is inside tri(U0,U1,U2) */ \
    a=U1[i1]-U0[i1];                          \
    b=-(U1[i0]-U0[i0]);                       \
    c=-a*U0[i0]-b*U0[i1];                     \
    d0=a*V0[i0]+b*V0[i1]+c;                   \
                                              \
    a=U2[i1]-U1[i1];                          \
    b=-(U2[i0]-U1[i0]);                       \
    c=-a*U1[i0]-b*U1[i1];                     \
    d1=a*V0[i0]+b*V0[i1]+c;                   \
                                              \
    a=U0[i1]-U2[i1];                          \
    b=-(U0[i0]-U2[i0]);                       \
    c=-a*U2[i0]-b*U2[i1];                     \
    d2=a*V0[i0]+b*V0[i1]+c;                   \
    if(d0*d1>0.0)                             \
    {                                         \
        if(d0*d2>0.0) return 1;               \
    }                                         \
}

int coplanar_tri_tri(float N[3],float V0[3],float V1[3],float V2[3], float U0[3],float U1[3],float U2[3])
{
    float A[3];
    short i0,i1;
    /* first project onto an axis-aligned plane, that maximizes the area */
    /* of the triangles, compute indices: i0,i1. */
    A[0]=fabs(N[0]);
    A[1]=fabs(N[1]);
    A[2]=fabs(N[2]);
    if(A[0]>A[1])
    {
        if(A[0]>A[2])  
        {
            i0=1;      /* A[0] is greatest */
            i1=2;
        }
        else
        {
            i0=0;      /* A[2] is greatest */
            i1=1;
        }
    }
    else   /* A[0]<=A[1] */
    {
        if(A[2]>A[1])
        {
            i0=0;      /* A[2] is greatest */
            i1=1;                                           
        }
        else
        {
            i0=0;      /* A[1] is greatest */
            i1=2;
        }
    }               
                
    /* test all edges of triangle 1 against the edges of triangle 2 */
    EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);
                
    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    POINT_IN_TRI(V0,U0,U1,U2);
    POINT_IN_TRI(U0,V0,V1,V2);

    return 0;
}

int tri_tri_intersect(float V0[3], float V1[3], float V2[3], float U0[3], float U1[3], float U2[3])
{
    float E1[3],E2[3];
    float N1[3],N2[3],d1,d2;
    float du0,du1,du2,dv0,dv1,dv2;
    float D[3];
    float isect1[2], isect2[2];
    float du0du1,du0du2,dv0dv1,dv0dv2;
    short index;
    float vp0,vp1,vp2;
    float up0,up1,up2;
    float b,c,max;

    /* compute plane equation of triangle(V0,V1,V2) */
    SUB(E1,V1,V0);
    SUB(E2,V2,V0);
    CROSS(N1,E1,E2);
    d1=-DOT(N1,V0);
    /* plane equation 1: N1.X+d1=0 */

    /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
    du0=DOT(N1,U0)+d1;
    du1=DOT(N1,U1)+d1;
    du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
    if(fabs(du0)<EPSILON) du0=0.0;
    if(fabs(du1)<EPSILON) du1=0.0;
    if(fabs(du2)<EPSILON) du2=0.0;
#endif
    du0du1=du0*du1;
    du0du2=du0*du2;

    if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
        return 0;                    /* no intersection occurs */

    /* compute plane of triangle (U0,U1,U2) */
    SUB(E1,U1,U0);
    SUB(E2,U2,U0);
    CROSS(N2,E1,E2);
    d2=-DOT(N2,U0);
    /* plane equation 2: N2.X+d2=0 */

    /* put V0,V1,V2 into plane equation 2 */
    dv0=DOT(N2,V0)+d2;
    dv1=DOT(N2,V1)+d2;
    dv2=DOT(N2,V2)+d2;

#if USE_EPSILON_TEST==TRUE
    if(fabs(dv0)<EPSILON) dv0=0.0;
    if(fabs(dv1)<EPSILON) dv1=0.0;
    if(fabs(dv2)<EPSILON) dv2=0.0;
#endif

    dv0dv1=dv0*dv1;
    dv0dv2=dv0*dv2;

    if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

    /* compute direction of intersection line */
    CROSS(D,N1,N2);

    /* compute and index to the largest component of D */
    max=fabs(D[0]);
    index=0;
    b=fabs(D[1]);
    c=fabs(D[2]);
    if(b>max) max=b,index=1;
    if(c>max) max=c,index=2;

    /* this is the simplified projection onto L*/
    vp0=V0[index];
    vp1=V1[index];
    vp2=V2[index];

    up0=U0[index];
    up1=U1[index];
    up2=U2[index];

    /* compute interval for triangle 1 */
    COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,isect1[0],isect1[1]);

    /* compute interval for triangle 2 */
    COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,isect2[0],isect2[1]);

    SORT(isect1[0],isect1[1]);
    SORT(isect2[0],isect2[1]);

    if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
    return 1;
}


static bool LoadModelMesh(const char* model_path, melt_params_t& melt_params)
{
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string error;
    bool obj_parsing_res = tinyobj::LoadObj(shapes, materials, error, model_path, NULL);

    if (!error.empty() || !obj_parsing_res) return false;
    
    melt_params.mesh.vertices.clear();
    melt_params.mesh.indices.clear();
    for (size_t i = 0; i < shapes.size(); i++) 
    {
        for (size_t f = 0; f < shapes[i].mesh.indices.size(); f++) 
            melt_params.mesh.indices.push_back(shapes[i].mesh.indices[f]);

        for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) 
        {
            melt_params.mesh.vertices.emplace_back();
            melt_params.mesh.vertices.back().x = shapes[i].mesh.positions[3 * v + 0];
            melt_params.mesh.vertices.back().y = shapes[i].mesh.positions[3 * v + 1];
            melt_params.mesh.vertices.back().z = shapes[i].mesh.positions[3 * v + 2];
        }
    }

    return true;
}

bool EnsureMeshExclusive(const melt_mesh_t& mesh0, const melt_mesh_t& mesh1)
{
    for (u32 i = 0; i < mesh0.indices.size(); i += 3)
    {
        float v0[3];
        float v1[3];
        float v2[3];

        memcpy(v0, &mesh0.vertices[mesh0.indices[i + 0]], sizeof(float) * 3);
        memcpy(v1, &mesh0.vertices[mesh0.indices[i + 1]], sizeof(float) * 3);
        memcpy(v2, &mesh0.vertices[mesh0.indices[i + 2]], sizeof(float) * 3);

        float e1[3], e2[3];
        float normal[3];

        SUB(e1, v1, v0);
        SUB(e2, v2, v0);

        CROSS(normal, e1, e2);

        for (u32 j = 0; j < mesh1.indices.size(); j += 3)
        {
            float u0[3];
            float u1[3];
            float u2[3];

            memcpy(u0, &mesh1.vertices[mesh1.indices[j + 0]], sizeof(float) * 3);
            memcpy(u1, &mesh1.vertices[mesh1.indices[j + 1]], sizeof(float) * 3);
            memcpy(u2, &mesh1.vertices[mesh1.indices[j + 2]], sizeof(float) * 3);

            float p0[3];
            // Point lies on plane
            SUB(p0, v0, u0);
            if (abs(DOT(normal, p0)) < EPSILON) continue;
            SUB(p0, v0, u1);
            if (abs(DOT(normal, p0)) < EPSILON) continue;
            SUB(p0, v0, u2);
            if (abs(DOT(normal, p0)) < EPSILON) continue;

            if (tri_tri_intersect(&v0[0], &v1[0], &v2[0], &u0[0], &u1[0], &u2[0]))
            {
                return false;
            }
        }
    }
    return true;
}

TEST_CASE("melt.bunny", "") 
{
    melt_params_t params;
    memset(&params.debug, 0, sizeof(melt_debug_params_t));
    params.voxel_size = 0.25f;
    params.fill_pct = 1.0f;
    params.box_type_flags = MELT_OCCLUDER_BOX_TYPE_REGULAR;

    melt_result_t result;

    REQUIRE(LoadModelMesh("models/bunny.obj", params));
    REQUIRE(melt_generate_occluder(params, result));
    REQUIRE(EnsureMeshExclusive(params.mesh, result.mesh));

    params.voxel_size = 0.05f;
    REQUIRE(!melt_generate_occluder(params, result));
}

TEST_CASE("melt.suzanne", "") 
{
    melt_params_t params;
    memset(&params.debug, 0, sizeof(melt_debug_params_t));
    params.voxel_size = 0.25f;
    params.fill_pct = 1.0f;
    params.box_type_flags = MELT_OCCLUDER_BOX_TYPE_REGULAR;

    melt_result_t result;

    REQUIRE(LoadModelMesh("models/suzanne.obj", params));
    REQUIRE(melt_generate_occluder(params, result));
    REQUIRE(EnsureMeshExclusive(params.mesh, result.mesh));

    params.voxel_size = 0.15f;
    REQUIRE(melt_generate_occluder(params, result));
    REQUIRE(EnsureMeshExclusive(params.mesh, result.mesh));
}

TEST_CASE("melt.cube", "") 
{
    melt_params_t params;
    memset(&params.debug, 0, sizeof(melt_debug_params_t));
    params.voxel_size = 0.25f;
    params.fill_pct = 1.0f;
    params.box_type_flags = MELT_OCCLUDER_BOX_TYPE_REGULAR;

    melt_result_t result;

    REQUIRE(LoadModelMesh("models/cube.obj", params));
    REQUIRE(melt_generate_occluder(params, result));
    REQUIRE(EnsureMeshExclusive(params.mesh, result.mesh));
}

TEST_CASE("melt.sphere", "") 
{
    melt_params_t params;
    memset(&params.debug, 0, sizeof(melt_debug_params_t));
    params.voxel_size = 0.25f;
    params.fill_pct = 1.0f;
    params.box_type_flags = MELT_OCCLUDER_BOX_TYPE_REGULAR;

    melt_result_t result;

    REQUIRE(LoadModelMesh("models/sphere.obj", params));
    REQUIRE(melt_generate_occluder(params, result));
    REQUIRE(EnsureMeshExclusive(params.mesh, result.mesh));
}

TEST_CASE("melt.teapot", "") 
{
    melt_params_t params;
    memset(&params.debug, 0, sizeof(melt_debug_params_t));
    params.voxel_size = 0.25f;
    params.fill_pct = 1.0f;
    params.box_type_flags = MELT_OCCLUDER_BOX_TYPE_REGULAR;

    melt_result_t result;

    REQUIRE(LoadModelMesh("models/teapot.obj", params));
    REQUIRE(!melt_generate_occluder(params, result));

    params.voxel_size = 0.5f;
    REQUIRE(melt_generate_occluder(params, result));
    REQUIRE(EnsureMeshExclusive(params.mesh, result.mesh));
}

TEST_CASE("melt.column", "") 
{
    melt_params_t params;
    memset(&params.debug, 0, sizeof(melt_debug_params_t));
    params.voxel_size = 0.25f;
    params.fill_pct = 1.0f;
    params.box_type_flags = MELT_OCCLUDER_BOX_TYPE_REGULAR;

    melt_result_t result;

    REQUIRE(LoadModelMesh("models/column.obj", params));
    REQUIRE(melt_generate_occluder(params, result));
    REQUIRE(EnsureMeshExclusive(params.mesh, result.mesh));  
}
