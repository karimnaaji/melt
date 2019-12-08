//
// The MIT License (MIT)
//
// Copyright (c) 2019 Karim Naaji, karim.naaji@gmail.com
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#ifndef MELT_H
#define MELT_H

#include <vector>
#include <stdint.h>

#ifndef MELT_ASSERT
#define MELT_ASSERT(stmt) (void)(stmt)
#endif
#ifndef MELT_PROFILE_BEGIN
#define MELT_PROFILE_BEGIN()
#endif
#ifndef MELT_PROFILE_END
#define MELT_PROFILE_END()
#endif

typedef struct {float x, y, z;} vec3_t;

struct MeltMesh
{
    std::vector<vec3_t> vertices;
    std::vector<unsigned short> indices;
};

enum MeltOccluderBoxType
{
    MeltOccluderBoxTypeNone      = 0,
    MeltOccluderBoxTypeDiagonals = 1 << 0,
    MeltOccluderBoxTypeTop       = 1 << 1,
    MeltOccluderBoxTypeBottom    = 1 << 2,
    MeltOccluderBoxTypeSides     = 1 << 3,
    MeltOccluderBoxTypeRegular   = MeltOccluderBoxTypeSides | MeltOccluderBoxTypeTop | MeltOccluderBoxTypeBottom,
};

typedef int MeltOccluderBoxTypeFlags;

enum MeltDebugType
{
    MeltDebugTypeShowInner          = 1 << 0,
    MeltDebugTypeShowExtent         = 1 << 1,
    MeltDebugTypeShowResult         = 1 << 2,
    MeltDebugTypeShowOuter          = 1 << 3,
    MeltDebugTypeShowMinDistance    = 1 << 4,
    MeltDebugTypeShowSliceSelection = 1 << 5,
};

typedef int MeltDebugTypeFlags;

struct MeltDebugParams
{    
    uint32_t _start_canary;
    MeltDebugTypeFlags flags;

    int sliceIndexX;
    int sliceIndexY;
    int sliceIndexZ;
    int voxelX;
    int voxelY;
    int voxelZ;
    int extentIndex;
    int extentMaxStep;

    float voxelScale;

    uint32_t _end_canary;
};

struct MeltParams
{
    MeltMesh mesh;
    MeltOccluderBoxTypeFlags boxTypeFlags;
    MeltDebugParams debug;
    float voxelSize;
    float fillPercentage;
};

struct MeltResult
{
    MeltMesh mesh;
    MeltMesh debugMesh;
};

bool MeltGenerateOccluder(const MeltParams& params, MeltResult& result);

#endif // MELT_H

#ifdef MELT_IMPLEMENTATION

#include <math.h>
#include <float.h>

typedef unsigned int    u32;
typedef unsigned short  u16;
typedef unsigned char   u8;
typedef signed long     s64;
typedef signed short    s16;
typedef signed int      s32;
typedef signed char     s8;
typedef float           f32;

typedef struct {f32 x, y;} vec2_t;

typedef struct {s32 x, y, z;} svec3_t;

typedef struct {u32 x, y;} uvec2_t;
typedef struct {u32 x, y, z;} uvec3_t;

vec3_t vec3_new(float x, float y, float z)
{
    vec3_t v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

uvec3_t uvec3_new(vec3_t v)
{ 
    uvec3_t out;
    out.x = (u32)v.x;
    out.y = (u32)v.y;
    out.z = (u32)v.z;
    return out; 
}


vec3_t vec3_new(uvec3_t v)
{
    vec3_t out;
    out.x = (float)v.x;
    out.y = (float)v.y;
    out.z = (float)v.z;
    return out;   
}


vec3_t vec3_mul(vec3_t v, float factor) {
    return vec3_new(v.x * factor, v.y * factor, v.z * factor);
}

vec3_t vec3_mul(vec3_t a, vec3_t b) {
    return vec3_new(a.x * b.x, a.y * b.y, a.z * b.z);
}


vec3_t vec3_div(vec3_t v, float divisor) {
    return vec3_mul(v, 1.0f / divisor);
}

uvec3_t uvec3_new(float x, float y, float z)
{
    uvec3_t v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;    
}

vec3_t vec3_sub(vec3_t a, vec3_t b)
{
    return vec3_new(a.x - b.x, a.y - b.y, a.z - b.z);
}

vec3_t vec3_add(vec3_t a, vec3_t b)
{
    return vec3_new(a.x + b.x, a.y + b.y, a.z + b.z);
}

vec3_t vec3_abs(vec3_t v)
{
    return vec3_new(fabsf(v.x), fabsf(v.y), fabsf(v.z));
}

vec3_t vec3_cross(vec3_t a, vec3_t b)
{
    float x = a.y * b.z - a.z * b.y;
    float y = a.z * b.x - a.x * b.z;
    float z = a.x * b.y - a.y * b.x;
    return vec3_new(x, y, z);
}

f32 vec3_dot(vec3_t a, vec3_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

uvec2_t uvec2_new(u32 x, u32 y)
{
    uvec2_t v;
    v.x = x;
    v.y = y;
    return v;
}

svec3_t svec3_new(s32 x, s32 y, s32 z)
{
    svec3_t v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

int svec3_equals(svec3_t a, svec3_t b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

int uvec3_equals(uvec3_t a, uvec3_t b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;   
}

float float_min(float a, float b) {
    return a < b ? a : b;
}

float float_max(float a, float b) {
    return a > b ? a : b;
}

u32 u32_min(u32 a, u32 b)
{
    return a > b ? a : b;
}

vec3_t vec3_min(vec3_t a, vec3_t b) {
    float x = float_min(a.x, b.x);
    float y = float_min(a.y, b.y);
    float z = float_min(a.z, b.z);
    return vec3_new(x, y, z);
}

vec3_t vec3_max(vec3_t a, vec3_t b) {
    float x = float_max(a.x, b.x);
    float y = float_max(a.y, b.y);
    float z = float_max(a.z, b.z);
    return vec3_new(x, y, z);
}

typedef uvec3_t Color3u8;

const Color3u8 Blueviolet     = { 138,  43, 226 };
const Color3u8 Darkseagreen   = { 143, 188, 143 };
const Color3u8 Lightslategray = { 119, 136, 153 };
const Color3u8 Powderblue     = { 176, 224, 230 };
const Color3u8 Springgreen    = {   0, 255, 127 };
const Color3u8 Steelblue      = {  70, 130, 180 };
const Color3u8 Teal           = {   0, 128, 128 };
const Color3u8 Whitesmoke     = { 245, 245, 245 };
const Color3u8 Floralwhite    = { 255, 250, 240 };
const Color3u8 Lightpink      = { 255, 182, 193 };

const Color3u8 Colors[] =
{
    Whitesmoke,
    Steelblue,
    Springgreen,
    Teal,
    Lightpink,
    Powderblue,
    Lightslategray,
    Darkseagreen,
    Floralwhite,
};

const u16 VoxelCubeIndices[36] =
{
    0, 1, 2,
    0, 2, 3,
    3, 2, 6,
    3, 6, 7,
    0, 7, 4,
    0, 3, 7,
    4, 7, 5,
    7, 6, 5,
    0, 4, 5,
    0, 5, 1,
    1, 5, 6,
    1, 6, 2,
};

const u16 VoxelCubeIndicesSides[24] =
{
    0, 1, 2,
    0, 2, 3,
    3, 2, 6,
    3, 6, 7,
    4, 7, 5,
    7, 6, 5,
    0, 4, 5,
    0, 5, 1,
};

const u16 VoxelCubeIndicesDiagonals[12] =
{
    0, 1, 6,
    0, 6, 7,
    4, 5, 2,
    4, 2, 3,
};

const u16 VoxelCubeIndicesBottom[6] =
{
    1, 5, 6,
    1, 6, 2,
};

const u16 VoxelCubeIndicesTop[6] =
{
    0, 7, 4,
    0, 3, 7,
};

const vec3_t VoxelCubeVertices[8] =
{
    {-1.0f,  1.0f,  1.0f},
    {-1.0f, -1.0f,  1.0f},
    { 1.0f, -1.0f,  1.0f},
    { 1.0f,  1.0f,  1.0f},
    {-1.0f,  1.0f, -1.0f},
    {-1.0f, -1.0f, -1.0f},
    { 1.0f, -1.0f, -1.0f},
    { 1.0f,  1.0f, -1.0f},
};

struct AABB
{
    vec3_t min;
    vec3_t max;
};

vec3_t aabb_center(AABB aabb)
{
    return vec3_mul(vec3_add(aabb.min, aabb.max), 0.5f);
}

struct Voxel
{
    AABB aabb;
    uvec3_t position;
};

struct Triangle
{
    vec3_t v0;
    vec3_t v1;
    vec3_t v2;
};

struct Plane
{
    vec3_t normal;
    f32 distance;
};

struct MinDistance
{
    MinDistance() {}
    svec3_t dist;
    union
    {
        struct
        {
            u32 x;
            u32 y;
            u32 z;
        };
        uvec3_t position;
    };
};

struct VoxelStatus
{
    u8 visibility : 6;
    u8 clipped : 1;
    u8 inner : 1;
};

enum Visibility
{
    VisibilityNull   =      0,
    VisibilityPlusX  = 1 << 0,
    VisibilityMinusX = 1 << 1,
    VisibilityPlusY  = 1 << 2,
    VisibilityMinusY = 1 << 3,
    VisibilityPlusZ  = 1 << 4,
    VisibilityMinusZ = 1 << 5,
    VisibilityAll    =   0x3f,
};

struct MaxExtent
{
    uvec3_t position;
    uvec3_t extent;
    u32 volume;
};

using VoxelIndices = std::vector<s32>;
using VoxelSet = std::vector<Voxel>;
using MinDistanceField = std::vector<MinDistance>;
using VoxelField =  std::vector<VoxelStatus>;
using MaxExtents =  std::vector<MaxExtent>;

struct Context
{
    uvec3_t dimension;
    u32 size;
};

struct VoxelSetPlanes
{   
    std::vector<VoxelSet> x;
    std::vector<VoxelSet> y;
    std::vector<VoxelSet> z;
};

#define ARRAY_LENGTH(array) ((int)(sizeof(array) / sizeof(*array)))

static bool AABBCollidesPlane(const Plane plane, const vec3_t& half_aabb_dim)
{
    vec3_t vmin;
    vec3_t vmax;

    if (plane.normal.x > 0.0f)
    {
        vmin.x = -half_aabb_dim.x;
        vmax.x =  half_aabb_dim.x;
    }
    else
    {
        vmin.x =  half_aabb_dim.x;
        vmax.x = -half_aabb_dim.x;
    }

    if (plane.normal.y > 0.0f)
    {
        vmin.y = -half_aabb_dim.y;
        vmax.y =  half_aabb_dim.y;
    }
    else
    {
        vmin.y =  half_aabb_dim.y;
        vmax.y = -half_aabb_dim.y;
    }

    if (plane.normal.z > 0.0f)
    {
        vmin.z = -half_aabb_dim.z;
        vmax.z =  half_aabb_dim.z;
    }
    else
    {
        vmin.z =  half_aabb_dim.z;
        vmax.z = -half_aabb_dim.z;
    }

    if (vec3_dot(plane.normal, vmin) + plane.distance > 0.0f)
        return false;

    if (vec3_dot(plane.normal, vmax) + plane.distance >= 0.0f)
        return true;

    return false;
}

#define AXISTEST_X01(a, b, fa, fb)                 \
p0 = a * v0.y - b * v0.z;                          \
p2 = a * v2.y - b * v2.z;                          \
if (p0 < p2) {                                     \
    min = p0; max = p2;                            \
} else {                                           \
    min = p2; max = p0;                            \
}                                                  \
rad = fa * half_aabb_dim.y + fb * half_aabb_dim.z; \
if (min > rad || max < -rad) {                     \
    return false;                                  \
}                                                  \

#define AXISTEST_X2(a, b, fa, fb)                  \
p0 = a * v0.y - b * v0.z;                          \
p1 = a * v1.y - b * v1.z;                          \
if (p0 < p1) {                                     \
    min = p0; max = p1;                            \
} else {                                           \
    min = p1; max = p0;                            \
}                                                  \
rad = fa * half_aabb_dim.y + fb * half_aabb_dim.z; \
if (min > rad || max < -rad) {                     \
    return false;                                  \
}                                                  \

#define AXISTEST_Y02(a, b, fa, fb)                 \
p0 = -a * v0.x + b * v0.z;                         \
p2 = -a * v2.x + b * v2.z;                         \
if (p0 < p2) {                                     \
    min = p0; max = p2;                            \
} else {                                           \
    min = p2; max = p0;                            \
}                                                  \
rad = fa * half_aabb_dim.x + fb * half_aabb_dim.z; \
if (min > rad || max < -rad) {                     \
    return false;                                  \
}                                                  \

#define AXISTEST_Y1(a, b, fa, fb)                  \
p0 = -a * v0.x + b * v0.z;                         \
p1 = -a * v1.x + b * v1.z;                         \
if (p0 < p1) {                                     \
    min = p0; max = p1;                            \
} else {                                           \
    min = p1; max = p0;                            \
}                                                  \
rad = fa * half_aabb_dim.x + fb * half_aabb_dim.z; \
if (min > rad || max < -rad) {                     \
    return false;                                  \
}

#define AXISTEST_Z12(a, b, fa, fb)                 \
p1 = a * v1.x - b * v1.y;                          \
p2 = a * v2.x - b * v2.y;                          \
if (p2 < p1) {                                     \
    min = p2; max = p1;                            \
} else {                                           \
    min = p1; max = p2;                            \
}                                                  \
rad = fa * half_aabb_dim.x + fb * half_aabb_dim.y; \
if (min > rad || max < -rad) {                     \
    return false;                                  \
}

#define AXISTEST_Z0(a, b, fa, fb)                  \
p0 = a * v0.x - b * v0.y;                          \
p1 = a * v1.x - b * v1.y;                          \
if (p0 < p1) {                                     \
    min = p0; max = p1;                            \
} else {                                           \
    min = p1; max = p0;                            \
}                                                  \
rad = fa * half_aabb_dim.x + fb * half_aabb_dim.y; \
if (min > rad || max < -rad) {                     \
    return false;                                  \
}

#define FINDMINMAX(x0, x1, x2, min, max)           \
min = max = x0;                                    \
if (x1 < min) min = x1;                            \
if (x1 > max) max = x1;                            \
if (x2 < min) min = x2;                            \
if (x2 > max) max = x2;

static bool AABBIntersectsTriangle(const Triangle& triangle, const vec3_t& aabb_center, const vec3_t& half_aabb_dim)
{
    vec3_t v0, v1, v2;
    vec3_t e0, e1, e2;
    vec3_t edge_abs;
    f32 min, max;
    f32 p0, p1, p2;
    f32 rad;

    v0 = vec3_sub(triangle.v0, aabb_center);
    v1 = vec3_sub(triangle.v1, aabb_center);
    v2 = vec3_sub(triangle.v2, aabb_center);

    e0 = vec3_sub(v1, v0);
    e1 = vec3_sub(v2, v1);
    e2 = vec3_sub(v0, v2);

    edge_abs = vec3_abs(e0);
    AXISTEST_X01(e0.z, e0.y, edge_abs.z, edge_abs.y);
    AXISTEST_Y02(e0.z, e0.x, edge_abs.z, edge_abs.x);
    AXISTEST_Z12(e0.y, e0.x, edge_abs.y, edge_abs.x);

    edge_abs = vec3_abs(e1);
    AXISTEST_X01(e1.z, e1.y, edge_abs.z, edge_abs.y);
    AXISTEST_Y02(e1.z, e1.x, edge_abs.z, edge_abs.x);
    AXISTEST_Z0 (e1.y, e1.x, edge_abs.y, edge_abs.x);

    edge_abs = vec3_abs(e2);
    AXISTEST_X2 (e2.z, e2.y, edge_abs.z, edge_abs.y);
    AXISTEST_Y1 (e2.z, e2.x, edge_abs.z, edge_abs.x);
    AXISTEST_Z12(e2.y, e2.x, edge_abs.y, edge_abs.x);

    FINDMINMAX(v0.x, v1.x, v2.x, min, max);
    if (min > half_aabb_dim.x || max < -half_aabb_dim.x)
        return false;

    FINDMINMAX(v0.y, v1.y, v2.y, min, max);
    if (min > half_aabb_dim.y || max < -half_aabb_dim.y)
        return false;

    FINDMINMAX(v0.z, v1.z, v2.z, min, max);
    if (min > half_aabb_dim.z || max < -half_aabb_dim.z)
        return false;

    Plane plane;
    plane.normal = vec3_cross(e0, e1);
    plane.distance = -vec3_dot(plane.normal, v0);

    if (!AABBCollidesPlane(plane, half_aabb_dim))
        return false;

    return true;
}

static inline u32 Flatten3d(const uvec3_t index, const uvec3_t dimension)
{
    u32 out_index = index.x + dimension.x * index.y + dimension.x * dimension.y * index.z;
    MELT_ASSERT(out_index < dimension.x * dimension.y * dimension.z);
    return out_index;
}
    
static inline u32 Flatten2d(const uvec2_t index, const uvec2_t dimension)
{
    u32 out_index = index.x + dimension.x * index.y;
    MELT_ASSERT(out_index < dimension.x * dimension.y);
    return out_index;
}

static inline uvec3_t UnFlatten3d(u32 position, const uvec3_t& dimension)
{
    uvec3_t out_index;

    u32 dim_xy = dimension.x * dimension.y;
    out_index.z = position / dim_xy;
    position -= out_index.z * dim_xy;
    out_index.y = position / dimension.x;
    out_index.x = position % dimension.x;

    MELT_ASSERT(out_index.x < dimension.x);
    MELT_ASSERT(out_index.y < dimension.y);
    MELT_ASSERT(out_index.z < dimension.z);

    return out_index;
}

static vec3_t MapToVoxelMaxBound(const vec3_t& position, f32 voxel_size)
{
    auto MapFunc = [&](f32 value)
    {
        f32 sign = value < 0.0f ? -1.0f : 1.0f;
        f32 result = value + sign * voxel_size * 0.5f;
        return ceil(result / voxel_size) * voxel_size;
    };

    return vec3_new(
        MapFunc(position.x),
        MapFunc(position.y),
        MapFunc(position.z));
}

static vec3_t MapToVoxelMinBound(const vec3_t& position, f32 voxel_size)
{
    auto MapFunc = [&](f32 value)
    {
        f32 sign = value < 0.0f ? -1.0f : 1.0f;
        f32 result = value + sign * voxel_size * 0.5f;
        return floorf(result / voxel_size) * voxel_size;
    };

    return vec3_new(
        MapFunc(position.x),
        MapFunc(position.y),
        MapFunc(position.z));
}

static AABB GenerateAABB(const Triangle& triangle)
{
    AABB aabb;

    aabb.min = vec3_new( FLT_MAX,  FLT_MAX,  FLT_MAX);
    aabb.max = vec3_new(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    aabb.min = vec3_min(aabb.min, triangle.v0);
    aabb.max = vec3_max(aabb.max, triangle.v0);

    aabb.min = vec3_min(aabb.min, triangle.v1);
    aabb.max = vec3_max(aabb.max, triangle.v1);

    aabb.min = vec3_min(aabb.min, triangle.v2);
    aabb.max = vec3_max(aabb.max, triangle.v2);

    return aabb;
}

static AABB GenerateAABB(const MeltMesh& mesh)
{
    AABB aabb;

    aabb.min = vec3_new( FLT_MAX,  FLT_MAX,  FLT_MAX);
    aabb.max = vec3_new(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (u32 i = 0; i < mesh.indices.size(); ++i)
    {
        aabb.min = vec3_min(aabb.min, mesh.vertices[mesh.indices[i]]);
        aabb.max = vec3_max(aabb.max, mesh.vertices[mesh.indices[i]]);
    }

    return aabb;
}

static void GeneratePerPlaneVoxelSet(const Context& context, const VoxelSet& voxel_set, const VoxelIndices& voxel_indices, VoxelSetPlanes& out_voxel_set_planes)
{
    MELT_PROFILE_BEGIN();

    out_voxel_set_planes.x.resize(context.dimension.y * context.dimension.z);
    out_voxel_set_planes.y.resize(context.dimension.x * context.dimension.z);
    out_voxel_set_planes.z.resize(context.dimension.x * context.dimension.y);
    
    uvec2_t dim_yz = uvec2_new(context.dimension.y, context.dimension.z);
    uvec2_t dim_xz = uvec2_new(context.dimension.x, context.dimension.z);
    uvec2_t dim_xy = uvec2_new(context.dimension.x, context.dimension.y);

    for (u32 x = 0; x < context.dimension.x; ++x)
    {
        for (u32 y = 0; y < context.dimension.y; ++y)
        {
            for (u32 z = 0; z < context.dimension.z; ++z)
            {
                uvec3_t position = uvec3_new(x, y, z);
                s32 voxel_index = voxel_indices[Flatten3d(position, context.dimension)];
                if (voxel_index != -1)
                {
                    u32 index_yz = Flatten2d(uvec2_new(y, z), dim_yz);
                    u32 index_xz = Flatten2d(uvec2_new(x, z), dim_xz);
                    u32 index_xy = Flatten2d(uvec2_new(x, y), dim_xy);

                    VoxelSet& voxels_x_planes = out_voxel_set_planes.x[index_yz];
                    VoxelSet& voxels_y_planes = out_voxel_set_planes.y[index_xz];
                    VoxelSet& voxels_z_planes = out_voxel_set_planes.z[index_xy];

                    voxels_x_planes.push_back(voxel_set[voxel_index]);
                    voxels_y_planes.push_back(voxel_set[voxel_index]);
                    voxels_z_planes.push_back(voxel_set[voxel_index]);
                }
            }
        }
    }

    MELT_PROFILE_END();
}

static void GetField(Context& context, const VoxelSetPlanes& voxel_set_planes, const VoxelIndices& voxel_indices, u32 x, u32 y, u32 z, MinDistance& out_min_distance, VoxelStatus& out_status)
{
    const svec3_t InfiniteDistance = svec3_new(INT_MAX, INT_MAX, INT_MAX);
    const svec3_t NullDistance = svec3_new(0, 0, 0);

    out_min_distance.dist = InfiniteDistance;
    out_min_distance.x = x;
    out_min_distance.y = y;
    out_min_distance.z = z;

    out_status.visibility = VisibilityNull;
    out_status.clipped = false;
    out_status.inner = false;

    uvec2_t dim_yz = uvec2_new(context.dimension.y, context.dimension.z);
    u32 index_yz = Flatten2d(uvec2_new(y, z), dim_yz);
    const auto& voxels_x_plane = voxel_set_planes.x[index_yz];
    for (const auto& voxel : voxels_x_plane)
    {
        s32 distance = voxel.position.x - x;
        if (distance > 0)
        {
            out_status.visibility |= VisibilityPlusX;
            out_min_distance.dist.x = u32_min(out_min_distance.dist.x, distance);
        }
        else if (distance < 0)
            out_status.visibility |= VisibilityMinusX;
        else
            out_min_distance.dist.x = 0;
    }
    uvec2_t dim_xz = uvec2_new(context.dimension.x, context.dimension.z);
    u32 index_xz = Flatten2d(uvec2_new(x, z), dim_xz);
    const auto& voxels_y_plane = voxel_set_planes.y[index_xz];
    for (const auto& voxel : voxels_y_plane)
    {
        s32 distance = voxel.position.y - y;
        if (distance > 0)
        {
            out_status.visibility |= VisibilityPlusY;
            out_min_distance.dist.y = u32_min(out_min_distance.dist.y, distance);
        }
        else if (distance < 0)
            out_status.visibility |= VisibilityMinusY;
        else
            out_min_distance.dist.y = 0;
    }
    uvec2_t dim_xy = uvec2_new(context.dimension.x, context.dimension.y);
    u32 index_xy = Flatten2d(uvec2_new(x, y), dim_xy);
    const auto& voxels_z_plane = voxel_set_planes.z[index_xy];
    for (const auto& voxel : voxels_z_plane)
    {
        s32 distance = voxel.position.z - z;
        if (distance > 0)
        {
            out_status.visibility |= VisibilityPlusZ;
            out_min_distance.dist.z = u32_min(out_min_distance.dist.z, distance);
        }
        else if (distance < 0)
            out_status.visibility |= VisibilityMinusZ;
        else
            out_min_distance.dist.z = 0;
    }
    if (out_status.visibility == VisibilityAll)
    {
        if (!svec3_equals(out_min_distance.dist, InfiniteDistance) &&
            !svec3_equals(out_min_distance.dist, NullDistance))
        {
            out_status.inner = true;
        }
    }
}

static void GenerateFields(Context& context, const VoxelSetPlanes& voxel_set_planes, const VoxelIndices& voxel_indices, VoxelField& out_voxel_field, MinDistanceField& out_distance_field)
{
    MELT_PROFILE_BEGIN();

    for (u32 i = 0; i < context.size; ++i)
    {
        MinDistance& min_distance = out_distance_field[i];
        VoxelStatus& voxel_status = out_voxel_field[i];
        const uvec3_t position = UnFlatten3d(i, context.dimension);

        GetField(context, voxel_set_planes, voxel_indices, position.x, position.y, position.z, min_distance, voxel_status);
    }

    MELT_PROFILE_END();
}

static inline bool InnerVoxel(VoxelStatus voxel_status)
{
    return voxel_status.inner && !voxel_status.clipped;
}

static uvec3_t GetMaxAABBExtent(const Context& context, const MinDistanceField& distance_field, const VoxelField& voxel_field, const MinDistance& min_distance)
{
    MELT_PROFILE_BEGIN();

    std::vector<uvec2_t> max_aabb_extents;

    for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
    {
        uvec3_t z_slice_position = uvec3_new(min_distance.x, min_distance.y, z);
        u32 z_slice_index = Flatten3d(z_slice_position, context.dimension);

        MELT_ASSERT(voxel_field[z_slice_index].inner);

        if (voxel_field[z_slice_index].clipped)
            continue;

        const MinDistance& sample_min_distance = distance_field[z_slice_index];

        uvec2_t max_extent = uvec2_new(sample_min_distance.dist.x, sample_min_distance.dist.y);

        u32 x = sample_min_distance.x + 1;
        u32 y = sample_min_distance.y + 1;
        u32 i = 1;
        while (x < sample_min_distance.x + sample_min_distance.dist.x &&
               y < sample_min_distance.y + sample_min_distance.dist.y)
        {
            const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
            if (InnerVoxel(voxel_field[index]))
            {
                const MinDistance& distance = distance_field[index];
                max_extent.x = u32_min(distance.dist.x + i, max_extent.x);
                max_extent.y = u32_min(distance.dist.y + i, max_extent.y);
            }
            else
            {
                max_extent.x = i;
                max_extent.y = i;
                break;
            }
            ++x;
            ++y;
            ++i;
        }

        max_aabb_extents.push_back(max_extent);
    }

    uvec2_t min_extent = uvec2_new(UINT_MAX, UINT_MAX);

    u32 z_slice = 1;
    u32 max_volume = 0;

    MELT_ASSERT(max_aabb_extents.size() > 0);

    for (const uvec2_t& extent : max_aabb_extents)
    {
        min_extent.x = u32_min(extent.x, min_extent.x);
        min_extent.y = u32_min(extent.y, min_extent.y);

        const u32 volume = min_extent.x * min_extent.y * z_slice;
        if (volume > max_volume)
            max_volume = volume;

        ++z_slice;
    }

    MELT_ASSERT(max_volume > 0);
    MELT_ASSERT(min_extent.x > 0);
    MELT_ASSERT(min_extent.y > 0);
    MELT_ASSERT(z_slice > 1);
    MELT_PROFILE_END();

    return uvec3_new(min_extent.x, min_extent.y, z_slice - 1);
}

static void ClipVoxelField(const Context& context, const uvec3_t start_position, const uvec3_t extent, VoxelField& voxel_field)
{
    MELT_PROFILE_BEGIN();

    for (u32 x = start_position.x; x < start_position.x + extent.x; ++x)
    {
        for (u32 y = start_position.y; y < start_position.y + extent.y; ++y)
        {
            for (u32 z = start_position.z; z < start_position.z + extent.z; ++z)
            {
                u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
                MELT_ASSERT(!voxel_field[index].clipped && "Clipping already clipped voxel field index");
                voxel_field[index].clipped = true;
            }
        }
    }

    MELT_PROFILE_END();
}
    
static bool WaterTightMesh(const Context& context, const VoxelField& voxel_field, const MinDistanceField& distance_field)
{
    for (const auto& min_distance : distance_field)
    {
        if (!InnerVoxel(voxel_field[Flatten3d(min_distance.position, context.dimension)]))
            continue;
        
        for (u32 x = min_distance.x; x < min_distance.x + min_distance.dist.x; ++x)
        {
            const u32 y = min_distance.y;
            const u32 z = min_distance.z;
            const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
            if (!InnerVoxel(voxel_field[index]))
            {
                return false;
            }
        }
        for (u32 y = min_distance.y; y < min_distance.y + min_distance.dist.y; ++y)
        {
            const u32 x = min_distance.x;
            const u32 z = min_distance.z;
            const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
            if (!InnerVoxel(voxel_field[index]))
            {
                return false;
            }
        }
        for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
        {
            const u32 x = min_distance.x;
            const u32 y = min_distance.y;
            const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
            if (!InnerVoxel(voxel_field[index]))
            {
                return false;
            }
        }
    }
    
    return true;
}

static void _Debug_ValidateMinDistanceField(const Context& context, const VoxelSet& outer_voxels, const VoxelField& voxel_field, const MinDistanceField& distance_field)
{
#if defined(MELT_DEBUG) && defined(MELT_ASSERT)
    for (const auto& min_distance : distance_field)
    {
        if (!InnerVoxel(voxel_field[Flatten3d(min_distance.position, context.dimension)]))
            continue;

        for (u32 x = min_distance.x; x < min_distance.x + min_distance.dist.x; ++x)
        {
            const u32 y = min_distance.y;
            const u32 z = min_distance.z;
            const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
            for (const Voxel& voxel : outer_voxels)
                MELT_ASSERT(!uvec3_equals(voxel.position, uvec3_new(x, y, z)));
            MELT_ASSERT(InnerVoxel(voxel_field[index]));
        }
        for (u32 y = min_distance.y; y < min_distance.y + min_distance.dist.y; ++y)
        {
            const u32 x = min_distance.x;
            const u32 z = min_distance.z;
            const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
            for (const Voxel& voxel : outer_voxels)
                MELT_ASSERT(!uvec3_equals(voxel.position, uvec3_new(x, y, z)));
            MELT_ASSERT(InnerVoxel(voxel_field[index]));
        }
        for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
        {
            const u32 x = min_distance.x;
            const u32 y = min_distance.y;
            const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
            for (const Voxel& voxel : outer_voxels)
                MELT_ASSERT(!uvec3_equals(voxel.position, uvec3_new(x, y, z)));
            MELT_ASSERT(InnerVoxel(voxel_field[index]));
        }
    }
#endif
}

static void _Debug_ValidateMaxExtents(const MaxExtents& max_extents, const VoxelSet& outer_voxels)
{
#if defined(MELT_DEBUG) && defined(MELT_ASSERT)
    for (u32 i = 0; i < max_extents.size(); ++i)
    {
        const auto& extent = max_extents[i];
        for (u32 x = extent.position.x; x < extent.position.x + extent.extent.x; ++x)
        {
            for (u32 y = extent.position.y; y < extent.position.y + extent.extent.y; ++y)
            {
                for (u32 z = extent.position.z; z < extent.position.z + extent.extent.z; ++z)
                {
                    for (const Voxel& voxel : outer_voxels)
                    {
                        MELT_ASSERT(!uvec3_equals(voxel.position, uvec3_new(x, y, z)));
                    }
                }
            }
        }
    }
#endif
}

static void UpdateMinDistanceField(const Context& context, const uvec3_t& start_position, const uvec3_t& extent, const VoxelField& voxel_field, MinDistanceField& distance_field)
{
    MELT_PROFILE_BEGIN();
    MELT_ASSERT(start_position.x - 1 != ~0U);
    MELT_ASSERT(start_position.y - 1 != ~0U);
    MELT_ASSERT(start_position.z - 1 != ~0U);

    for (u32 x = start_position.x - 1; x != ~0U; --x)
    {
        for (u32 y = start_position.y; y < start_position.y + extent.y; ++y)
        {
            for (u32 z = start_position.z; z < start_position.z + extent.z; ++z)
            {
                const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
                if (InnerVoxel(voxel_field[index]))
                {
                    MinDistance& min_distance = distance_field[index];
                    const u32 updated_distance_x = start_position.x - min_distance.x;
                    min_distance.dist.x = u32_min(updated_distance_x, min_distance.dist.x);
                }
            }
        }
    }
    for (u32 x = start_position.x; x < start_position.x + extent.x; ++x)
    {
        for (u32 y = start_position.y - 1; y != ~0U; --y)
        {
            for (u32 z = start_position.z; z < start_position.z + extent.z; ++z)
            {
                const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
                if (InnerVoxel(voxel_field[index]))
                {
                    MinDistance& min_distance = distance_field[index];
                    const u32 updated_distance_y = start_position.y - min_distance.y;
                    min_distance.dist.y = u32_min(updated_distance_y, min_distance.dist.y);
                }
            }
        }
    }
    for (u32 x = start_position.x; x < start_position.x + extent.x; ++x)
    {
        for (u32 y = start_position.y; y < start_position.y + extent.y; ++y)
        {
            for (u32 z = start_position.z - 1; z != ~0U; --z)
            {
                const u32 index = Flatten3d(uvec3_new(x, y, z), context.dimension);
                if (InnerVoxel(voxel_field[index]))
                {
                    MinDistance& min_distance = distance_field[index];
                    const u32 updated_distance_z = start_position.z - min_distance.z;
                    min_distance.dist.z = u32_min(updated_distance_z, min_distance.dist.z);
                }
            }
        }
    }

    MELT_PROFILE_END();
}

static MeltOccluderBoxType SelectVoxelCubeIndices(MeltOccluderBoxTypeFlags box_type_flags, const u16*& out_indices, u16& out_index_length)
{
#define CHECK_FLAG(flag) (box_type_flags & flag) == flag

    if (CHECK_FLAG(MeltOccluderBoxTypeRegular))
    {
        out_indices = static_cast<const u16*>(VoxelCubeIndices);
        out_index_length = ARRAY_LENGTH(VoxelCubeIndices);
        return MeltOccluderBoxTypeRegular;
    }
    else if (CHECK_FLAG(MeltOccluderBoxTypeSides))
    {
        out_indices = static_cast<const u16*>(VoxelCubeIndicesSides);
        out_index_length = ARRAY_LENGTH(VoxelCubeIndicesSides);
        return MeltOccluderBoxTypeSides;
    }
    else if (CHECK_FLAG(MeltOccluderBoxTypeBottom))
    {
        out_indices = static_cast<const u16*>(VoxelCubeIndicesBottom);
        out_index_length = ARRAY_LENGTH(VoxelCubeIndicesBottom);
        return MeltOccluderBoxTypeBottom;
    }
    else if (CHECK_FLAG(MeltOccluderBoxTypeTop))
    {
        out_indices = static_cast<const u16*>(VoxelCubeIndicesTop);
        out_index_length = ARRAY_LENGTH(VoxelCubeIndicesTop);
        return MeltOccluderBoxTypeTop;
    }
    else if (CHECK_FLAG(MeltOccluderBoxTypeDiagonals))
    {
        out_indices = static_cast<const u16*>(VoxelCubeIndicesDiagonals);
        out_index_length = ARRAY_LENGTH(VoxelCubeIndicesDiagonals);
        return MeltOccluderBoxTypeDiagonals;
    }

#undef CHECK_FLAG

    return MeltOccluderBoxTypeNone;
}

static void AddVoxelToMesh(vec3_t voxel_center, vec3_t half_voxel_size, MeltMesh& mesh, MeltOccluderBoxTypeFlags box_type_flags)
{
    u16 index_offset = (u16)mesh.vertices.size();
        
    for (u32 i = 0; i < ARRAY_LENGTH(VoxelCubeVertices); ++i)
    {
        vec3_t vertex = vec3_add(vec3_mul(half_voxel_size, VoxelCubeVertices[i]), voxel_center);
        mesh.vertices.push_back(vertex);
    }

    while (box_type_flags != MeltOccluderBoxTypeNone)
    {
        const u16* indices = nullptr;
        u16 indices_length = 0;

        MeltOccluderBoxType selected_type = SelectVoxelCubeIndices(box_type_flags, indices, indices_length);

        MELT_ASSERT(indices && indices_length > 0);
        for (u32 i = 0; i < indices_length; ++i)
        {
            mesh.indices.push_back(indices[i] + index_offset);
        }

        box_type_flags &= ~selected_type;
    }
}

#if defined(MELT_DEBUG)
static void AddVoxelToMeshColor(vec3_t voxel_center, vec3_t half_voxel_size, MeltMesh& mesh, MeltOccluderBoxTypeFlags box_type_flags = MeltOccluderBoxTypeRegular, const Color3u8& color = Blueviolet)
{
    u16 index_offset = (u16)mesh.vertices.size() / 2;

    for (u32 i = 0; i < ARRAY_LENGTH(VoxelCubeVertices); ++i)
    {
        vec3_t vertex = vec3_add(vec3_mul(half_voxel_size, VoxelCubeVertices[i]), voxel_center);
        mesh.vertices.push_back(vertex);
        mesh.vertices.push_back(vec3_div(vec3_new(color), 255.0f));
    }

    while (box_type_flags != MeltOccluderBoxTypeNone)
    {
        const u16* indices = nullptr;
        u16 indices_length = 0;

        MeltOccluderBoxType selected_type = SelectVoxelCubeIndices(box_type_flags, indices, indices_length);

        MELT_ASSERT(indices && indices_length > 0);
        for (u32 i = 0; i < indices_length; ++i)
        {
            mesh.indices.push_back(indices[i] + index_offset);
        }

        box_type_flags &= ~selected_type;
    }
}

static void AddVoxelSetToMesh(const VoxelSet& voxel_set, const vec3_t& half_voxel_extent, MeltMesh& mesh)
{
    for (const Voxel& voxel : voxel_set)
    {
        AddVoxelToMeshColor(aabb_center(voxel.aabb), half_voxel_extent, mesh);
    }
}
#endif

static MaxExtent GetMaxExtent(const Context& context, const VoxelField& voxel_field, const MinDistanceField& distance_field)
{
    MELT_PROFILE_BEGIN();

    MaxExtent max_extent;
    max_extent.extent = uvec3_new(0, 0, 0);
    max_extent.position = uvec3_new(0, 0, 0);
    max_extent.volume = 0;

    for (auto& min_distance : distance_field)
    {
        VoxelStatus voxel_status = voxel_field[Flatten3d(min_distance.position, context.dimension)];
        if (InnerVoxel(voxel_status))
        {
            uvec3_t extent = GetMaxAABBExtent(context, distance_field, voxel_field, min_distance);
            u32 volume = extent.x * extent.y * extent.z;
            if (volume > max_extent.volume)
            {
                max_extent.extent = extent;
                max_extent.position = min_distance.position;
                max_extent.volume = max_extent.extent.x * max_extent.extent.y * max_extent.extent.z;
            }
        }
    }

    MELT_PROFILE_END();

    return max_extent;
}

bool MeltGenerateOccluder(const MeltParams& params, MeltResult& out_result)
{
    out_result.debugMesh.vertices.clear();
    out_result.debugMesh.indices.clear();
    out_result.mesh.vertices.clear();
    out_result.mesh.indices.clear();

    vec3_t voxel_extent = vec3_new(params.voxelSize, params.voxelSize, params.voxelSize);
    vec3_t half_voxel_extent = vec3_mul(voxel_extent, 0.5f);

    AABB mesh_aabb = GenerateAABB(params.mesh);

    mesh_aabb.min = vec3_sub(MapToVoxelMinBound(mesh_aabb.min, params.voxelSize), voxel_extent);
    mesh_aabb.max = vec3_add(MapToVoxelMaxBound(mesh_aabb.max, params.voxelSize), voxel_extent);

    vec3_t mesh_extent = vec3_sub(mesh_aabb.max, mesh_aabb.min);
    vec3_t inv_mesh_extent = vec3_new(1.0f / mesh_extent.x, 1.0f / mesh_extent.y, 1.0f / mesh_extent.z);
    vec3_t voxel_count = vec3_div(mesh_extent, params.voxelSize);
    vec3_t voxel_resolution = vec3_mul(voxel_count, inv_mesh_extent);

    Context context;
    context.dimension = uvec3_new(voxel_count);
    context.size = u32(voxel_count.x) * u32(voxel_count.y) * u32(voxel_count.z);

    VoxelIndices voxel_indices(context.size, -1);
    VoxelSet outer_voxels;
    outer_voxels.reserve(context.size);

    // Perform shell voxelization
    for (u32 i = 0; i < params.mesh.indices.size(); i += 3)
    {
        MELT_PROFILE_BEGIN();

        Triangle triangle;

        triangle.v0 = params.mesh.vertices[params.mesh.indices[i + 0]];
        triangle.v1 = params.mesh.vertices[params.mesh.indices[i + 1]];
        triangle.v2 = params.mesh.vertices[params.mesh.indices[i + 2]];

        AABB triangle_aabb = GenerateAABB(triangle);

        // Voxel snapping, snap the triangle extent to find the 3d grid to iterate on.
        triangle_aabb.min = vec3_sub(MapToVoxelMinBound(triangle_aabb.min, params.voxelSize), voxel_extent);
        triangle_aabb.max = vec3_add(MapToVoxelMaxBound(triangle_aabb.max, params.voxelSize), voxel_extent);

        for (f32 x = triangle_aabb.min.x; x <= triangle_aabb.max.x; x += params.voxelSize)
        {
            for (f32 y = triangle_aabb.min.y; y <= triangle_aabb.max.y; y += params.voxelSize)
            {
                for (f32 z = triangle_aabb.min.z; z <= triangle_aabb.max.z; z += params.voxelSize)
                {
                    Voxel voxel;

                    voxel.aabb.min = vec3_sub(vec3_new(x, y, z), half_voxel_extent);
                    voxel.aabb.max = vec3_add(vec3_new(x, y, z), half_voxel_extent);

                    MELT_ASSERT(voxel.aabb.min.x >= mesh_aabb.min.x - half_voxel_extent.x);
                    MELT_ASSERT(voxel.aabb.min.y >= mesh_aabb.min.y - half_voxel_extent.y);
                    MELT_ASSERT(voxel.aabb.min.z >= mesh_aabb.min.z - half_voxel_extent.z);

                    MELT_ASSERT(voxel.aabb.max.x <= mesh_aabb.max.x + half_voxel_extent.x);
                    MELT_ASSERT(voxel.aabb.max.y <= mesh_aabb.max.y + half_voxel_extent.y);
                    MELT_ASSERT(voxel.aabb.max.z <= mesh_aabb.max.z + half_voxel_extent.z);

                    vec3_t voxel_center = aabb_center(voxel.aabb);
                    vec3_t relative_to_origin = vec3_sub(vec3_sub(voxel_center, mesh_aabb.min), half_voxel_extent);

                    if (!AABBIntersectsTriangle(triangle, voxel_center, half_voxel_extent))
                        continue;
                    voxel.position = uvec3_new(vec3_mul(relative_to_origin, voxel_resolution));

                    const u32 index = Flatten3d(voxel.position, context.dimension);
                    if (voxel_indices[index] != -1)
                        continue;

                    voxel_indices[index] = (s32)outer_voxels.size();
                    outer_voxels.push_back(voxel);
                }
            }
        }

        MELT_PROFILE_END();
    }

    VoxelSetPlanes voxel_set_planes;

    // Generate a flat voxel list per plane (x,y), (x,z), (y,z)
    GeneratePerPlaneVoxelSet(context, outer_voxels, voxel_indices, voxel_set_planes);

    // The minimum distance field is a data structure representing, for each voxel,
    // the minimum distance that we can go in each of the positive directions x, y,
    // z until we collide with a shell voxel. The voxel field is a data structure
    // representing the state of the voxels. The clip status is a representation of
    // whether a voxel is in the clip state, whether a voxel can 'see' a voxel in
    // each of the directions +x,-x,+y,-y,+z,-z, and whether the voxel is an 'inner'
    // voxel (contained within the shell voxels).

    MinDistanceField min_distance_field(context.size);

    VoxelField voxel_field(context.size);

    // Generate the minimum distance field, and voxel status from the initial shell.

    GenerateFields(context, voxel_set_planes, voxel_indices, voxel_field, min_distance_field);

    if (!WaterTightMesh(context, voxel_field, min_distance_field))
    {
        return false;
    }

    _Debug_ValidateMinDistanceField(context, outer_voxels, voxel_field, min_distance_field);

    std::vector<MaxExtent> max_extents;
    if (params.debug.extentMaxStep != 0)
    {
        for (s32 i = 0; i < params.debug.extentMaxStep; ++i)
        {
            MaxExtent max_extent = GetMaxExtent(context, voxel_field, min_distance_field);

            ClipVoxelField(context, max_extent.position, max_extent.extent, voxel_field);

            UpdateMinDistanceField(context, max_extent.position, max_extent.extent, voxel_field, min_distance_field);

            _Debug_ValidateMinDistanceField(context, outer_voxels, voxel_field, min_distance_field);

            max_extents.push_back(max_extent);
        }
    }
    else
    {
        u32 volume = 0;
        u32 total_volume = 0;
        f32 fillPercentage = 0.0f;

        // Approximate the volume of the mesh by the number of voxels that can fit within.
        for (auto& min_distance : min_distance_field)
        {
            VoxelStatus voxel_status = voxel_field[Flatten3d(min_distance.position, context.dimension)];

            // Each inner voxel adds one unit to the volume.
            if (InnerVoxel(voxel_status))
                ++total_volume;
        }

        // One iteration to find an extent does the following:
        // . Get the extent that maximizes the volume considering the minimum distance
        //    field
        // . Clip the max extent found to the set of inner voxels
        // . Update the minimum distance field by adjusting the distances on the set
        //    of inner voxels. This is done by extending the extent cube to infinity
        //    on each of the axes +x, +y, +z
        while (fillPercentage < params.fillPercentage && volume != total_volume)
        {
            MaxExtent max_extent = GetMaxExtent(context, voxel_field, min_distance_field);

            ClipVoxelField(context, max_extent.position, max_extent.extent, voxel_field);

            UpdateMinDistanceField(context, max_extent.position, max_extent.extent, voxel_field, min_distance_field);

            _Debug_ValidateMinDistanceField(context, outer_voxels, voxel_field, min_distance_field);

            max_extents.push_back(max_extent);

            fillPercentage += f32(max_extent.volume) / total_volume;
            volume += max_extent.volume;
        }
    }

    for (u32 i = 0; i < max_extents.size(); ++i)
    {
        const MaxExtent& extent = max_extents[i];

        vec3_t half_extent = vec3_mul(vec3_new(extent.extent), half_voxel_extent);
        vec3_t voxel_position = vec3_mul(vec3_new(extent.position), voxel_extent);
        vec3_t voxel_position_biased_to_center = vec3_add(voxel_position, half_extent);
        vec3_t aabb_center = vec3_add(mesh_aabb.min, voxel_position_biased_to_center);

        AddVoxelToMesh(vec3_add(aabb_center, half_voxel_extent), half_extent, out_result.mesh, params.boxTypeFlags);
    }

    _Debug_ValidateMaxExtents(max_extents, outer_voxels);

#if defined(MELT_DEBUG)
    MELT_ASSERT(params.debug._start_canary == 0 && params.debug._end_canary == 0 && "Make sure to memset MeltDebugParams to 0 before use");

    if (params.debug.flags > 0)
    {
        // AddVoxelToMeshColor(aabb_center(mesh_aabb), (mesh_aabb.max - mesh_aabb.min) * 0.5f, out_result.debugMesh, Colors[0]);

        if (params.debug.flags & MeltDebugTypeShowOuter)
        {
            AddVoxelSetToMesh(outer_voxels, vec3_mul(half_voxel_extent, params.debug.voxelScale), out_result.debugMesh);
        }
        if (params.debug.flags & MeltDebugTypeShowSliceSelection)
        {
            if (params.debug.voxelY > 0 && params.debug.voxelZ > 0)
            {
                u32 index = Flatten2d(uvec2_new(params.debug.voxelY, params.debug.voxelZ), uvec2_new(context.dimension.y, context.dimension.z));
                VoxelSet& voxels_x_planes = voxel_set_planes.x[index];
                AddVoxelSetToMesh(voxels_x_planes, vec3_mul(half_voxel_extent, params.debug.voxelScale), out_result.debugMesh);
            }
            if (params.debug.voxelX > 0 && params.debug.voxelZ > 0)
            {
                u32 index = Flatten2d(uvec2_new(params.debug.voxelX, params.debug.voxelZ), uvec2_new(context.dimension.x, context.dimension.z));
                VoxelSet& voxels_y_planes = voxel_set_planes.y[index];
                AddVoxelSetToMesh(voxels_y_planes, vec3_mul(half_voxel_extent, params.debug.voxelScale), out_result.debugMesh);
            }
            if (params.debug.voxelX > 0 && params.debug.voxelY > 0)
            {
                u32 index = Flatten2d(uvec2_new(params.debug.voxelX, params.debug.voxelY), uvec2_new(context.dimension.x, context.dimension.y));
                VoxelSet& voxels_z_planes = voxel_set_planes.z[index];
                AddVoxelSetToMesh(voxels_z_planes, vec3_mul(half_voxel_extent, params.debug.voxelScale), out_result.debugMesh);
            }
        }
        if (params.debug.flags & MeltDebugTypeShowInner)
        {
            for (const auto& min_distance : min_distance_field)
            {
                const u32 index = Flatten3d(min_distance.position, context.dimension);
                if (!voxel_field[index].inner)
                    continue;

                vec3_t voxel_position = vec3_mul(vec3_new(min_distance.position), voxel_extent);
                vec3_t voxel_center = vec3_add(mesh_aabb.min, voxel_position);
                if (params.debug.voxelX < 0 ||
                    params.debug.voxelY < 0 ||
                    params.debug.voxelZ < 0)
                {
                    AddVoxelToMeshColor(vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debugMesh);
                }
            }
        }
        if (params.debug.flags & MeltDebugTypeShowMinDistance)
        {
            for (const auto& min_distance : min_distance_field)
            {
                vec3_t voxel_center = vec3_add(mesh_aabb.min, vec3_mul(vec3_new(min_distance.position), voxel_extent));
                if (u32(params.debug.voxelX) == min_distance.x &&
                    u32(params.debug.voxelY) == min_distance.y &&
                    u32(params.debug.voxelZ) == min_distance.z)
                {
                    AddVoxelToMeshColor(vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debugMesh);

                    for (u32 x = min_distance.x; x < min_distance.x + min_distance.dist.x; ++x)
                    {
                        vec3_t voxel_center = vec3_add(mesh_aabb.min, vec3_mul(vec3_new(x, min_distance.y, min_distance.z), voxel_extent));
                        AddVoxelToMeshColor(vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debugMesh);
                    }
                    for (u32 y = min_distance.y; y < min_distance.y + min_distance.dist.y; ++y)
                    {
                        vec3_t voxel_center = vec3_add(mesh_aabb.min, vec3_mul(vec3_new(min_distance.x, y, min_distance.z), voxel_extent));
                        AddVoxelToMeshColor(vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debugMesh);
                    }
                    for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
                    {
                        vec3_t voxel_center = vec3_add(mesh_aabb.min, vec3_mul(vec3_new(min_distance.x, min_distance.y, z), voxel_extent));
                        AddVoxelToMeshColor(vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debugMesh);
                    }
                }
            }
        }
        if (params.debug.flags & MeltDebugTypeShowExtent)
        {
            for (const auto& min_distance : min_distance_field)
            {
                uvec3_t max_extent = GetMaxAABBExtent(context, min_distance_field, voxel_field, min_distance);
                for (u32 x = min_distance.x; x < min_distance.x + max_extent.x; ++x)
                {
                    for (u32 y = min_distance.y; y < min_distance.y + max_extent.y; ++y)
                    {
                        for (u32 z = min_distance.z; z < min_distance.z + max_extent.z; ++z)
                        {
                            vec3_t voxel_center = vec3_add(mesh_aabb.min, vec3_mul(vec3_new(x, y, z), voxel_extent));
                            AddVoxelToMeshColor(vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debugMesh);
                        }
                    }
                }
            }
        }
        if (params.debug.flags & MeltDebugTypeShowResult)
        {
            for (size_t i = 0; i < max_extents.size(); ++i)
            {
                const MaxExtent& extent = max_extents[i];
                if (s32(i) == params.debug.extentIndex || params.debug.extentIndex < 0)
                {
                    vec3_t half_extent = vec3_mul(vec3_new(extent.extent), half_voxel_extent);
                    vec3_t voxel_position = vec3_mul(vec3_new(extent.position), voxel_extent);
                    vec3_t voxel_position_biased_to_center = vec3_add(voxel_position, half_extent);
                    vec3_t aabb_center = vec3_add(mesh_aabb.min, voxel_position_biased_to_center);
                    Color3u8 color = Colors[i % ARRAY_LENGTH(Colors)];
                    AddVoxelToMeshColor(vec3_add(aabb_center, half_voxel_extent), half_extent, out_result.debugMesh, params.boxTypeFlags, color);
                }
            }
        }
    }
#endif

    return true;
}

#endif // MELT_IMPLEMENTATION
