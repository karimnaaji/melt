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

#include "melt.h"

namespace Melt
{

typedef unsigned int    u32;
typedef unsigned short  u16;
typedef unsigned char   u8;
typedef signed long     s64;
typedef signed short    s16;
typedef signed int      s32;
typedef signed char     s8;
typedef float           f32;

typedef glm::vec3       vec3;
typedef glm::ivec3      svec3;
typedef glm::uvec2      uvec2;
typedef glm::uvec3      uvec3;

template<typename T>
static T Min(T a, T b)
{
    return a < b ? a : b;
}

typedef uvec3 Color3u8;

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

struct AABB
{
    vec3 min;
    vec3 max;

    inline vec3 Center() const
    {
        return (min + max) * 0.5f;
    }
};

struct Voxel
{
    AABB aabb;
    uvec3 position;
};

struct Triangle
{
    vec3 v0;
    vec3 v1;
    vec3 v2;
};

struct Plane
{
    vec3 normal;
    float distance;
};

struct MinDistance
{
    MinDistance() {}
    svec3 dist;
    union
    {
        struct
        {
            u32 x;
            u32 y;
            u32 z;
        };
        uvec3 position;
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

struct Context
{
    uvec3 dimension;
    u32 size;
};

struct MaxExtent
{
    uvec3 position;
    uvec3 extent;
    u32 volume;
};

using VoxelIndices = std::vector<s32>;
using VoxelSet = std::vector<Voxel>;
using MinDistanceField = std::vector<MinDistance>;
using VoxelField = std::vector<VoxelStatus>;
using MaxExtents = std::vector<MaxExtent>;

#define ARRAY_LENGTH(array) ((int)(sizeof(array) / sizeof(*array)))

static bool AABBCollidesPlane(const Plane plane, const vec3& half_aabb_dim)
{
    vec3 vmin;
    vec3 vmax;

    for (u32 dimension = 0; dimension < 3; ++dimension)
    {
        if (plane.normal[dimension] > 0.0f)
        {
            vmin[dimension] = -half_aabb_dim[dimension];
            vmax[dimension] =  half_aabb_dim[dimension];
        }
        else
        {
            vmin[dimension] =  half_aabb_dim[dimension];
            vmax[dimension] = -half_aabb_dim[dimension];
        }
    }
    if (glm::dot(plane.normal, vmin) + plane.distance > 0.0f)
        return false;

    if (glm::dot(plane.normal, vmax) + plane.distance >= 0.0f)
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

static bool AABBIntersectsTriangle(const Triangle& triangle, const vec3& aabb_center, const vec3& half_aabb_dim)
{
    vec3 v0, v1, v2;
    vec3 e0, e1, e2;
    vec3 edge_abs;
    f32 min, max;
    f32 p0, p1, p2;
    f32 rad;

    v0 = triangle.v0 - aabb_center;
    v1 = triangle.v1 - aabb_center;
    v2 = triangle.v2 - aabb_center;

    e0 = v1 - v0;
    e1 = v2 - v1;
    e2 = v0 - v2;

    edge_abs = glm::abs(e0);
    AXISTEST_X01(e0.z, e0.y, edge_abs.z, edge_abs.y);
    AXISTEST_Y02(e0.z, e0.x, edge_abs.z, edge_abs.x);
    AXISTEST_Z12(e0.y, e0.x, edge_abs.y, edge_abs.x);

    edge_abs = glm::abs(e1);
    AXISTEST_X01(e1.z, e1.y, edge_abs.z, edge_abs.y);
    AXISTEST_Y02(e1.z, e1.x, edge_abs.z, edge_abs.x);
    AXISTEST_Z0 (e1.y, e1.x, edge_abs.y, edge_abs.x);

    edge_abs = glm::abs(e2);
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
    plane.normal = glm::cross(e0, e1);
    plane.distance = -glm::dot(plane.normal, v0);

    if (!AABBCollidesPlane(plane, half_aabb_dim))
        return false;

    return true;
}

static inline u32 Flatten(const uvec3& index, const uvec3& dimension)
{
    u32 out_index = index.x + dimension.x * index.y + dimension.x * dimension.y * index.z;
    return out_index;
}

static inline uvec3 UnFlatten(u32 index, const uvec3& dimension)
{
    uvec3 out_index;

    u32 dim_xy = dimension.x * dimension.y;
    out_index.z = index / dim_xy;
    index -= out_index.z * dim_xy;
    out_index.y = index / dimension.x;
    out_index.x = index % dimension.x;
    return out_index;
}

static vec3 MapToVoxelMaxBound(const vec3& position, f32 voxel_size)
{
    auto MapFunc = [&](f32 value)
    {
        f32 sign = value < 0.0f ? -1.0f : 1.0f;
        f32 result = value + sign * voxel_size * 0.5f;
        return ceil(result / voxel_size) * voxel_size;
    };

    vec3 mapped_position;

    mapped_position.x = MapFunc(position.x);
    mapped_position.y = MapFunc(position.y);
    mapped_position.z = MapFunc(position.z);

    return mapped_position;
}

static vec3 MapToVoxelMinBound(const vec3& position, f32 voxel_size)
{
    auto MapFunc = [&](f32 value)
    {
        f32 sign = value < 0.0f ? -1.0f : 1.0f;
        f32 result = value + sign * voxel_size * 0.5f;
        return floorf(result / voxel_size) * voxel_size;
    };

    vec3 mapped_position;

    mapped_position.x = MapFunc(position.x);
    mapped_position.y = MapFunc(position.y);
    mapped_position.z = MapFunc(position.z);

    return mapped_position;
}

static AABB GenerateAABB(const Triangle& triangle)
{
    AABB aabb;

    aabb.min = vec3( FLT_MAX);
    aabb.max = vec3(-FLT_MAX);

    aabb.min = glm::min(aabb.min, triangle.v0);
    aabb.max = glm::max(aabb.max, triangle.v0);

    aabb.min = glm::min(aabb.min, triangle.v1);
    aabb.max = glm::max(aabb.max, triangle.v1);

    aabb.min = glm::min(aabb.min, triangle.v2);
    aabb.max = glm::max(aabb.max, triangle.v2);

    return aabb;
}

static AABB GenerateAABB(const Mesh& mesh)
{
    AABB aabb;

    aabb.min = vec3( FLT_MAX);
    aabb.max = vec3(-FLT_MAX);

    for (u32 i = 0; i < mesh.indices.size(); ++i)
    {
        aabb.min = glm::min(aabb.min, mesh.vertices[mesh.indices[i]]);
        aabb.max = glm::max(aabb.max, mesh.vertices[mesh.indices[i]]);
    }

    return aabb;
}

static VoxelSet SelectVoxelsOnXAxis(const VoxelSet& voxel_set, u32 y, u32 z)
{
    VoxelSet result;
    for (u32 i = 0; i < voxel_set.size(); ++i)
    {
        if (voxel_set[i].position.y == y && voxel_set[i].position.z == z)
            result.push_back(voxel_set[i]);
    }
    return result;
}

static VoxelSet SelectVoxelsOnYAxis(const VoxelSet& voxel_set, u32 x, u32 z)
{
    VoxelSet result;
    for (u32 i = 0; i < voxel_set.size(); ++i)
    {
        if (voxel_set[i].position.x == x && voxel_set[i].position.z == z)
            result.push_back(voxel_set[i]);
    }
    return result;
}

static VoxelSet SelectVoxelsOnZAxis(const VoxelSet& voxel_set, u32 x, u32 y)
{
    VoxelSet result;
    for (u32 i = 0; i < voxel_set.size(); ++i)
    {
        if (voxel_set[i].position.x == x && voxel_set[i].position.y == y)
            result.push_back(voxel_set[i]);
    }
    return result;
}

static void GetField(const VoxelSet& voxel_set, u32 x, u32 y, u32 z, MinDistance& out_min_distance, VoxelStatus& out_status)
{
    const svec3 InfiniteDistance = svec3(INT_MAX);
    const svec3 NullDistance = svec3(0);

    out_min_distance.dist = InfiniteDistance;
    out_min_distance.x = x;
    out_min_distance.y = y;
    out_min_distance.z = z;

    out_status.visibility = VisibilityNull;
    out_status.clipped = false;
    out_status.inner = false;

    const VoxelSet voxels_x = SelectVoxelsOnXAxis(voxel_set, y, z);
    for (const auto& voxel : voxels_x)
    {
        s32 distance = voxel.position.x - x;
        if (distance > 0)
        {
            out_status.visibility |= VisibilityPlusX;
            out_min_distance.dist.x = Min(out_min_distance.dist.x, distance);
        }
        else if (distance < 0)
            out_status.visibility |= VisibilityMinusX;
        else
            out_min_distance.dist.x = 0;
    }
    const VoxelSet voxels_y = SelectVoxelsOnYAxis(voxel_set, x, z);
    for (const auto& voxel : voxels_y)
    {
        s32 distance = voxel.position.y - y;
        if (distance > 0)
        {
            out_status.visibility |= VisibilityPlusY;
            out_min_distance.dist.y = Min(out_min_distance.dist.y, distance);
        }
        else if (distance < 0)
            out_status.visibility |= VisibilityMinusY;
        else
            out_min_distance.dist.y = 0;
    }
    const VoxelSet voxels_z = SelectVoxelsOnZAxis(voxel_set, x, y);
    for (const auto& voxel : voxels_z)
    {
        s32 distance = voxel.position.z - z;
        if (distance > 0)
        {
            out_status.visibility |= VisibilityPlusZ;
            out_min_distance.dist.z = Min(out_min_distance.dist.z, distance);
        }
        else if (distance < 0)
            out_status.visibility |= VisibilityMinusZ;
        else
            out_min_distance.dist.z = 0;
    }
    if (out_status.visibility == VisibilityAll)
    {
        if (out_min_distance.dist != InfiniteDistance &&
            out_min_distance.dist != NullDistance)
        {
            out_status.inner = true;
        }
    }
}

static void GenerateFields(const Context& context, const VoxelSet& voxel_set, VoxelField& out_voxel_field, MinDistanceField& out_distance_field)
{
    for (u32 i = 0; i < context.size; ++i)
    {
        MinDistance& min_distance = out_distance_field[i];
        VoxelStatus& voxel_status = out_voxel_field[i];
        const uvec3 position = UnFlatten(i, context.dimension);

        GetField(voxel_set, position.x, position.y, position.z, min_distance, voxel_status);
    }
}

static inline bool InnerVoxel(VoxelStatus voxel_status)
{
    return voxel_status.inner && !voxel_status.clipped;
}

static uvec3 GetMaxAABBExtent(const Context& context, const MinDistanceField& distance_field, const VoxelField& voxel_field, const MinDistance& min_distance)
{
    std::vector<uvec2> max_aabb_extents;

    for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
    {
        uvec3 z_slice_position = uvec3(min_distance.x, min_distance.y, z);
        u32 z_slice_index = Flatten(z_slice_position, context.dimension);

        MELT_ASSERT(voxel_field[z_slice_index].inner);

        if (voxel_field[z_slice_index].clipped)
            continue;

        const MinDistance& sample_min_distance = distance_field[z_slice_index];

        uvec2 max_extent;
        max_extent.x = sample_min_distance.dist.x;
        max_extent.y = sample_min_distance.dist.y;

        u32 x = sample_min_distance.x + 1;
        u32 y = sample_min_distance.y + 1;
        u32 i = 1;
        while (x < sample_min_distance.x + sample_min_distance.dist.x &&
               y < sample_min_distance.y + sample_min_distance.dist.y)
        {
            const u32 index = Flatten(uvec3(x, y, z), context.dimension);
            if (InnerVoxel(voxel_field[index]))
            {
                const MinDistance& distance = distance_field[index];
                max_extent.x = Min(distance.dist.x + i, max_extent.x);
                max_extent.y = Min(distance.dist.y + i, max_extent.y);
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

    uvec2 min_extent(UINT_MAX);

    u32 z_slice = 1;
    u32 max_volume_slice = z_slice;
    u32 max_volume = 0;

    MELT_ASSERT(max_aabb_extents.size() > 0);

    for (const uvec2& extent : max_aabb_extents)
    {
        min_extent.x = Min(extent.x, min_extent.x);
        min_extent.y = Min(extent.y, min_extent.y);

        const u32 volume = min_extent.x * min_extent.y * z_slice;
        if (volume > max_volume)
        {
            max_volume_slice = z_slice;
            max_volume = volume;
        }

        ++z_slice;
    }

    MELT_ASSERT(max_volume > 0);
    MELT_ASSERT(min_extent.x > 0);
    MELT_ASSERT(min_extent.y > 0);
    MELT_ASSERT(z_slice > 1);

    return uvec3(min_extent.x, min_extent.y, z_slice - 1);
}

static void ClipVoxelField(const Context& context, const uvec3& start_position, const uvec3& extent, VoxelField& voxel_field)
{
    for (u32 x = start_position.x; x < start_position.x + extent.x; ++x)
    {
        for (u32 y = start_position.y; y < start_position.y + extent.y; ++y)
        {
            for (u32 z = start_position.z; z < start_position.z + extent.z; ++z)
            {
                u32 index = Flatten(uvec3(x, y, z), context.dimension);
                MELT_ASSERT(!voxel_field[index].clipped && "Clipping already clipped voxel field index");
                voxel_field[index].clipped = true;
            }
        }
    }
}

static void _Debug_ValidateMinDistanceField(const Context& context, const VoxelSet& outer_voxels, const VoxelField& voxel_field, const MinDistanceField& distance_field)
{
#ifndef MELT_ASSERT
    return;
#endif
    for (const auto& min_distance : distance_field)
    {
        if (!InnerVoxel(voxel_field[Flatten(min_distance.position, context.dimension)]))
            continue;

        for (u32 x = min_distance.x; x < min_distance.x + min_distance.dist.x; ++x)
        {
            const u32 y = min_distance.y;
            const u32 z = min_distance.z;
            const u32 index = Flatten(uvec3(x, y, z), context.dimension);
            for (const Voxel& voxel : outer_voxels)
                MELT_ASSERT(voxel.position != uvec3(x, y, z));
            MELT_ASSERT(InnerVoxel(voxel_field[index]));
        }
        for (u32 y = min_distance.y; y < min_distance.y + min_distance.dist.y; ++y)
        {
            const u32 x = min_distance.x;
            const u32 z = min_distance.z;
            const u32 index = Flatten(uvec3(x, y, z), context.dimension);
            for (const Voxel& voxel : outer_voxels)
                MELT_ASSERT(voxel.position != uvec3(x, y, z));
            MELT_ASSERT(InnerVoxel(voxel_field[index]));
        }
        for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
        {
            const u32 x = min_distance.x;
            const u32 y = min_distance.y;
            const u32 index = Flatten(uvec3(x, y, z), context.dimension);
            for (const Voxel& voxel : outer_voxels)
                MELT_ASSERT(voxel.position != uvec3(x, y, z));
            MELT_ASSERT(InnerVoxel(voxel_field[index]));
        }
    }
}

static void _Debug_ValidateMaxExtents(const MaxExtents& max_extents, const VoxelSet& outer_voxels)
{
#ifndef MELT_ASSERT
    return;
#endif
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
                        MELT_ASSERT(voxel.position != uvec3(x, y, z));
                    }
                }
            }
        }
    }
}

static void UpdateMinDistanceField(const Context& context, const uvec3& start_position, const uvec3& extent, const VoxelField& voxel_field, MinDistanceField& distance_field)
{
    MELT_ASSERT(start_position.x - 1 != ~0U);
    MELT_ASSERT(start_position.y - 1 != ~0U);
    MELT_ASSERT(start_position.z - 1 != ~0U);

    for (u32 x = start_position.x - 1; x != ~0U; --x)
    {
        for (u32 y = start_position.y; y < start_position.y + extent.y; ++y)
        {
            for (u32 z = start_position.z; z < start_position.z + extent.z; ++z)
            {
                const u32 index = Flatten(uvec3(x, y, z), context.dimension);
                if (InnerVoxel(voxel_field[index]))
                {
                    MinDistance& min_distance = distance_field[index];
                    const u32 updated_distance_x = start_position.x - min_distance.x;
                    min_distance.dist.x = Min<u32>(updated_distance_x, min_distance.dist.x);
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
                const u32 index = Flatten(uvec3(x, y, z), context.dimension);
                if (InnerVoxel(voxel_field[index]))
                {
                    MinDistance& min_distance = distance_field[index];
                    const u32 updated_distance_y = start_position.y - min_distance.y;
                    min_distance.dist.y = Min<u32>(updated_distance_y, min_distance.dist.y);
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
                const u32 index = Flatten(uvec3(x, y, z), context.dimension);
                if (InnerVoxel(voxel_field[index]))
                {
                    MinDistance& min_distance = distance_field[index];
                    const u32 updated_distance_z = start_position.z - min_distance.z;
                    min_distance.dist.z = Min<u32>(updated_distance_z, min_distance.dist.z);
                }
            }
        }
    }
}

static void AddVoxelToMesh(vec3 voxel_center, vec3 half_voxel_size, Mesh& mesh, const Color3u8& color = Blueviolet)
{
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

    const vec3 VoxelCubeVertices[8] =
    {
        vec3(-1.0f,  1.0f,  1.0f),
        vec3(-1.0f, -1.0f,  1.0f),
        vec3( 1.0f, -1.0f,  1.0f),
        vec3( 1.0f,  1.0f,  1.0f),
        vec3(-1.0f,  1.0f, -1.0f),
        vec3(-1.0f, -1.0f, -1.0f),
        vec3( 1.0f, -1.0f, -1.0f),
        vec3( 1.0f,  1.0f, -1.0f),
    };

    u16 index_offset = (u16)mesh.vertices.size() / 2;
    for (u32 i = 0; i < 8; ++i)
    {
        mesh.vertices.push_back(half_voxel_size * VoxelCubeVertices[i] + voxel_center);
        mesh.vertices.push_back(vec3(color) / 255.0f);
    }
    for (u32 i = 0; i < ARRAY_LENGTH(VoxelCubeIndices); ++i)
    {
        mesh.indices.push_back(VoxelCubeIndices[i] + index_offset);
    }
}

static void AddVoxelSetToMesh(const VoxelSet& voxel_set, const vec3& half_voxel_extent, Mesh& mesh)
{
    for (const Voxel& voxel : voxel_set)
    {
        AddVoxelToMesh(voxel.aabb.Center(), half_voxel_extent, mesh);
    }
}

static MaxExtent GetMaxExtent(const Context& context, const VoxelField& voxel_field, const MinDistanceField& distance_field)
{
    MaxExtent max_extent;
    max_extent.extent = uvec3(0);
    max_extent.position = uvec3(0);
    max_extent.volume = 0;

    for (auto& min_distance : distance_field)
    {
        VoxelStatus voxel_status = voxel_field[Flatten(min_distance.position, context.dimension)];
        if (InnerVoxel(voxel_status))
        {
            uvec3 extent = GetMaxAABBExtent(context, distance_field, voxel_field, min_distance);
            u32 volume = extent.x * extent.y * extent.z;
            if (volume > max_extent.volume)
            {
                max_extent.extent = extent;
                max_extent.position = min_distance.position;
                max_extent.volume = max_extent.extent.x * max_extent.extent.y * max_extent.extent.z;
            }
        }
    }

    return max_extent;
}

void GenerateConservativeOccluder(const Mesh& mesh, const OccluderGenerationParams& gen_params, const DebugParams& debug_params, Mesh& out_mesh)
{
    out_mesh.vertices.clear();
    out_mesh.indices.clear();

    vec3 voxel_extent(gen_params.voxelSize);
    vec3 half_voxel_extent = voxel_extent * 0.5f;

    AABB mesh_aabb = GenerateAABB(mesh);

    mesh_aabb.min = MapToVoxelMinBound(mesh_aabb.min, gen_params.voxelSize) - voxel_extent;
    mesh_aabb.max = MapToVoxelMaxBound(mesh_aabb.max, gen_params.voxelSize) + voxel_extent;

    vec3 mesh_extent = mesh_aabb.max - mesh_aabb.min;
    vec3 voxel_count = mesh_extent / gen_params.voxelSize;

    Context context;
    context.dimension = uvec3(voxel_count);
    context.size = u32(voxel_count.x * voxel_count.y * voxel_count.z);

    VoxelIndices voxel_indices(context.size, -1);
    VoxelSet outer_voxels;
    outer_voxels.reserve(context.size);

    // Perform shell voxelization
    for (u32 i = 0; i < mesh.indices.size(); i += 3)
    {
        Triangle triangle;

        triangle.v0 = mesh.vertices[mesh.indices[i + 0]];
        triangle.v1 = mesh.vertices[mesh.indices[i + 1]];
        triangle.v2 = mesh.vertices[mesh.indices[i + 2]];

        AABB triangle_aabb = GenerateAABB(triangle);

        // Voxel snapping, snap the triangle extent to find the 3d grid to iterate on.
        triangle_aabb.min = MapToVoxelMinBound(triangle_aabb.min, gen_params.voxelSize) - voxel_extent;
        triangle_aabb.max = MapToVoxelMaxBound(triangle_aabb.max, gen_params.voxelSize) + voxel_extent;

        for (f32 x = triangle_aabb.min.x; x <= triangle_aabb.max.x; x += gen_params.voxelSize)
        {
            for (f32 y = triangle_aabb.min.y; y <= triangle_aabb.max.y; y += gen_params.voxelSize)
            {
                for (f32 z = triangle_aabb.min.z; z <= triangle_aabb.max.z; z += gen_params.voxelSize)
                {
                    Voxel voxel;

                    voxel.aabb.min = vec3(x, y, z) - half_voxel_extent;
                    voxel.aabb.max = vec3(x, y, z) + half_voxel_extent;

                    vec3 voxel_center = voxel.aabb.Center();
                    vec3 relative_to_origin = voxel_center - mesh_aabb.min - half_voxel_extent;

                    if (!AABBIntersectsTriangle(triangle, voxel_center, half_voxel_extent))
                        continue;

                    voxel.position = uvec3((relative_to_origin / mesh_extent) * voxel_count);

                    const u32 index = Flatten(voxel.position, context.dimension);
                    if (voxel_indices[index] != -1)
                        continue;

                    voxel_indices[index] = (s32)outer_voxels.size();
                    outer_voxels.push_back(voxel);
                }
            }
        }
    }

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

    GenerateFields(context, outer_voxels, voxel_field, min_distance_field);

    _Debug_ValidateMinDistanceField(context, outer_voxels, voxel_field, min_distance_field);

    std::vector<MaxExtent> max_extents;
    if (debug_params.extentMaxStep != 0)
    {
        for (u32 i = 0; i < debug_params.extentMaxStep; ++i)
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
            VoxelStatus voxel_status = voxel_field[Flatten(min_distance.position, context.dimension)];

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
        while (fillPercentage < gen_params.fillPercentage && volume != total_volume)
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

    _Debug_ValidateMaxExtents(max_extents, outer_voxels);

    if (debug_params.flags > 0)
    {
        if (debug_params.flags & DebugTypeShowOuter)
        {
            AddVoxelSetToMesh(outer_voxels, half_voxel_extent * debug_params.voxelScale, out_mesh);
        }
        if (debug_params.flags & DebugTypeShowSliceSelection)
        {
            if (debug_params.voxelY > 0 && debug_params.voxelZ > 0)
            {
                VoxelSet voxels_x = SelectVoxelsOnXAxis(outer_voxels, debug_params.voxelY, debug_params.voxelZ);
                AddVoxelSetToMesh(voxels_x, half_voxel_extent * debug_params.voxelScale, out_mesh);
            }
            if (debug_params.voxelX > 0 && debug_params.voxelZ > 0)
            {
                VoxelSet voxels_y = SelectVoxelsOnYAxis(outer_voxels, debug_params.voxelX, debug_params.voxelZ);
                AddVoxelSetToMesh(voxels_y, half_voxel_extent * debug_params.voxelScale, out_mesh);
            }
            if (debug_params.voxelX > 0 && debug_params.voxelY > 0)
            {
                VoxelSet voxels_z = SelectVoxelsOnZAxis(outer_voxels, debug_params.voxelX, debug_params.voxelY);
                AddVoxelSetToMesh(voxels_z, half_voxel_extent * debug_params.voxelScale, out_mesh);
            }
        }
        if (debug_params.flags & DebugTypeShowInner)
        {
            for (const auto& min_distance : min_distance_field)
            {
                const u32 index = Flatten(min_distance.position, context.dimension);
                if (!voxel_field[index].inner)
                    continue;

                vec3 voxel_center = mesh_aabb.min + vec3(min_distance.position) * voxel_extent;
                if (debug_params.voxelX < 0 ||
                    debug_params.voxelY < 0 ||
                    debug_params.voxelZ < 0)
                {
                    AddVoxelToMesh(voxel_center + voxel_extent, half_voxel_extent, out_mesh);
                }
            }
        }
        if (debug_params.flags & DebugTypeShowMinDistance)
        {
            for (const auto& min_distance : min_distance_field)
            {
                vec3 voxel_center = mesh_aabb.min + vec3(min_distance.position) * voxel_extent;
                if (debug_params.voxelX == min_distance.x &&
                    debug_params.voxelY == min_distance.y &&
                    debug_params.voxelZ == min_distance.z)
                {
                    AddVoxelToMesh(voxel_center + voxel_extent, half_voxel_extent, out_mesh);

                    for (u32 x = min_distance.x; x < min_distance.x + min_distance.dist.x; ++x)
                    {
                        vec3 voxel_center = mesh_aabb.min + vec3(x, min_distance.y, min_distance.z) * voxel_extent;
                        AddVoxelToMesh(voxel_center + voxel_extent, half_voxel_extent, out_mesh);
                    }
                    for (u32 y = min_distance.y; y < min_distance.y + min_distance.dist.y; ++y)
                    {
                        vec3 voxel_center = mesh_aabb.min + vec3(min_distance.x, y, min_distance.z) * voxel_extent;
                        AddVoxelToMesh(voxel_center + voxel_extent, half_voxel_extent, out_mesh);
                    }
                    for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
                    {
                        vec3 voxel_center = mesh_aabb.min + vec3(min_distance.x, min_distance.y, z) * voxel_extent;
                        AddVoxelToMesh(voxel_center + voxel_extent, half_voxel_extent, out_mesh);
                    }
                }
            }
        }
        if (debug_params.flags & DebugTypeShowExtent)
        {
            for (const auto& min_distance : min_distance_field)
            {
                uvec3 max_extent = GetMaxAABBExtent(context, min_distance_field, voxel_field, min_distance);
                for (u32 x = min_distance.x; x < min_distance.x + max_extent.x; ++x)
                {
                    for (u32 y = min_distance.y; y < min_distance.y + max_extent.y; ++y)
                    {
                        for (u32 z = min_distance.z; z < min_distance.z + max_extent.z; ++z)
                        {
                            vec3 voxel_center = mesh_aabb.min + vec3(x, y, z) * voxel_extent;
                            AddVoxelToMesh(voxel_center + voxel_extent, half_voxel_extent, out_mesh);
                        }
                    }
                }
            }
        }
        if (debug_params.flags & DebugTypeShowResult)
        {
            for (u32 i = 0; i < max_extents.size(); ++i)
            {
                const auto& extent = max_extents[i];
                if (i == debug_params.extentIndex || debug_params.extentIndex < 0)
                {
                    vec3 half_extent = vec3(extent.extent) * voxel_extent * 0.5f;
                    vec3 aabb_center = mesh_aabb.min + vec3(extent.position) * voxel_extent + half_extent;
                    Color3u8 color = Colors[i % ARRAY_LENGTH(Colors)];
                    AddVoxelToMesh(aabb_center + half_voxel_extent, half_extent, out_mesh, color);
                }
            }
        }
    }
}

} // namespace Melt
