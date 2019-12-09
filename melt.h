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

typedef struct
{
    float x;
    float y;
    float z;
} melt_vec3_t;

typedef struct
{
    std::vector<melt_vec3_t> vertices;
    std::vector<unsigned short> indices;
}  melt_mesh_t;

typedef enum melt_occluder_box_type_t
{
    MELT_OCCLUDER_BOX_TYPE_NONE      = 0,
    MELT_OCCLUDER_BOX_TYPE_DIAGONALS = 1 << 0,
    MELT_OCCLUDER_BOX_TYPE_TOP       = 1 << 1,
    MELT_OCCLUDER_BOX_TYPE_BOTTOM    = 1 << 2,
    MELT_OCCLUDER_BOX_TYPE_SIDES     = 1 << 3,
    MELT_OCCLUDER_BOX_TYPE_REGULAR   = MELT_OCCLUDER_BOX_TYPE_SIDES | MELT_OCCLUDER_BOX_TYPE_TOP | MELT_OCCLUDER_BOX_TYPE_BOTTOM,
} melt_occluder_box_type_t;

typedef int32_t melt_occluder_box_type_flags_t;

typedef enum melt_debug_type_t
{
    MELT_DEBUG_TYPE_SHOW_INNER           = 1 << 0,
    MELT_DEBUG_TYPE_SHOW_EXTENT          = 1 << 1,
    MELT_DEBUG_TYPE_SHOW_RESULT          = 1 << 2,
    MELT_DEBUG_TYPE_SHOW_OUTER           = 1 << 3,
    MELT_DEBUG_TYPE_SHOW_MIN_DISTANCE    = 1 << 4,
    MELT_DEBUG_TYPE_SHOW_SLICE_SELECTION = 1 << 5,
} melt_debug_type_t;

typedef int32_t melt_debug_type_flags_t;

typedef struct
{
    uint32_t _start_canary;
    melt_debug_type_flags_t flags;
    int32_t voxel_x;
    int32_t voxel_y;
    int32_t voxel_z;
    int32_t extent_index;
    int32_t extent_max_step;
    float voxelScale;
    uint32_t _end_canary;
} melt_debug_params_t;

typedef struct
{
    melt_mesh_t mesh;
    melt_occluder_box_type_flags_t box_type_flags;
    melt_debug_params_t debug;
    float voxel_size;
    float fill_pct;
} melt_params_t;

typedef struct
{
    melt_mesh_t mesh;
    melt_mesh_t debug_mesh;
} melt_result_t;

int melt_generate_occluder(const melt_params_t& params, melt_result_t& result);

#ifndef MELT_ASSERT
#define MELT_ASSERT(stmt) (void)(stmt)
#endif
#ifndef MELT_PROFILE_BEGIN
#define MELT_PROFILE_BEGIN()
#endif
#ifndef MELT_PROFILE_END
#define MELT_PROFILE_END()
#endif
#ifndef MELT_MALLOC
#include <stdlib.h>
#define MELT_MALLOC(T, N) ((T*)malloc(N * sizeof(T)))
#define MELT_FREE(T) free(T)
#endif

#endif // MELT_H

#ifdef MELT_IMPLEMENTATION

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4244) // 'unsigned int' to 'float', possible loss of data
#pragma warning(disable:4201) // nonstandard extension used: nameless struct/union
#endif

#include <math.h>  // fabsf
#include <float.h> // FLT_MAX

#define MELT_ARRAY_LENGTH(array) ((int)(sizeof(array) / sizeof(*array)))

typedef unsigned int    u32;
typedef unsigned short  u16;
typedef unsigned char   u8;
typedef signed long     s64;
typedef signed short    s16;
typedef signed int      s32;
typedef signed char     s8;
typedef float           f32;

typedef melt_vec3_t vec3_t;
typedef struct
{
    f32 x, y;
} vec2_t;

typedef struct
{
    s32 x, y, z;
} svec3_t;

typedef struct
{
    u32 x, y;
} uvec2_t;

typedef struct
{
    u32 x, y, z;
} uvec3_t;

typedef uvec3_t color_3u8_t;

typedef struct
{
    vec3_t min;
    vec3_t max;
} _aabb_t;

typedef struct
{
    _aabb_t aabb;
    uvec3_t position;
} _voxel_t;

typedef struct
{
    vec3_t v0;
    vec3_t v1;
    vec3_t v2;
} _triangle_t;

typedef struct
{
    vec3_t normal;
    f32 distance;
} _plane_t;

typedef struct
{
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
} _min_distance_t;

typedef struct
{
    u8 visibility : 6;
    u8 clipped : 1;
    u8 inner : 1;
} _voxel_status_t;

typedef enum _visibility_t
{
    MELT_AXIS_VISIBILITY_NULL   =      0,
    MELT_AXIS_VISIBILITY_PLUS_X  = 1 << 0,
    MELT_AXIS_VISIBILITY_MINUS_X = 1 << 1,
    MELT_AXIS_VISIBILITY_PLUS_Y  = 1 << 2,
    MELT_AXIS_VISIBILITY_MINUS_Y = 1 << 3,
    MELT_AXIS_VISIBILITY_PLUS_Z  = 1 << 4,
    MELT_AXIS_VISIBILITY_MINUS_Z = 1 << 5,
    MELT_AXIS_VISIBILITY_ALL    =   0x3f,
} _visibility_t;

typedef struct
{
    uvec3_t position;
    uvec3_t extent;
    u32 volume;
} _max_extent_t;

typedef struct
{
    u32 element_count;
    u32 element_byte_size;
    void* data;
} _array_t;

typedef std::vector<_voxel_t> _voxel_set_t;
typedef std::vector<_max_extent_t> _max_extents_t;

typedef struct
{
    uvec3_t dimension;
    u32 size;

    s32* voxel_indices;    
    _voxel_status_t* voxel_field;
    _min_distance_t* min_distance_field;

    _voxel_t* voxel_set;
    u32 voxel_set_count;

    _max_extent_t* max_extents;
    u32 max_extents_count;
} _context_t;

typedef struct
{   
    std::vector<_voxel_set_t> x;
    std::vector<_voxel_set_t> y;
    std::vector<_voxel_set_t> z;
} _voxel_set_planes_t;

static const color_3u8_t _blue_violet       = { 138,  43, 226 };
static const color_3u8_t _dark_see_green    = { 143, 188, 143 };
static const color_3u8_t _light_slated_gray = { 119, 136, 153 };
static const color_3u8_t _powder_blue       = { 176, 224, 230 };
static const color_3u8_t _spring_green      = {   0, 255, 127 };
static const color_3u8_t _steel_blue        = {  70, 130, 180 };
static const color_3u8_t _teal              = {   0, 128, 128 };
static const color_3u8_t _whitesmoke        = { 245, 245, 245 };
static const color_3u8_t _floral_white      = { 255, 250, 240 };
static const color_3u8_t _light_pink        = { 255, 182, 193 };

static const color_3u8_t _colors[] =
{
    _whitesmoke,
    _steel_blue,
    _spring_green,
    _teal,
    _light_pink,
    _powder_blue,
    _light_slated_gray,
    _dark_see_green,
    _floral_white,
};

static const u16 _voxel_cube_indices[36] =
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

static const u16 _voxel_cube_indices_sides[24] =
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

static const u16 _voxel_cube_indices_diagonals[12] =
{
    0, 1, 6,
    0, 6, 7,
    4, 5, 2,
    4, 2, 3,
};

static const u16 _voxel_cube_indices_bottom[6] =
{
    1, 5, 6,
    1, 6, 2,
};

static const u16 _voxel_cube_indices_top[6] =
{
    0, 7, 4,
    0, 3, 7,
};

static const vec3_t _voxel_cube_vertices[8] =
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

static vec3_t _vec3_init(float x, float y, float z)
{
    vec3_t v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

static uvec3_t _uvec3_init(vec3_t v)
{ 
    uvec3_t out;
    out.x = (u32)v.x;
    out.y = (u32)v.y;
    out.z = (u32)v.z;
    return out; 
}

static vec3_t _vec3_init(uvec3_t v)
{
    vec3_t out;
    out.x = (float)v.x;
    out.y = (float)v.y;
    out.z = (float)v.z;
    return out;
}

static vec3_t _vec3_mul(vec3_t v, float factor)
{
    return _vec3_init(v.x * factor, v.y * factor, v.z * factor);
}

static vec3_t _vec3_mul(vec3_t a, vec3_t b) 
{
    return _vec3_init(a.x * b.x, a.y * b.y, a.z * b.z);
}

static vec3_t _vec3_div(vec3_t v, float divisor)
{
    return _vec3_mul(v, 1.0f / divisor);
}

static uvec3_t _uvec3_init(float x, float y, float z)
{
    uvec3_t v;
    v.x = (u32)x;
    v.y = (u32)y;
    v.z = (u32)z;
    return v;
}

static vec3_t _vec3_sub(vec3_t a, vec3_t b)
{
    return _vec3_init(a.x - b.x, a.y - b.y, a.z - b.z);
}

static vec3_t _vec3_add(vec3_t a, vec3_t b)
{
    return _vec3_init(a.x + b.x, a.y + b.y, a.z + b.z);
}

static vec3_t _vec3_abs(vec3_t v)
{
    return _vec3_init(fabsf(v.x), fabsf(v.y), fabsf(v.z));
}

static vec3_t vec3_cross(vec3_t a, vec3_t b)
{
    float x = a.y * b.z - a.z * b.y;
    float y = a.z * b.x - a.x * b.z;
    float z = a.x * b.y - a.y * b.x;
    return _vec3_init(x, y, z);
}

static f32 vec3_dot(vec3_t a, vec3_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static uvec2_t uvec2_new(u32 x, u32 y)
{
    uvec2_t v;
    v.x = x;
    v.y = y;
    return v;
}

static svec3_t _svec3_init(s32 x, s32 y, s32 z)
{
    svec3_t v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

static int _svec3_equals(svec3_t a, svec3_t b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

static int _uvec3_equals(uvec3_t a, uvec3_t b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;   
}

static float _float_min(float a, float b) {
    return a < b ? a : b;
}

static float _float_max(float a, float b) {
    return a > b ? a : b;
}

static u32 _u32_min(u32 a, u32 b)
{
    return a < b ? a : b;
}

static vec3_t _vec3_min(vec3_t a, vec3_t b)
{
    float x = _float_min(a.x, b.x);
    float y = _float_min(a.y, b.y);
    float z = _float_min(a.z, b.z);
    return _vec3_init(x, y, z);
}

static vec3_t _vec3_max(vec3_t a, vec3_t b)
{
    float x = _float_max(a.x, b.x);
    float y = _float_max(a.y, b.y);
    float z = _float_max(a.z, b.z);
    return _vec3_init(x, y, z);
}

static vec3_t aabb_center(_aabb_t aabb)
{
    return _vec3_mul(_vec3_add(aabb.min, aabb.max), 0.5f);
}

static bool _aabb_intersects_plane(const _plane_t plane, const vec3_t& half_aabb_dim)
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

static bool _aabb_intersects_triangle(const _triangle_t& triangle, const vec3_t& aabb_center, const vec3_t& half_aabb_dim)
{
    vec3_t v0, v1, v2;
    vec3_t e0, e1, e2;
    vec3_t edge_abs;
    f32 min, max;
    f32 p0, p1, p2;
    f32 rad;

    v0 = _vec3_sub(triangle.v0, aabb_center);
    v1 = _vec3_sub(triangle.v1, aabb_center);
    v2 = _vec3_sub(triangle.v2, aabb_center);

    e0 = _vec3_sub(v1, v0);
    e1 = _vec3_sub(v2, v1);
    e2 = _vec3_sub(v0, v2);

    edge_abs = _vec3_abs(e0);
    AXISTEST_X01(e0.z, e0.y, edge_abs.z, edge_abs.y);
    AXISTEST_Y02(e0.z, e0.x, edge_abs.z, edge_abs.x);
    AXISTEST_Z12(e0.y, e0.x, edge_abs.y, edge_abs.x);

    edge_abs = _vec3_abs(e1);
    AXISTEST_X01(e1.z, e1.y, edge_abs.z, edge_abs.y);
    AXISTEST_Y02(e1.z, e1.x, edge_abs.z, edge_abs.x);
    AXISTEST_Z0 (e1.y, e1.x, edge_abs.y, edge_abs.x);

    edge_abs = _vec3_abs(e2);
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

    _plane_t plane;
    plane.normal = vec3_cross(e0, e1);
    plane.distance = -vec3_dot(plane.normal, v0);

    if (!_aabb_intersects_plane(plane, half_aabb_dim))
        return false;

    return true;
}

static inline u32 _flatten_3d(const uvec3_t index, const uvec3_t dimension)
{
    u32 out_index = index.x + dimension.x * index.y + dimension.x * dimension.y * index.z;
    MELT_ASSERT(out_index < dimension.x * dimension.y * dimension.z);
    return out_index;
}
    
static inline u32 _flatten_2d(const uvec2_t index, const uvec2_t dimension)
{
    u32 out_index = index.x + dimension.x * index.y;
    MELT_ASSERT(out_index < dimension.x * dimension.y);
    return out_index;
}

static inline uvec3_t _unflatten_3d(u32 position, const uvec3_t& dimension)
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

static vec3_t _map_to_voxel_max_bound(const vec3_t& position, f32 voxel_size)
{
    auto map_func = [&](f32 value)
    {
        f32 sign = value < 0.0f ? -1.0f : 1.0f;
        f32 result = value + sign * voxel_size * 0.5f;
        return ceil(result / voxel_size) * voxel_size;
    };

    return _vec3_init(map_func(position.x), map_func(position.y), map_func(position.z));
}

static vec3_t _map_to_voxel_min_bound(const vec3_t& position, f32 voxel_size)
{
    auto map_func = [&](f32 value)
    {
        f32 sign = value < 0.0f ? -1.0f : 1.0f;
        f32 result = value + sign * voxel_size * 0.5f;
        return floorf(result / voxel_size) * voxel_size;
    };

    return _vec3_init(map_func(position.x), map_func(position.y), map_func(position.z));
}

static _aabb_t _generate_aabb(const _triangle_t& triangle)
{
    _aabb_t aabb;

    aabb.min = _vec3_init( FLT_MAX,  FLT_MAX,  FLT_MAX);
    aabb.max = _vec3_init(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    aabb.min = _vec3_min(aabb.min, triangle.v0);
    aabb.max = _vec3_max(aabb.max, triangle.v0);

    aabb.min = _vec3_min(aabb.min, triangle.v1);
    aabb.max = _vec3_max(aabb.max, triangle.v1);

    aabb.min = _vec3_min(aabb.min, triangle.v2);
    aabb.max = _vec3_max(aabb.max, triangle.v2);

    return aabb;
}

static _aabb_t _generate_aabb(const melt_mesh_t& mesh)
{
    _aabb_t aabb;

    aabb.min = _vec3_init( FLT_MAX,  FLT_MAX,  FLT_MAX);
    aabb.max = _vec3_init(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (u32 i = 0; i < mesh.indices.size(); ++i)
    {
        aabb.min = _vec3_min(aabb.min, mesh.vertices[mesh.indices[i]]);
        aabb.max = _vec3_max(aabb.max, mesh.vertices[mesh.indices[i]]);
    }

    return aabb;
}

static void _generate_per_plane_voxel_set(const _context_t& context, const _voxel_set_t& voxel_set, _voxel_set_planes_t& out_voxel_set_planes)
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
                uvec3_t position = _uvec3_init(x, y, z);
                s32 voxel_index = context.voxel_indices[_flatten_3d(position, context.dimension)];
                if (voxel_index != -1)
                {
                    u32 index_yz = _flatten_2d(uvec2_new(y, z), dim_yz);
                    u32 index_xz = _flatten_2d(uvec2_new(x, z), dim_xz);
                    u32 index_xy = _flatten_2d(uvec2_new(x, y), dim_xy);

                    _voxel_set_t& voxels_x_planes = out_voxel_set_planes.x[index_yz];
                    _voxel_set_t& voxels_y_planes = out_voxel_set_planes.y[index_xz];
                    _voxel_set_t& voxels_z_planes = out_voxel_set_planes.z[index_xy];

                    voxels_x_planes.push_back(voxel_set[voxel_index]);
                    voxels_y_planes.push_back(voxel_set[voxel_index]);
                    voxels_z_planes.push_back(voxel_set[voxel_index]);
                }
            }
        }
    }

    MELT_PROFILE_END();
}

static void _get_field(_context_t& context, const _voxel_set_planes_t& voxel_set_planes, u32 x, u32 y, u32 z, _min_distance_t& out_min_distance, _voxel_status_t& out_status)
{
    const svec3_t InfiniteDistance = _svec3_init(INT_MAX, INT_MAX, INT_MAX);
    const svec3_t NullDistance = _svec3_init(0, 0, 0);

    out_min_distance.dist = InfiniteDistance;
    out_min_distance.x = x;
    out_min_distance.y = y;
    out_min_distance.z = z;

    out_status.visibility = MELT_AXIS_VISIBILITY_NULL;
    out_status.clipped = false;
    out_status.inner = false;

    uvec2_t dim_yz = uvec2_new(context.dimension.y, context.dimension.z);
    u32 index_yz = _flatten_2d(uvec2_new(y, z), dim_yz);
    const auto& voxels_x_plane = voxel_set_planes.x[index_yz];
    for (const auto& voxel : voxels_x_plane)
    {
        s32 distance = voxel.position.x - x;
        if (distance > 0)
        {
            out_status.visibility |= MELT_AXIS_VISIBILITY_PLUS_X;
            out_min_distance.dist.x = _u32_min(out_min_distance.dist.x, distance);
        }
        else if (distance < 0)
            out_status.visibility |= MELT_AXIS_VISIBILITY_MINUS_X;
        else
            out_min_distance.dist.x = 0;
    }
    uvec2_t dim_xz = uvec2_new(context.dimension.x, context.dimension.z);
    u32 index_xz = _flatten_2d(uvec2_new(x, z), dim_xz);
    const auto& voxels_y_plane = voxel_set_planes.y[index_xz];
    for (const auto& voxel : voxels_y_plane)
    {
        s32 distance = voxel.position.y - y;
        if (distance > 0)
        {
            out_status.visibility |= MELT_AXIS_VISIBILITY_PLUS_Y;
            out_min_distance.dist.y = _u32_min(out_min_distance.dist.y, distance);
        }
        else if (distance < 0)
            out_status.visibility |= MELT_AXIS_VISIBILITY_MINUS_Y;
        else
            out_min_distance.dist.y = 0;
    }
    uvec2_t dim_xy = uvec2_new(context.dimension.x, context.dimension.y);
    u32 index_xy = _flatten_2d(uvec2_new(x, y), dim_xy);
    const auto& voxels_z_plane = voxel_set_planes.z[index_xy];
    for (const auto& voxel : voxels_z_plane)
    {
        s32 distance = voxel.position.z - z;
        if (distance > 0)
        {
            out_status.visibility |= MELT_AXIS_VISIBILITY_PLUS_Z;
            out_min_distance.dist.z = _u32_min(out_min_distance.dist.z, distance);
        }
        else if (distance < 0)
            out_status.visibility |= MELT_AXIS_VISIBILITY_MINUS_Z;
        else
            out_min_distance.dist.z = 0;
    }
    if (out_status.visibility == MELT_AXIS_VISIBILITY_ALL)
    {
        if (!_svec3_equals(out_min_distance.dist, InfiniteDistance) &&
            !_svec3_equals(out_min_distance.dist, NullDistance))
        {
            out_status.inner = true;
        }
    }
}

static void _generate_fields(_context_t& context, const _voxel_set_planes_t& voxel_set_planes)
{
    MELT_PROFILE_BEGIN();

    for (u32 i = 0; i < context.size; ++i)
    {
        _min_distance_t& min_distance = context.min_distance_field[i];
        _voxel_status_t& voxel_status = context.voxel_field[i];
        const uvec3_t position = _unflatten_3d(i, context.dimension);
        _get_field(context, voxel_set_planes, position.x, position.y, position.z, min_distance, voxel_status);
    }

    MELT_PROFILE_END();
}

static inline bool _inner_voxel(_voxel_status_t voxel_status)
{
    return voxel_status.inner && !voxel_status.clipped;
}

static uvec3_t _get_max_aabb_extent(const _context_t& context, const _min_distance_t& min_distance)
{
    MELT_PROFILE_BEGIN();

    std::vector<uvec2_t> max_aabb_extents;

    for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
    {
        uvec3_t z_slice_position = _uvec3_init(min_distance.x, min_distance.y, z);
        u32 z_slice_index = _flatten_3d(z_slice_position, context.dimension);

        MELT_ASSERT(context.voxel_field[z_slice_index].inner);

        if (context.voxel_field[z_slice_index].clipped)
            continue;

        const _min_distance_t& sample_min_distance = context.min_distance_field[z_slice_index];

        uvec2_t max_extent = uvec2_new(sample_min_distance.dist.x, sample_min_distance.dist.y);

        u32 x = sample_min_distance.x + 1;
        u32 y = sample_min_distance.y + 1;
        u32 i = 1;
        while (x < sample_min_distance.x + sample_min_distance.dist.x &&
               y < sample_min_distance.y + sample_min_distance.dist.y)
        {
            const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
            if (_inner_voxel(context.voxel_field[index]))
            {
                const _min_distance_t& distance = context.min_distance_field[index];
                max_extent.x = _u32_min(distance.dist.x + i, max_extent.x);
                max_extent.y = _u32_min(distance.dist.y + i, max_extent.y);
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
        min_extent.x = _u32_min(extent.x, min_extent.x);
        min_extent.y = _u32_min(extent.y, min_extent.y);

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

    return _uvec3_init(min_extent.x, min_extent.y, z_slice - 1);
}

static void _clip_voxel_field(const _context_t& context, const uvec3_t start_position, const uvec3_t extent)
{
    MELT_PROFILE_BEGIN();

    for (u32 x = start_position.x; x < start_position.x + extent.x; ++x)
    {
        for (u32 y = start_position.y; y < start_position.y + extent.y; ++y)
        {
            for (u32 z = start_position.z; z < start_position.z + extent.z; ++z)
            {
                u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
                MELT_ASSERT(!context.voxel_field[index].clipped && "Clipping already clipped voxel field index");
                context.voxel_field[index].clipped = true;
            }
        }
    }

    MELT_PROFILE_END();
}
    
static bool _water_tight_mesh(const _context_t& context)
{
    for (u32 i = 0; i < context.size; ++i)
    {
        const _min_distance_t& min_distance = context.min_distance_field[i];
        if (!_inner_voxel(context.voxel_field[_flatten_3d(min_distance.position, context.dimension)]))
            continue;
        
        for (u32 x = min_distance.x; x < min_distance.x + min_distance.dist.x; ++x)
        {
            const u32 y = min_distance.y;
            const u32 z = min_distance.z;
            const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
            if (!_inner_voxel(context.voxel_field[index]))
            {
                return false;
            }
        }
        for (u32 y = min_distance.y; y < min_distance.y + min_distance.dist.y; ++y)
        {
            const u32 x = min_distance.x;
            const u32 z = min_distance.z;
            const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
            if (!_inner_voxel(context.voxel_field[index]))
            {
                return false;
            }
        }
        for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
        {
            const u32 x = min_distance.x;
            const u32 y = min_distance.y;
            const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
            if (!_inner_voxel(context.voxel_field[index]))
            {
                return false;
            }
        }
    }
    
    return true;
}

static void _debug_validate_min_distance_field(const _context_t& context, const _voxel_set_t& outer_voxels)
{
#if defined(MELT_DEBUG) && defined(MELT_ASSERT)
    for (u32 i = 0; i < context.size; ++i)
    {
        const _min_distance_t& min_distance = context.min_distance_field[i];
        if (!_inner_voxel(context.voxel_field[_flatten_3d(min_distance.position, context.dimension)]))
            continue;

        for (u32 x = min_distance.x; x < min_distance.x + min_distance.dist.x; ++x)
        {
            const u32 y = min_distance.y;
            const u32 z = min_distance.z;
            const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
            for (const _voxel_t& voxel : outer_voxels)
                MELT_ASSERT(!_uvec3_equals(voxel.position, _uvec3_init(x, y, z)));
            MELT_ASSERT(_inner_voxel(context.voxel_field[index]));
        }
        for (u32 y = min_distance.y; y < min_distance.y + min_distance.dist.y; ++y)
        {
            const u32 x = min_distance.x;
            const u32 z = min_distance.z;
            const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
            for (const _voxel_t& voxel : outer_voxels)
                MELT_ASSERT(!_uvec3_equals(voxel.position, _uvec3_init(x, y, z)));
            MELT_ASSERT(_inner_voxel(context.voxel_field[index]));
        }
        for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
        {
            const u32 x = min_distance.x;
            const u32 y = min_distance.y;
            const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
            for (const _voxel_t& voxel : outer_voxels)
                MELT_ASSERT(!_uvec3_equals(voxel.position, _uvec3_init(x, y, z)));
            MELT_ASSERT(_inner_voxel(context.voxel_field[index]));
        }
    }
#endif
}

static void _debug_validate_max_extents(const _max_extents_t& max_extents, const _voxel_set_t& outer_voxels)
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
                    for (const _voxel_t& voxel : outer_voxels)
                    {
                        MELT_ASSERT(!_uvec3_equals(voxel.position, _uvec3_init(x, y, z)));
                    }
                }
            }
        }
    }
#endif
}

static void _update_min_distance_field(const _context_t& context, const uvec3_t& start_position, const uvec3_t& extent)
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
                const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
                if (_inner_voxel(context.voxel_field[index]))
                {
                    _min_distance_t& min_distance = context.min_distance_field[index];
                    const u32 updated_distance_x = start_position.x - min_distance.x;
                    min_distance.dist.x = _u32_min(updated_distance_x, min_distance.dist.x);
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
                const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
                if (_inner_voxel(context.voxel_field[index]))
                {
                    _min_distance_t& min_distance = context.min_distance_field[index];
                    const u32 updated_distance_y = start_position.y - min_distance.y;
                    min_distance.dist.y = _u32_min(updated_distance_y, min_distance.dist.y);
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
                const u32 index = _flatten_3d(_uvec3_init(x, y, z), context.dimension);
                if (_inner_voxel(context.voxel_field[index]))
                {
                    _min_distance_t& min_distance = context.min_distance_field[index];
                    const u32 updated_distance_z = start_position.z - min_distance.z;
                    min_distance.dist.z = _u32_min(updated_distance_z, min_distance.dist.z);
                }
            }
        }
    }

    MELT_PROFILE_END();
}

static melt_occluder_box_type_t _select_voxel_indices(melt_occluder_box_type_flags_t box_type_flags, const u16*& out_indices, u16& out_index_length)
{
#define CHECK_FLAG(flag) (box_type_flags & flag) == flag

    if (CHECK_FLAG(MELT_OCCLUDER_BOX_TYPE_REGULAR))
    {
        out_indices = static_cast<const u16*>(_voxel_cube_indices);
        out_index_length = MELT_ARRAY_LENGTH(_voxel_cube_indices);
        return MELT_OCCLUDER_BOX_TYPE_REGULAR;
    }
    else if (CHECK_FLAG(MELT_OCCLUDER_BOX_TYPE_SIDES))
    {
        out_indices = static_cast<const u16*>(_voxel_cube_indices_sides);
        out_index_length = MELT_ARRAY_LENGTH(_voxel_cube_indices_sides);
        return MELT_OCCLUDER_BOX_TYPE_SIDES;
    }
    else if (CHECK_FLAG(MELT_OCCLUDER_BOX_TYPE_BOTTOM))
    {
        out_indices = static_cast<const u16*>(_voxel_cube_indices_bottom);
        out_index_length = MELT_ARRAY_LENGTH(_voxel_cube_indices_bottom);
        return MELT_OCCLUDER_BOX_TYPE_BOTTOM;
    }
    else if (CHECK_FLAG(MELT_OCCLUDER_BOX_TYPE_TOP))
    {
        out_indices = static_cast<const u16*>(_voxel_cube_indices_top);
        out_index_length = MELT_ARRAY_LENGTH(_voxel_cube_indices_top);
        return MELT_OCCLUDER_BOX_TYPE_TOP;
    }
    else if (CHECK_FLAG(MELT_OCCLUDER_BOX_TYPE_DIAGONALS))
    {
        out_indices = static_cast<const u16*>(_voxel_cube_indices_diagonals);
        out_index_length = MELT_ARRAY_LENGTH(_voxel_cube_indices_diagonals);
        return MELT_OCCLUDER_BOX_TYPE_DIAGONALS;
    }

#undef CHECK_FLAG

    return MELT_OCCLUDER_BOX_TYPE_NONE;
}

static void _add_voxel_to_mesh(vec3_t voxel_center, vec3_t half_voxel_size, melt_mesh_t& mesh, melt_occluder_box_type_flags_t box_type_flags)
{
    u16 index_offset = (u16)mesh.vertices.size();
        
    for (u32 i = 0; i < MELT_ARRAY_LENGTH(_voxel_cube_vertices); ++i)
    {
        vec3_t vertex = _vec3_add(_vec3_mul(half_voxel_size, _voxel_cube_vertices[i]), voxel_center);
        mesh.vertices.push_back(vertex);
    }

    while (box_type_flags != MELT_OCCLUDER_BOX_TYPE_NONE)
    {
        const u16* indices = nullptr;
        u16 indices_length = 0;

        melt_occluder_box_type_t selected_type = _select_voxel_indices(box_type_flags, indices, indices_length);

        MELT_ASSERT(indices && indices_length > 0);
        for (u32 i = 0; i < indices_length; ++i)
        {
            mesh.indices.push_back(indices[i] + index_offset);
        }

        box_type_flags &= ~selected_type;
    }
}

#if defined(MELT_DEBUG)
static void _add_voxel_with_color_to_mesh(vec3_t voxel_center, vec3_t half_voxel_size, melt_mesh_t& mesh, melt_occluder_box_type_flags_t box_type_flags = MELT_OCCLUDER_BOX_TYPE_REGULAR, const color_3u8_t& color = _blue_violet)
{
    u16 index_offset = (u16)mesh.vertices.size() / 2;

    for (u32 i = 0; i < MELT_ARRAY_LENGTH(_voxel_cube_vertices); ++i)
    {
        vec3_t vertex = _vec3_add(_vec3_mul(half_voxel_size, _voxel_cube_vertices[i]), voxel_center);
        mesh.vertices.push_back(vertex);
        mesh.vertices.push_back(_vec3_div(_vec3_init(color), 255.0f));
    }

    while (box_type_flags != MELT_OCCLUDER_BOX_TYPE_NONE)
    {
        const u16* indices = nullptr;
        u16 indices_length = 0;

        melt_occluder_box_type_t selected_type = _select_voxel_indices(box_type_flags, indices, indices_length);

        MELT_ASSERT(indices && indices_length > 0);
        for (u32 i = 0; i < indices_length; ++i)
        {
            mesh.indices.push_back(indices[i] + index_offset);
        }

        box_type_flags &= ~selected_type;
    }
}

static void _add_voxel_set_to_mesh(const _voxel_set_t& voxel_set, const vec3_t& half_voxel_extent, melt_mesh_t& mesh)
{
    for (const _voxel_t& voxel : voxel_set)
    {
        _add_voxel_with_color_to_mesh(aabb_center(voxel.aabb), half_voxel_extent, mesh);
    }
}
#endif

static _max_extent_t _get_max_extent(const _context_t& context)
{
    MELT_PROFILE_BEGIN();

    _max_extent_t max_extent;
    max_extent.extent = _uvec3_init(0, 0, 0);
    max_extent.position = _uvec3_init(0, 0, 0);
    max_extent.volume = 0;

    for (u32 i = 0; i < context.size; ++i)
    {
        const _min_distance_t& min_distance = context.min_distance_field[i];
        _voxel_status_t voxel_status = context.voxel_field[_flatten_3d(min_distance.position, context.dimension)];
        if (_inner_voxel(voxel_status))
        {
            uvec3_t extent = _get_max_aabb_extent(context, min_distance);
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

int melt_generate_occluder(const melt_params_t& params, melt_result_t& out_result)
{
    out_result.debug_mesh.vertices.clear();
    out_result.debug_mesh.indices.clear();
    out_result.mesh.vertices.clear();
    out_result.mesh.indices.clear();

    vec3_t voxel_extent = _vec3_init(params.voxel_size, params.voxel_size, params.voxel_size);
    vec3_t half_voxel_extent = _vec3_mul(voxel_extent, 0.5f);

    _aabb_t mesh_aabb = _generate_aabb(params.mesh);

    mesh_aabb.min = _vec3_sub(_map_to_voxel_min_bound(mesh_aabb.min, params.voxel_size), voxel_extent);
    mesh_aabb.max = _vec3_add(_map_to_voxel_max_bound(mesh_aabb.max, params.voxel_size), voxel_extent);

    vec3_t mesh_extent = _vec3_sub(mesh_aabb.max, mesh_aabb.min);
    vec3_t inv_mesh_extent = _vec3_init(1.0f / mesh_extent.x, 1.0f / mesh_extent.y, 1.0f / mesh_extent.z);
    vec3_t voxel_count = _vec3_div(mesh_extent, params.voxel_size);
    vec3_t voxel_resolution = _vec3_mul(voxel_count, inv_mesh_extent);

    _context_t context;
    memset(&context, 0, sizeof(_context_t));
    context.dimension = _uvec3_init(voxel_count);
    context.size = u32(voxel_count.x) * u32(voxel_count.y) * u32(voxel_count.z);
    context.voxel_field = MELT_MALLOC(_voxel_status_t, context.size);
    context.min_distance_field = MELT_MALLOC(_min_distance_t, context.size);
    context.voxel_indices = MELT_MALLOC(s32, context.size);
    for (u32 i = 0; i < context.size; ++i)
        context.voxel_indices[i] = -1;

    _voxel_set_t outer_voxels;
    outer_voxels.reserve(context.size);

    // Perform shell voxelization
    for (u32 i = 0; i < params.mesh.indices.size(); i += 3)
    {
        MELT_PROFILE_BEGIN();

        _triangle_t triangle;

        triangle.v0 = params.mesh.vertices[params.mesh.indices[i + 0]];
        triangle.v1 = params.mesh.vertices[params.mesh.indices[i + 1]];
        triangle.v2 = params.mesh.vertices[params.mesh.indices[i + 2]];

        _aabb_t triangle_aabb = _generate_aabb(triangle);

        // Voxel snapping, snap the triangle extent to find the 3d grid to iterate on.
        triangle_aabb.min = _vec3_sub(_map_to_voxel_min_bound(triangle_aabb.min, params.voxel_size), voxel_extent);
        triangle_aabb.max = _vec3_add(_map_to_voxel_max_bound(triangle_aabb.max, params.voxel_size), voxel_extent);

        for (f32 x = triangle_aabb.min.x; x <= triangle_aabb.max.x; x += params.voxel_size)
        {
            for (f32 y = triangle_aabb.min.y; y <= triangle_aabb.max.y; y += params.voxel_size)
            {
                for (f32 z = triangle_aabb.min.z; z <= triangle_aabb.max.z; z += params.voxel_size)
                {
                    _voxel_t voxel;

                    voxel.aabb.min = _vec3_sub(_vec3_init(x, y, z), half_voxel_extent);
                    voxel.aabb.max = _vec3_add(_vec3_init(x, y, z), half_voxel_extent);

                    MELT_ASSERT(voxel.aabb.min.x >= mesh_aabb.min.x - half_voxel_extent.x);
                    MELT_ASSERT(voxel.aabb.min.y >= mesh_aabb.min.y - half_voxel_extent.y);
                    MELT_ASSERT(voxel.aabb.min.z >= mesh_aabb.min.z - half_voxel_extent.z);

                    MELT_ASSERT(voxel.aabb.max.x <= mesh_aabb.max.x + half_voxel_extent.x);
                    MELT_ASSERT(voxel.aabb.max.y <= mesh_aabb.max.y + half_voxel_extent.y);
                    MELT_ASSERT(voxel.aabb.max.z <= mesh_aabb.max.z + half_voxel_extent.z);

                    vec3_t voxel_center = aabb_center(voxel.aabb);
                    vec3_t relative_to_origin = _vec3_sub(_vec3_sub(voxel_center, mesh_aabb.min), half_voxel_extent);

                    if (!_aabb_intersects_triangle(triangle, voxel_center, half_voxel_extent))
                        continue;
                    voxel.position = _uvec3_init(_vec3_mul(relative_to_origin, voxel_resolution));

                    const u32 index = _flatten_3d(voxel.position, context.dimension);
                    if (context.voxel_indices[index] != -1)
                        continue;

                    context.voxel_indices[index] = (s32)outer_voxels.size();
                    outer_voxels.push_back(voxel);
                }
            }
        }

        MELT_PROFILE_END();
    }

    _voxel_set_planes_t voxel_set_planes;

    // Generate a flat voxel list per plane (x,y), (x,z), (y,z)
    _generate_per_plane_voxel_set(context, outer_voxels, voxel_set_planes);

    // The minimum distance field is a data structure representing, for each voxel,
    // the minimum distance that we can go in each of the positive directions x, y,
    // z until we collide with a shell voxel. The voxel field is a data structure
    // representing the state of the voxels. The clip status is a representation of
    // whether a voxel is in the clip state, whether a voxel can 'see' a voxel in
    // each of the directions +x,-x,+y,-y,+z,-z, and whether the voxel is an 'inner'
    // voxel (contained within the shell voxels).

    // Generate the minimum distance field, and voxel status from the initial shell.

    _generate_fields(context, voxel_set_planes);

    if (!_water_tight_mesh(context))
    {
        return 0;
    }

    _debug_validate_min_distance_field(context, outer_voxels);

    std::vector<_max_extent_t> max_extents;
    if (params.debug.extent_max_step != 0)
    {
        for (s32 i = 0; i < params.debug.extent_max_step; ++i)
        {
            _max_extent_t max_extent = _get_max_extent(context);

            _clip_voxel_field(context, max_extent.position, max_extent.extent);

            _update_min_distance_field(context, max_extent.position, max_extent.extent);

            _debug_validate_min_distance_field(context, outer_voxels);

            max_extents.push_back(max_extent);
        }
    }
    else
    {
        u32 volume = 0;
        u32 total_volume = 0;
        f32 fill_pct = 0.0f;

        // Approximate the volume of the mesh by the number of voxels that can fit within.
        for (u32 i = 0; i < context.size; ++i)
        {
            const _min_distance_t& min_distance = context.min_distance_field[i];
            _voxel_status_t voxel_status = context.voxel_field[_flatten_3d(min_distance.position, context.dimension)];

            // Each inner voxel adds one unit to the volume.
            if (_inner_voxel(voxel_status))
                ++total_volume;
        }

        // One iteration to find an extent does the following:
        // . Get the extent that maximizes the volume considering the minimum distance
        //    field
        // . Clip the max extent found to the set of inner voxels
        // . Update the minimum distance field by adjusting the distances on the set
        //    of inner voxels. This is done by extending the extent cube to infinity
        //    on each of the axes +x, +y, +z
        while (fill_pct < params.fill_pct && volume != total_volume)
        {
            _max_extent_t max_extent = _get_max_extent(context);

            _clip_voxel_field(context, max_extent.position, max_extent.extent);

            _update_min_distance_field(context, max_extent.position, max_extent.extent);

            _debug_validate_min_distance_field(context, outer_voxels);

            max_extents.push_back(max_extent);

            fill_pct += f32(max_extent.volume) / total_volume;
            volume += max_extent.volume;
        }
    }

    for (u32 i = 0; i < max_extents.size(); ++i)
    {
        const _max_extent_t& extent = max_extents[i];

        vec3_t half_extent = _vec3_mul(_vec3_init(extent.extent), half_voxel_extent);
        vec3_t voxel_position = _vec3_mul(_vec3_init(extent.position), voxel_extent);
        vec3_t voxel_position_biased_to_center = _vec3_add(voxel_position, half_extent);
        vec3_t aabb_center = _vec3_add(mesh_aabb.min, voxel_position_biased_to_center);

        _add_voxel_to_mesh(_vec3_add(aabb_center, half_voxel_extent), half_extent, out_result.mesh, params.box_type_flags);
    }

    _debug_validate_max_extents(max_extents, outer_voxels);

#if defined(MELT_DEBUG)
    MELT_ASSERT(params.debug._start_canary == 0 && params.debug._end_canary == 0 && "Make sure to memset melt_debug_params_t to 0 before use");

    if (params.debug.flags > 0)
    {
        // _add_voxel_with_color_to_mesh(aabb_center(mesh_aabb), (mesh_aabb.max - mesh_aabb.min) * 0.5f, out_result.debug_mesh, _colors[0]);

        if (params.debug.flags & MELT_DEBUG_TYPE_SHOW_OUTER)
        {
            _add_voxel_set_to_mesh(outer_voxels, _vec3_mul(half_voxel_extent, params.debug.voxelScale), out_result.debug_mesh);
        }
        if (params.debug.flags & MELT_DEBUG_TYPE_SHOW_SLICE_SELECTION)
        {
            if (params.debug.voxel_y > 0 && params.debug.voxel_z > 0)
            {
                u32 index = _flatten_2d(uvec2_new(params.debug.voxel_y, params.debug.voxel_z), uvec2_new(context.dimension.y, context.dimension.z));
                _voxel_set_t& voxels_x_planes = voxel_set_planes.x[index];
                _add_voxel_set_to_mesh(voxels_x_planes, _vec3_mul(half_voxel_extent, params.debug.voxelScale), out_result.debug_mesh);
            }
            if (params.debug.voxel_x > 0 && params.debug.voxel_z > 0)
            {
                u32 index = _flatten_2d(uvec2_new(params.debug.voxel_x, params.debug.voxel_z), uvec2_new(context.dimension.x, context.dimension.z));
                _voxel_set_t& voxels_y_planes = voxel_set_planes.y[index];
                _add_voxel_set_to_mesh(voxels_y_planes, _vec3_mul(half_voxel_extent, params.debug.voxelScale), out_result.debug_mesh);
            }
            if (params.debug.voxel_x > 0 && params.debug.voxel_y > 0)
            {
                u32 index = _flatten_2d(uvec2_new(params.debug.voxel_x, params.debug.voxel_y), uvec2_new(context.dimension.x, context.dimension.y));
                _voxel_set_t& voxels_z_planes = voxel_set_planes.z[index];
                _add_voxel_set_to_mesh(voxels_z_planes, _vec3_mul(half_voxel_extent, params.debug.voxelScale), out_result.debug_mesh);
            }
        }
        if (params.debug.flags & MELT_DEBUG_TYPE_SHOW_INNER)
        {
            for (u32 i = 0; i < context.size; ++i)
            {
                const _min_distance_t& min_distance = context.min_distance_field[i];
                const u32 index = _flatten_3d(min_distance.position, context.dimension);
                if (!context.voxel_field[index].inner)
                    continue;

                vec3_t voxel_position = _vec3_mul(_vec3_init(min_distance.position), voxel_extent);
                vec3_t voxel_center = _vec3_add(mesh_aabb.min, voxel_position);
                if (params.debug.voxel_x < 0 ||
                    params.debug.voxel_y < 0 ||
                    params.debug.voxel_z < 0)
                {
                    _add_voxel_with_color_to_mesh(_vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debug_mesh);
                }
            }
        }
        if (params.debug.flags & MELT_DEBUG_TYPE_SHOW_MIN_DISTANCE)
        {
            for (u32 i = 0; i < context.size; ++i)
            {
                const _min_distance_t& min_distance = context.min_distance_field[i];
                vec3_t voxel_center = _vec3_add(mesh_aabb.min, _vec3_mul(_vec3_init(min_distance.position), voxel_extent));
                if (u32(params.debug.voxel_x) == min_distance.x &&
                    u32(params.debug.voxel_y) == min_distance.y &&
                    u32(params.debug.voxel_z) == min_distance.z)
                {
                    _add_voxel_with_color_to_mesh(_vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debug_mesh);

                    for (u32 x = min_distance.x; x < min_distance.x + min_distance.dist.x; ++x)
                    {
                        vec3_t voxel_center_x = _vec3_add(mesh_aabb.min, _vec3_mul(_vec3_init(x, min_distance.y, min_distance.z), voxel_extent));
                        _add_voxel_with_color_to_mesh(_vec3_add(voxel_center_x, voxel_extent), half_voxel_extent, out_result.debug_mesh);
                    }
                    for (u32 y = min_distance.y; y < min_distance.y + min_distance.dist.y; ++y)
                    {
                        vec3_t voxel_center_y = _vec3_add(mesh_aabb.min, _vec3_mul(_vec3_init(min_distance.x, y, min_distance.z), voxel_extent));
                        _add_voxel_with_color_to_mesh(_vec3_add(voxel_center_y, voxel_extent), half_voxel_extent, out_result.debug_mesh);
                    }
                    for (u32 z = min_distance.z; z < min_distance.z + min_distance.dist.z; ++z)
                    {
                        vec3_t voxel_center_z = _vec3_add(mesh_aabb.min, _vec3_mul(_vec3_init(min_distance.x, min_distance.y, z), voxel_extent));
                        _add_voxel_with_color_to_mesh(_vec3_add(voxel_center_z, voxel_extent), half_voxel_extent, out_result.debug_mesh);
                    }
                }
            }
        }
        if (params.debug.flags & MELT_DEBUG_TYPE_SHOW_EXTENT)
        {
            for (u32 i = 0; i < context.size; ++i)
            {
                const _min_distance_t& min_distance = context.min_distance_field[i];
                uvec3_t max_extent = _get_max_aabb_extent(context, min_distance);
                for (u32 x = min_distance.x; x < min_distance.x + max_extent.x; ++x)
                {
                    for (u32 y = min_distance.y; y < min_distance.y + max_extent.y; ++y)
                    {
                        for (u32 z = min_distance.z; z < min_distance.z + max_extent.z; ++z)
                        {
                            vec3_t voxel_center = _vec3_add(mesh_aabb.min, _vec3_mul(_vec3_init(x, y, z), voxel_extent));
                            _add_voxel_with_color_to_mesh(_vec3_add(voxel_center, voxel_extent), half_voxel_extent, out_result.debug_mesh);
                        }
                    }
                }
            }
        }
        if (params.debug.flags & MELT_DEBUG_TYPE_SHOW_RESULT)
        {
            for (size_t i = 0; i < max_extents.size(); ++i)
            {
                const _max_extent_t& extent = max_extents[i];
                if (s32(i) == params.debug.extent_index || params.debug.extent_index < 0)
                {
                    vec3_t half_extent = _vec3_mul(_vec3_init(extent.extent), half_voxel_extent);
                    vec3_t voxel_position = _vec3_mul(_vec3_init(extent.position), voxel_extent);
                    vec3_t voxel_position_biased_to_center = _vec3_add(voxel_position, half_extent);
                    vec3_t aabb_center = _vec3_add(mesh_aabb.min, voxel_position_biased_to_center);
                    static color_3u8_t color = _colors[i % MELT_ARRAY_LENGTH(_colors)];
                    _add_voxel_with_color_to_mesh(_vec3_add(aabb_center, half_voxel_extent), half_extent, out_result.debug_mesh, params.box_type_flags, color);
                }
            }
        }
    }
#endif

    MELT_FREE(context.voxel_indices);
    MELT_FREE(context.voxel_field);
    MELT_FREE(context.min_distance_field);
    return 1;
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // MELT_IMPLEMENTATION
