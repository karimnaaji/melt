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

#pragma once

#include "glm/glm.hpp"

#include <vector>

#ifndef MELT_ASSERT
#include <cassert>
#define MELT_ASSERT(stmt) assert(stmt)
#endif

namespace Melt
{

struct Mesh
{
    std::vector<glm::vec3> vertices;
    std::vector<unsigned short> indices;
};

enum OccluderBoxType
{
    OccluderBoxTypeDiagonals = 1 << 0,
    OccluderBoxTypeTop       = 1 << 1,
    OccluderBoxTypeBottom    = 1 << 2,
    OccluderBoxTypeSides     = 1 << 3,
    OccluderBoxTypeRegular   = OccluderBoxTypeSides | OccluderBoxTypeTop | OccluderBoxTypeBottom,
};

typedef int OccluderBoxTypeFlags;

struct OccluderGenerationParams
{
    OccluderBoxTypeFlags flags;

    float voxelSize;
};

enum DebugType
{
    DebugTypeShowInner          = 1 << 0,
    DebugTypeShowExtent         = 1 << 1,
    DebugTypeShowResult         = 1 << 2,
    DebugTypeShowOuter          = 1 << 3,
    DebugTypeShowMinDistance    = 1 << 4,
    DebugTypeShowSliceSelection = 1 << 5,
};

typedef int DebugTypeFlags;

struct DebugParams
{
    DebugTypeFlags flags;

    int sliceIndexX;
    int sliceIndexY;
    int sliceIndexZ;
    int voxelX;
    int voxelY;
    int voxelZ;
    int extentIndex;
    int extentMaxStep;

    float voxelScale;
};

struct OccluderBox
{
    glm::vec3 center;
    glm::vec3 halfExtent;
};

struct Result
{
    Mesh debugMesh;
    std::vector<OccluderBox> occluderBoxes;
};

Mesh GenerateConservativeOccluder(const Mesh& mesh, const OccluderGenerationParams& gen_params, const DebugParams& debug_params);

} // namespace Melt
