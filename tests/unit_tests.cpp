#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#define MELT_DEBUG
#define MELT_ASSERT(stmt) assert(stmt)
#define MELT_IMPLEMENTATION
#include "melt.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

static bool LoadModelMesh(const char* model_path, MeltParams& melt_params)
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

TEST_CASE("melt.flatten3d", "") 
{
#define EXPAND_DIMENSION0(DIM) \
	for (int x = 0; x < DIM.x; ++x) \
	for (int y = 0; y < DIM.y; ++y) \
	for (int z = 0; z < DIM.z; ++z) \
		REQUIRE(UnFlatten3d(Flatten3d(uvec3(x, y, z), DIM), DIM) == uvec3(x, y, z));

#define EXPAND_DIMENSION1(DIM) \
	for (int i = 0; i < DIM.x * DIM.y * DIM.z; ++i) \
		REQUIRE(Flatten3d(UnFlatten3d(i, DIM), DIM) == i);

	EXPAND_DIMENSION0(uvec3(10, 10, 10));
	EXPAND_DIMENSION0(uvec3(1, 7, 36));
	EXPAND_DIMENSION0(uvec3(56, 43, 36));

	EXPAND_DIMENSION1(uvec3(10, 10, 10));
	EXPAND_DIMENSION1(uvec3(1, 7, 36));
	EXPAND_DIMENSION1(uvec3(56, 43, 36));
}

TEST_CASE("melt.bunny", "") 
{
	MeltParams params;
    params.voxelSize = 0.25f;
    params.fillPercentage = 1.0f;

    MeltResult result;

	REQUIRE(LoadModelMesh("models/bunny.obj", params));
	REQUIRE(MeltGenerateOccluder(params, result));

    params.voxelSize = 0.05f;
	REQUIRE(!MeltGenerateOccluder(params, result));
}

TEST_CASE("melt.suzanne", "") 
{
	MeltParams params;
    params.voxelSize = 0.25f;
    params.fillPercentage = 1.0f;

    MeltResult result;

	REQUIRE(LoadModelMesh("models/suzanne.obj", params));
	REQUIRE(MeltGenerateOccluder(params, result));

    params.voxelSize = 0.15f;
    REQUIRE(MeltGenerateOccluder(params, result));
}

TEST_CASE("melt.cube", "") 
{
	MeltParams params;
    params.voxelSize = 0.25f;
    params.fillPercentage = 1.0f;

    MeltResult result;

	REQUIRE(LoadModelMesh("models/cube.obj", params));
	REQUIRE(MeltGenerateOccluder(params, result));
}

TEST_CASE("melt.sphere", "") 
{
	MeltParams params;
    params.voxelSize = 0.25f;
    params.fillPercentage = 1.0f;

    MeltResult result;

	REQUIRE(LoadModelMesh("models/sphere.obj", params));
	REQUIRE(MeltGenerateOccluder(params, result));
}

TEST_CASE("melt.teapot", "") 
{
	MeltParams params;
    params.voxelSize = 0.25f;
    params.fillPercentage = 1.0f;

    MeltResult result;

	REQUIRE(LoadModelMesh("models/teapot.obj", params));
	REQUIRE(!MeltGenerateOccluder(params, result));

    params.voxelSize = 0.5f;
    REQUIRE(MeltGenerateOccluder(params, result));
}

TEST_CASE("melt.column", "") 
{
	MeltParams params;
    params.voxelSize = 0.25f;
    params.fillPercentage = 1.0f;

    MeltResult result;

	REQUIRE(LoadModelMesh("models/column.obj", params));
	REQUIRE(MeltGenerateOccluder(params, result));
}