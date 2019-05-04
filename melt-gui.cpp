#define _USE_MATH_DEFINES
#include "glad/glad.h"
#include "GLFW/glfw3.h"
#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "imgui.cpp"
#include "imgui_draw.cpp"
#include "imgui_impl_glfw_gl3.cpp"

#define YO_IMPLEMENTATION
#define YO_NOIMG
#include "deps/yocto_obj.h"

#include "melt.h"

struct Scene
{
    GLuint program;
    GLuint vao;
    GLuint vbo;
    GLuint occ_vao;
    GLuint occ_vbo;
    GLuint occ_indices;
    GLint u_modelViewProjection;
    GLint u_alpha;
    int nvertices;
    Melt::Mesh occluder_mesh;
};

struct Camera
{
    glm::mat4 view;
    glm::vec3 position;
};

static bool LoadScene(const char* mode_path, Scene& out_scene)
{
    GLchar* vertex_source = (GLchar*)R"END(
        #version 150
        in vec3 position;
        in vec3 color;
        uniform mat4 u_modelViewProjection;
        out vec3 f_normal;
        out vec3 f_color;
        void main(void) {
            gl_Position = u_modelViewProjection * vec4(position, 1.0);
            f_color = color;
        }
    )END";

    GLchar* fragment_source = (GLchar*)R"END(
        #version 150
        in vec3 f_color;
        out vec4 color;
        uniform float u_alpha;
        void main(void) {
            color = vec4(f_color, u_alpha);
        }
    )END";

    GLint is_compiled, is_linked;

    GLuint vertex_id = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_id, 1, &vertex_source, NULL);
    glCompileShader(vertex_id);
    glGetShaderiv(vertex_id, GL_COMPILE_STATUS, &is_compiled);
    assert(is_compiled && "Vertex shader complation failed");

    GLuint fragment_id = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_id, 1, &fragment_source, NULL);
    glCompileShader(fragment_id);
    glGetShaderiv(fragment_id, GL_COMPILE_STATUS, &is_compiled);
    assert(is_compiled && "Fragment shader complation failed");

    out_scene.program = glCreateProgram();

    glAttachShader(out_scene.program, vertex_id);
    glAttachShader(out_scene.program, fragment_id);

    glBindAttribLocation(out_scene.program, 0, "position");
    glBindAttribLocation(out_scene.program, 1, "color");

    glLinkProgram(out_scene.program);
    glGetProgramiv(out_scene.program, GL_LINK_STATUS, &is_linked);
    assert(is_linked && "Program link failed");

    glDeleteShader(vertex_id);
    glDeleteShader(fragment_id);

    out_scene.u_modelViewProjection = glGetUniformLocation(out_scene.program, "u_modelViewProjection");
    out_scene.u_alpha = glGetUniformLocation(out_scene.program, "u_alpha");

    yo_scene* yo = yo_load_obj(mode_path, true, false);
    if (!yo || !yo->nshapes)
    {
        fprintf(stderr, "Error loading obj file\n");
        return false;
    }

    out_scene.nvertices = 0;
    for (int i = 0; i < yo->nshapes; i++)
        out_scene.nvertices += yo->shapes[i].nelems * 3;

    const size_t vertex_size = sizeof(glm::vec3) * 2;
    glm::vec3* buffer_data = (glm::vec3*)calloc(out_scene.nvertices, vertex_size);

    int output_index = 0;
    int input_index = 0;
    for (int i = 0; i < yo->nshapes; i++)
    {
        yo_shape* shape = yo->shapes + i;
        for (int j = 0; j < shape->nelems * 3; j++)
        {
            glm::vec3 position = *(glm::vec3*)&shape->pos[shape->elem[j] * 3];
            buffer_data[output_index++] = position;
            buffer_data[output_index++] = glm::vec3(1.0f, 0.5f, 0.5f);
            out_scene.occluder_mesh.vertices.push_back(position);
            out_scene.occluder_mesh.indices.push_back(input_index + j);
        }
        input_index += shape->nelems * 3;
    }
    yo_free_scene(yo);

    glGenVertexArrays(1, &out_scene.vao);
    glBindVertexArray(out_scene.vao);
    glGenBuffers(1, &out_scene.vbo);
    glBindBuffer(GL_ARRAY_BUFFER, out_scene.vbo);

    glBufferData(GL_ARRAY_BUFFER, out_scene.nvertices * vertex_size, buffer_data, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (void*)(sizeof(float) * 3));

    free(buffer_data);

    return true;
}

static Camera FPSCameraViewMatrix(GLFWwindow* window, bool ignore_input)
{
    static glm::vec3 position(-4.0f, 5.0f, 1.0f);
    static float rotation[] = { 0.0f, 0.0f };
    static double lastMouse[] = { 0.0, 0.0 };
    double mouse[2];

    glfwGetCursorPos(window, &mouse[0], &mouse[1]);
    if (!ignore_input && glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
    {
        rotation[0] += (float)(mouse[1] - lastMouse[1]) * 0.005f;
        rotation[1] += (float)(mouse[0] - lastMouse[0]) * 0.005f;
    }
    lastMouse[0] = mouse[0];
    lastMouse[1] = mouse[1];

    float pitch = -rotation[0];
    float yaw = rotation[1];

    if (pitch >= (float) M_PI_2) pitch = (float) M_PI_2;
    if (pitch <= (float)-M_PI_2) pitch = (float)-M_PI_2;

    glm::vec3 forward;
    const glm::vec3 up = glm::vec3(0.0, 1.0, 0.0);
    forward.x = cos(pitch) * cos(yaw);
    forward.y = sin(pitch);
    forward.z = cos(pitch) * sin(yaw);
    forward = glm::normalize(forward);

    float speed = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) ? 0.1f : 0.01f;
    if (!ignore_input)
    {
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            position -= forward * speed;
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            position += forward * speed;
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            position -= glm::normalize(glm::cross(forward, up)) * speed;
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            position += glm::normalize(glm::cross(forward, up)) * speed;
    }

    glm::mat4 view = glm::lookAt(position, position + forward, up);

    return {view, position};
}

static void ErrorCallback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

static void SetupIMGUIStyle()
{
    ImGuiStyle& style = ImGui::GetStyle();

    style.WindowMinSize     = ImVec2(320, 5000);
    style.FramePadding      = ImVec2(6, 6);
    style.ItemSpacing       = ImVec2(6, 6);
    style.ItemInnerSpacing  = ImVec2(6, 6);
    style.Alpha             = 1.0f;
    style.WindowRounding    = 0.0f;
    style.FrameRounding     = 0.0f;
    style.IndentSpacing     = 6.0f;
    style.ItemInnerSpacing  = ImVec2(6, 6);
    style.ColumnsMinSpacing = 50.0f;
    style.GrabMinSize       = 14.0f;
    style.GrabRounding      = 0.0f;
    style.ScrollbarSize     = 12.0f;
    style.ScrollbarRounding = 0.0f;

    style.Colors[ImGuiCol_Text]                  = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
    style.Colors[ImGuiCol_TextDisabled]          = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
    style.Colors[ImGuiCol_WindowBg]              = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
    style.Colors[ImGuiCol_ChildWindowBg]         = ImVec4(0.20f, 0.20f, 0.20f, 0.58f);
    style.Colors[ImGuiCol_Border]                = ImVec4(0.31f, 0.31f, 0.31f, 0.00f);
    style.Colors[ImGuiCol_BorderShadow]          = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    style.Colors[ImGuiCol_FrameBg]               = ImVec4(0.20f, 0.20f, 0.20f, 0.60f);
    style.Colors[ImGuiCol_FrameBgHovered]        = ImVec4(0.92f, 0.18f, 0.29f, 0.78f);
    style.Colors[ImGuiCol_FrameBgActive]         = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_TitleBg]               = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    style.Colors[ImGuiCol_TitleBgCollapsed]      = ImVec4(0.20f, 0.22f, 0.27f, 0.75f);
    style.Colors[ImGuiCol_TitleBgActive]         = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_MenuBarBg]             = ImVec4(0.20f, 0.22f, 0.27f, 0.47f);
    style.Colors[ImGuiCol_ScrollbarBg]           = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarGrab]         = ImVec4(0.47f, 0.47f, 0.47f, 0.21f);
    style.Colors[ImGuiCol_ScrollbarGrabHovered]  = ImVec4(0.92f, 0.18f, 0.29f, 0.78f);
    style.Colors[ImGuiCol_ScrollbarGrabActive]   = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_ComboBg]               = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
    style.Colors[ImGuiCol_CheckMark]             = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
    style.Colors[ImGuiCol_SliderGrab]            = ImVec4(0.47f, 0.47f, 0.47f, 0.14f);
    style.Colors[ImGuiCol_SliderGrabActive]      = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_Button]                = ImVec4(0.47f, 0.47f, 0.47f, 0.14f);
    style.Colors[ImGuiCol_ButtonHovered]         = ImVec4(0.92f, 0.18f, 0.29f, 0.86f);
    style.Colors[ImGuiCol_ButtonActive]          = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_Header]                = ImVec4(0.92f, 0.18f, 0.29f, 0.76f);
    style.Colors[ImGuiCol_HeaderHovered]         = ImVec4(0.92f, 0.18f, 0.29f, 0.86f);
    style.Colors[ImGuiCol_HeaderActive]          = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_Column]                = ImVec4(0.47f, 0.77f, 0.83f, 0.32f);
    style.Colors[ImGuiCol_ColumnHovered]         = ImVec4(0.92f, 0.18f, 0.29f, 0.78f);
    style.Colors[ImGuiCol_ColumnActive]          = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_ResizeGrip]            = ImVec4(0.47f, 0.77f, 0.83f, 0.04f);
    style.Colors[ImGuiCol_ResizeGripHovered]     = ImVec4(0.92f, 0.18f, 0.29f, 0.78f);
    style.Colors[ImGuiCol_ResizeGripActive]      = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_CloseButton]           = ImVec4(0.86f, 0.93f, 0.89f, 0.16f);
    style.Colors[ImGuiCol_CloseButtonHovered]    = ImVec4(0.86f, 0.93f, 0.89f, 0.39f);
    style.Colors[ImGuiCol_CloseButtonActive]     = ImVec4(0.86f, 0.93f, 0.89f, 1.00f);
    style.Colors[ImGuiCol_PlotLines]             = ImVec4(0.86f, 0.93f, 0.89f, 0.63f);
    style.Colors[ImGuiCol_PlotLinesHovered]      = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_PlotHistogram]         = ImVec4(0.86f, 0.86f, 0.86f, 0.63f);
    style.Colors[ImGuiCol_PlotHistogramHovered]  = ImVec4(0.92f, 0.18f, 0.29f, 1.00f);
    style.Colors[ImGuiCol_TextSelectedBg]        = ImVec4(0.92f, 0.18f, 0.29f, 0.43f);
    style.Colors[ImGuiCol_ModalWindowDarkening]  = ImVec4(0.20f, 0.22f, 0.27f, 0.73f);
    style.Colors[ImGuiCol_PopupBg]               = ImVec4(0.20f, 0.20f, 0.20f, 0.50f);
    ImGuiIO& io = ImGui::GetIO();

#if __APPLE__
    io.Fonts->AddFontFromFileTTF("/Library/Fonts/Arial Narrow Bold.ttf", 16);
#endif
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Please provide an obj model path\n");
        return 1;
    }

    glfwSetErrorCallback(ErrorCallback);

    if (!glfwInit())
    {
        fprintf(stderr, "glfwInit failed\n");
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    GLFWwindow *window = glfwCreateWindow(920, 720, "", NULL, NULL);
    if (!window)
    {
        fprintf(stderr, "Window creation failed\n");
        return 1;
    }

    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    glfwSwapInterval(1);

    ImGui_ImplGlfwGL3_Init(window, true);

    Scene scene;
    bool success = LoadScene(argv[1], scene);
    if (!success)
    {
        fprintf(stderr, "Failed to load scene\n");
        return 1;
    }

    Melt::DebugParams debug_params;
    std::memset(&debug_params, 0, sizeof(Melt::DebugParams));
    debug_params.voxelScale = 0.8f;

    Melt::OccluderGenerationParams gen_params;
    gen_params.voxelSize = 0.15f;
    gen_params.fillPercentage = 0.5f;

    Melt::Mesh occluder = Melt::GenerateConservativeOccluder(scene.occluder_mesh, gen_params, debug_params);

    glGenVertexArrays(1, &scene.occ_vao);
    glBindVertexArray(scene.occ_vao);

    glGenBuffers(1, &scene.occ_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, scene.occ_vbo);
    glBufferData(GL_ARRAY_BUFFER, occluder.vertices.size() * sizeof(glm::vec3), occluder.vertices.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &scene.occ_indices);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, scene.occ_indices);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, occluder.indices.size() * sizeof(GLushort), occluder.indices.data(), GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (void*)(sizeof(float) * 3));

    SetupIMGUIStyle();

    glClearColor(0.64f, 0.76f, 0.91f , 1.0f);
    glDepthFunc(GL_LEQUAL);

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        ImGui_ImplGlfwGL3_NewFrame();

        bool reload = false;
        static float alpha = 0.25f;
        static bool depth_test = false;
        static bool show_slice_selection = false;
        static bool show_inner = true;
        static bool show_outer = true;
        static bool show_dist = false;
        static bool show_extent = false;
        static bool show_result = true;
        static bool show_debug_gui = false;
        {
            debug_params.flags = 0;

            ImGui::SetNextWindowPos(ImVec2(0, 0));

            ImGuiWindowFlags options =
                ImGuiWindowFlags_NoTitleBar |
                ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoSavedSettings;

            ImGui::Begin("Fixed Overlay", nullptr, ImVec2(0,0), 0.3f, options);

            ImGui::Checkbox("Show Debug Controls", &show_debug_gui);
            ImGui::InputFloat("Voxel Size", &gen_params.voxelSize);
            ImGui::DragFloat("Fill Percentage", &gen_params.fillPercentage, 0.01f, 0.0f, 1.0f);

            if (ImGui::Button("Generate"))
            {
                reload = true;
            }

            if (show_debug_gui)
            {
                ImGui::InputFloat("Voxel Scale", &debug_params.voxelScale);
                ImGui::DragFloat("Alpha", &alpha, 0.01f, 0.0f, 1.0f);

                ImGui::InputInt("Slice X", &debug_params.sliceIndexX);
                ImGui::InputInt("Slice Y", &debug_params.sliceIndexY);
                ImGui::InputInt("Slice Z", &debug_params.sliceIndexZ);
                ImGui::InputInt("Voxel X", &debug_params.voxelX);
                ImGui::InputInt("Voxel Y", &debug_params.voxelY);
                ImGui::InputInt("Voxel Z", &debug_params.voxelZ);
                ImGui::InputInt("Extent Index", &debug_params.extentIndex);
                ImGui::InputInt("Extent Max Step", &debug_params.extentMaxStep);
                ImGui::Checkbox("Show Slice Selection", &show_slice_selection);
                ImGui::Checkbox("Show Inner", &show_inner);
                ImGui::Checkbox("Show Outer", &show_outer);
                ImGui::Checkbox("Show Dist", &show_dist);
                ImGui::Checkbox("Show Extent", &show_extent);
                ImGui::Checkbox("Show Result", &show_result);

                if (show_inner) debug_params.flags |= Melt::DebugTypeShowInner;
                if (show_slice_selection) debug_params.flags |= Melt::DebugTypeShowSliceSelection;
                if (show_outer) debug_params.flags |= Melt::DebugTypeShowOuter;
                if (show_dist) debug_params.flags |= Melt::DebugTypeShowMinDistance;
                if (show_extent) debug_params.flags |= Melt::DebugTypeShowExtent;
                if (show_result) debug_params.flags |= Melt::DebugTypeShowResult;

                ImGui::Checkbox("Depth Test", &depth_test);
                if (ImGui::Button("Next Diagonal"))
                {
                    debug_params.voxelX++;
                    debug_params.voxelY++;
                    debug_params.voxelZ++;
                }
                if (ImGui::Button("Previous Diagonal"))
                {
                    debug_params.voxelX--;
                    debug_params.voxelY--;
                    debug_params.voxelZ--;
                }
            }
            else
            {
                debug_params.extentIndex = -1;
                debug_params.flags |= Melt::DebugTypeShowResult;
            }
            ImGui::End();
        }

        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);

        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

        {
            if (depth_test)
            {
                glDisable(GL_BLEND);
                glEnable(GL_DEPTH_TEST);
            }
            else
            {
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                glDisable(GL_DEPTH_TEST);
            }

            Camera camera = FPSCameraViewMatrix(window, ImGui::IsAnyItemActive());

            float fov = 55.0f;
            glm::mat4 projection = glm::perspective(glm::radians(fov), float(width) / height, 0.01f, 100.0f);
            glm::mat4 viewProjection = projection * camera.view;

            glUseProgram(scene.program);

            glUniformMatrix4fv(scene.u_modelViewProjection, 1, GL_FALSE, glm::value_ptr(viewProjection));
            glUniform1f(scene.u_alpha, alpha);

            glBindVertexArray(scene.vao);
            glDrawArrays(GL_TRIANGLES, 0, scene.nvertices);

            glBindVertexArray(scene.occ_vao);
            glDrawElements(GL_TRIANGLES, occluder.indices.size(), GL_UNSIGNED_SHORT, 0);
        }

        if (reload)
        {
            occluder = Melt::GenerateConservativeOccluder(scene.occluder_mesh, gen_params, debug_params);

            glBindBuffer(GL_ARRAY_BUFFER, scene.occ_vbo);
            glBufferData(GL_ARRAY_BUFFER, occluder.vertices.size() * sizeof(glm::vec3), occluder.vertices.data(), GL_STATIC_DRAW);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, scene.occ_indices);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, occluder.indices.size() * sizeof(GLushort), occluder.indices.data(), GL_STATIC_DRAW);
        }

        ImGui::Render();

        glfwSwapBuffers(window);
    }

    ImGui_ImplGlfwGL3_Shutdown();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
