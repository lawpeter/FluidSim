#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/random.hpp>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw_gl3.h>
#include <limits.h>
#include <iostream>
#include <vector>
#include "particle.hpp"

#include "shaders/fragment.glsl"
#include "shaders/vertex.glsl"
#include "shaders/updateParticleCompute.glsl"
#include "shaders/gravityCompute.glsl"
#include "shaders/calculateDensitiesCompute.glsl"
#include "shaders/calculatePressureForceCompute.glsl"
#include "shaders/calculateViscosityCompute.glsl"

#define gridWidth ceil(SCREEN_WIDTH/smoothingRadius)
#define gridHeight ceil(SCREEN_HEIGHT/smoothingRadius)

int timeCount = 0;


const float smoothingRadius = 4.0f * Particle::radius;
float targetDensity = 3.75f;
float pressureMultiplier = 15.0f;
float gravity = 98.1f;
float Particle::restitution = 0.9f;
float mass = 1.0f;
float viscosityStrength = 100.0f;
float mouseInteractionRadius = 200.0f;
float mouseForceStrength = 500.0f;

bool spacePressed = false;
bool periodPressed = false;
bool periodHeld = false;
bool commaPressed = false;
bool rPressed = false;
int frameCount = 0;

float cursorX;
float cursorY;
bool leftMousePressed = false;
bool rightMousePressed = false;

int numRows = 30;
int numColumns = 30;

GLFWwindow* loadSim()
{
    if (!glfwInit())
    {
        std::cout << "failed to initialize glfw" << std::endl;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  
    GLFWwindow* window;
    window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "FluidSim", NULL, NULL);

    if (window == NULL)
    {
        std::cout << "could not create window" << std::endl;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "failed to load glad" << std::endl;
    }

    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);

    return window;
}

GLuint createShaderProgram(const char* vertexSource, const char* fragmentSource)
{
    //compile shaders
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    glCompileShader(vertexShader);

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
    glCompileShader(fragmentShader);

    //create program and attach
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);

    //link
    glLinkProgram(shaderProgram);

    //cleanup shaders
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return shaderProgram;
}

void initializeParticles(std::vector<Particle>& particles)
{
    for (int i = 0; i < numColumns; i++)
    {
        for (int j = 0; j < numRows; j++)
        {
            Particle temp(glm::vec2(50 + (j * 2.5f * Particle::radius), (50 + i * 2.5f * Particle::radius)), glm::vec2(0, 0));
            particles.push_back(temp);
        }
    }
}

int positionToCellArrayIndex(glm::vec2 position)
{
    // makes a grid based on smoothing radius
    int cellX = (int) (position.x / smoothingRadius);
    int cellY = (int) (position.y / smoothingRadius);

    return (cellY * gridWidth) + cellX;
}

void checkNearbyParticles(std::vector<int>& results, std::vector<int>& cellOffsets, Particle& particle, std::vector<Particle>& particles, std::vector<std::vector<int>>& grid)
{
    results.clear();

    int cell = positionToCellArrayIndex(particle.position + particle.velocity * 0.0008f);

    // Check all nearby grid cells and add their particles to the results array if they are within the grid
    for (int offset : cellOffsets)
    {
        int cellToCheck = cell + offset;

        if (cellToCheck < grid.size() && cellToCheck >= 0)
        {
            for (int otherParticleIndex : grid[cellToCheck])
            {
                results.push_back(otherParticleIndex);
            }
        }
    }
}

void updateParticleCells(std::vector<Particle>& particles, std::vector<std::vector<int>>& grid)
{
    grid.clear();
    grid.resize(gridWidth * gridHeight);

    for (int i = 0; i < particles.size(); i++)
    {
        // Calculate one dimensional array index from position, and add particle to respective grid cell if it is within bounds
        int cellIndex = positionToCellArrayIndex(particles[i].position + particles[i].velocity * 0.0008f);
        if (cellIndex < grid.size() && cellIndex >= 0) 
        {
            grid[cellIndex].push_back(i);
        }
        else
        {
            // std::cout << "ROGUE PARTICLE" << std::endl;
        }
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
    {
        spacePressed = !spacePressed;
    }
    if (key == GLFW_KEY_PERIOD && action == GLFW_PRESS)
    {
        periodPressed = true;
    }
    if (key == GLFW_KEY_PERIOD && action == GLFW_REPEAT)
    {
        periodHeld = true;
    }
    if (key == GLFW_KEY_PERIOD && action == GLFW_RELEASE)
    {
        periodHeld = false;
    }
    if (key == GLFW_KEY_COMMA && action == GLFW_PRESS)
    {
        commaPressed = true;
    }
    if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
        rPressed = true;
    }
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    cursorX = xpos;
    cursorY = ypos;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if (action == GLFW_PRESS)
            leftMousePressed = true;
        else if (action == GLFW_RELEASE)
            leftMousePressed = false;
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
        if (action == GLFW_PRESS)
        {
            rightMousePressed = true;
        }
        else if (action == GLFW_RELEASE)
            rightMousePressed = false;
    }
}

GLuint compileComputeShader(const char* computeShaderSource)
{
    GLuint computeShader = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(computeShader, 1, &computeShaderSource, NULL);
    glCompileShader(computeShader);

    int result;
    glGetShaderiv(computeShader, GL_COMPILE_STATUS, &result);
    if (result == GL_FALSE)
    {
        int length;
        glGetShaderiv(computeShader, GL_INFO_LOG_LENGTH, &length);
        char message[length];
        glGetShaderInfoLog(computeShader, length, &length, message);
        std::cout << message << std::endl;
    }

    //create program and attach
    GLuint computeShaderProgram = glCreateProgram();
    glAttachShader(computeShaderProgram, computeShader);
    //link
    glLinkProgram(computeShaderProgram);

    glGetProgramiv(computeShaderProgram, GL_LINK_STATUS, &result);
    if (result == GL_FALSE)
    {
        int length;
        glGetProgramiv(computeShaderProgram, GL_INFO_LOG_LENGTH, &length);
        char message[length];
        glGetProgramInfoLog(computeShaderProgram, length, &length, message);
        std::cout << message << std::endl;
    }

    glDeleteShader(computeShader);

    return computeShaderProgram;
}

void useComputeProgram(GLuint computeShaderProgram, int numParticles, GLuint particleSSBO)
{
    glUseProgram(computeShaderProgram);
    glDispatchCompute((numParticles + 127) / 128, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT); // Make sure GPU writes are visible
}

void constructOneDimensionalGrid(std::vector<std::vector<int>>& grid, int* oneDimensionalGrid, int* startingIndices, int* endIndices)
{
    int count = 0;

    for (int i = 0; i < grid.size(); i++)
    {
        startingIndices[i] = count;

        //if (grid[i].size() == 0) startingIndices[i] = -1;

        for (int j = 0; j < grid[i].size(); j++)
        {
            oneDimensionalGrid[count] = grid[i][j];
            count++;
        }

        endIndices[i] = count;
    }

    //if (count < numColumns * numRows) oneDimensionalGrid[count] = -1;

    /*
    if (timeCount > 220)
    {
        for (int i = 0; i < count; i++) 
        {
            std::cout << "Start of cell " << i << ": " << startingIndices[i] << ", End of cell " << i << ": " << endIndices[i] << std::endl;
        }
    }
    */
}

int main(int argc, char* argv[])
{
    // Orthographic projection for shaders
    glm::mat4 projection = glm::ortho(0.0f, SCREEN_WIDTH, 0.0f, SCREEN_HEIGHT, -1.0f, 1.0f);

    // Initialize window to render
    GLFWwindow* window = loadSim();

    if (window == NULL)
    {
        return -1;
    }

    // Initialize settings control window
    ImGui::CreateContext();
    ImGui_ImplGlfwGL3_Init(window, true);
    ImGui::StyleColorsDark();

    GLuint VAO;
    GLuint quadVBO;
    GLuint indexBuffer;
    GLuint particleSSBO;
    GLuint oneDSSBO;
    GLuint startIndicesSSBO;
    GLuint endIndicesSSBO;

    // Call method that creates and compiles shaders with shader sources that are #included from external files
    GLuint renderingShaderProgram = createShaderProgram(vertexShaderSource, fragmentShaderSource);
    
    //compile compute shaders
    GLuint updateParticleCompute = compileComputeShader(updateParticleComputeShaderSource);
    GLuint gravityCompute = compileComputeShader(gravityComputeShaderSource);
    GLuint densitiesCompute = compileComputeShader(calculateDensitiesComputeShaderSource);
    GLuint pressureCompute = compileComputeShader(pressureForceComputeShaderSource);
    GLuint viscosityCompute = compileComputeShader(viscosityForceComputeShaderSource);


    // Activate the shader program
    glUseProgram(renderingShaderProgram);

    // Place projection matrix inside shader uniform
    GLint projectionUniformLocation = glGetUniformLocation(renderingShaderProgram, "projection");
    glUniformMatrix4fv(projectionUniformLocation, 1, GL_FALSE, &projection[0][0]);

    // Place radius inside shader uniform
    GLint radiusUniformLocation = glGetUniformLocation(renderingShaderProgram, "radius");
    glUniform1f(radiusUniformLocation, Particle::radius);

    // Construct the corners of our quad used for rendering
    float quadVertices[8] = {
        -Particle::radius, -Particle::radius, // 0
         Particle::radius, -Particle::radius, // 1
         Particle::radius,  Particle::radius, // 2
        -Particle::radius,  Particle::radius  // 3
    };

    // Construct list of index order used to decide triangles for quad rendering
    GLuint indices[6] = {
        0, 1, 2,
        2, 3, 0
    };

    std::vector<Particle> particles;
    initializeParticles(particles);

    // Pre-size grid array to the max size
    std::vector<std::vector<int>> grid;
    grid.resize(gridWidth * gridHeight);

    //std::vector<int> cellOffsets;
    //constructCellOffsets(cellOffsets);

    int cellOffsets[9] = {-1 * gridWidth - 1, -gridWidth, -gridWidth + 1, -1, 0, 1, gridWidth - 1, gridWidth, gridWidth + 1};

    std::vector<int> results;

    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glGenBuffers(1, &quadVBO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, 8 * sizeof(float), quadVertices, GL_STATIC_DRAW);

    glGenBuffers(1, &indexBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * sizeof(GLuint), indices, GL_STATIC_DRAW);

    glGenBuffers(1, &particleSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, particleSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, particles.size() * sizeof(Particle), particles.data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particleSSBO);

    glGenBuffers(1, &oneDSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, oneDSSBO);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, oneDSSBO);

    glGenBuffers(1, &startIndicesSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, startIndicesSSBO);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, startIndicesSSBO);
    
    glGenBuffers(1, &endIndicesSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, endIndicesSSBO);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, endIndicesSSBO);

    // Assign quad position instructions for vertex shader at layout position 0
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Initialize variables used for deltaTime and average FPS
    double lastTime = glfwGetTime();
    timeCount = 0;
    float deltaTimeSum = 0.0f;

    // Put in positions to predicted particle array (used for calculations)
    std::vector<Particle> predictedParticles;
    for (Particle p : particles)
    {
        predictedParticles.push_back(p);
    }

    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Runs every frame until it is closed
    while (!glfwWindowShouldClose(window))
    {
        // Updates deltaTime and average FPS
        double nowTime = glfwGetTime();
        float deltaTime = (float) nowTime - lastTime;
        lastTime = nowTime;
        //float deltaTime = 1 / 240.0f;

        deltaTimeSum += deltaTime;
        timeCount++;

        // Keypress handling
        if (commaPressed) frameCount++;
        if (spacePressed)
        {
            deltaTime *= 0;

            if (periodHeld)
            {
                deltaTime = 1 / 240.0f;
            }
            else if (periodPressed)
            {
                deltaTime = 1 / 240.0f;
                periodPressed = false;
            }
        }
        if (commaPressed)
        {
            if (frameCount < 10)
                spacePressed = false;
            else
            {
                spacePressed = true;
                frameCount = 0;
                commaPressed = false;
            }
        }
        if (rPressed)
        {
            particles.clear();
            initializeParticles(particles);
            rPressed = false;
        }

        // Clear screen each frame
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glGetNamedBufferSubData(particleSSBO, 0, sizeof(Particle) * particles.size(), particles.data());

        ImGui_ImplGlfwGL3_NewFrame();
        
        // Repopulate grid cells with particles
        updateParticleCells(particles, grid);

        int oneDimensionalGrid[particles.size()];
        int startIndices[grid.size()];
        int endIndices[grid.size()];

        constructOneDimensionalGrid(grid, oneDimensionalGrid, startIndices, endIndices);

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, oneDSSBO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, particles.size() * sizeof(int), oneDimensionalGrid, GL_STATIC_DRAW);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, startIndicesSSBO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, grid.size() * sizeof(int), startIndices, GL_STATIC_DRAW);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, endIndicesSSBO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, grid.size() * sizeof(int), endIndices, GL_STATIC_DRAW);        

        useComputeProgram(gravityCompute, particles.size(), particleSSBO);
        useComputeProgram(updateParticleCompute, particles.size(), particleSSBO);
        useComputeProgram(densitiesCompute, particles.size(), particleSSBO);
        useComputeProgram(pressureCompute, particles.size(), particleSSBO);
        useComputeProgram(viscosityCompute, particles.size(), particleSSBO);

        glUseProgram(renderingShaderProgram);
        
        // Draw particles using instanced data from quads
        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr, numRows * numColumns);

        // Set settings to be controlled using settings window
        ImGui::Begin("Controls");
        ImGui::SliderFloat("PressureMultiplier", &pressureMultiplier, 0.0f, 50.0f);
        ImGui::SliderFloat("Gravity", &gravity, 0.0f, 200.0f);
        ImGui::SliderFloat("targetDensity", &targetDensity, 0.5f, 10.0f);
        ImGui::SliderFloat("restitution", &Particle::restitution, 0.0f, 1.0f);
        ImGui::SliderFloat("viscosityStrength", &viscosityStrength, 0.0f, 100.0f);
        ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
        ImGui::End();

        ImGui::Render();
        ImGui_ImplGlfwGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
        glfwPollEvents();

    }

    std::cout << "Average FPS: " << 1 / (deltaTimeSum / timeCount) << " " << timeCount << std::endl;

    // Cleanup
    ImGui_ImplGlfwGL3_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}