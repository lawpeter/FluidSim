#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/random.hpp>
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
float Particle::restitution = 0.4f;
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

float cursorX = 0.0f;
float cursorY = 0.0f;
bool leftMousePressed = false;
bool rightMousePressed = false;

int numRows = 200;
int numColumns = 200;

float fixedDeltaTime = 0.0022222222222f;

GLuint VAO;
GLuint quadVBO;
GLuint indexBuffer;
GLuint particleSSBO;
GLuint oneDSSBO;
GLuint startIndicesSSBO;
GLuint endIndicesSSBO;
GLuint mouseinfoSSBO;

std::vector<Particle> particles;
std::vector<std::vector<int>> grid;

GLuint renderingShaderProgram;
GLuint updateParticleCompute;
GLuint gravityCompute;
GLuint densitiesCompute;
GLuint pressureCompute;
GLuint viscosityCompute;


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
            Particle temp(glm::vec2(320 + (j * 2.0f * Particle::radius), i * 2.0f * Particle::radius), glm::vec2(0, 0));
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

void updateParticleCells(std::vector<Particle>& particles, std::vector<std::vector<int>>& grid)
{
    grid.clear();
    grid.resize(gridWidth * gridHeight);

    for (int i = 0; i < particles.size(); i++)
    {
        // Calculate one dimensional array index from position, and add particle to respective grid cell if it is within bounds
        int cellIndex = positionToCellArrayIndex(particles[i].position + particles[i].velocity * fixedDeltaTime);
        if (cellIndex < grid.size() && cellIndex >= 0) 
        {
            grid[cellIndex].push_back(i);
        }
        else
        {
            //std::cout << "ROGUE PARTICLE" << std::endl;
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

void useComputeProgram(GLuint computeShaderProgram, int numParticles)
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

        for (int j = 0; j < grid[i].size(); j++)
        {
            oneDimensionalGrid[count] = grid[i][j];
            count++;
        }

        endIndices[i] = count;
    }
}

void runFrame()
{
    // Clear screen each frame
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glGetNamedBufferSubData(particleSSBO, 0, sizeof(Particle) * particles.size(), particles.data());
    
    // Repopulate grid cells with particles
    updateParticleCells(particles, grid);

    int oneDimensionalGrid[particles.size()];
    int startIndices[grid.size()];
    int endIndices[grid.size()];

    constructOneDimensionalGrid(grid, oneDimensionalGrid, startIndices, endIndices);

    float mouseinfo[4] = {(leftMousePressed == 1 ? 0.0f : 1.0f), (rightMousePressed == 1 ? 0.0f : 1.0f), cursorX, cursorY};

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, oneDSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, particles.size() * sizeof(int), oneDimensionalGrid, GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, startIndicesSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, grid.size() * sizeof(int), startIndices, GL_STATIC_DRAW);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, endIndicesSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, grid.size() * sizeof(int), endIndices, GL_STATIC_DRAW);   
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, mouseinfoSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, 4 * sizeof(float), mouseinfo, GL_STATIC_DRAW); 

    useComputeProgram(gravityCompute, particles.size());
    useComputeProgram(updateParticleCompute, particles.size());
    useComputeProgram(densitiesCompute, particles.size());
    useComputeProgram(pressureCompute, particles.size());
    useComputeProgram(viscosityCompute, particles.size());

    glUseProgram(renderingShaderProgram);
    
    // Draw particles using instanced data from quads
    glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr, numRows * numColumns);
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

    // Call method that creates and compiles shaders with shader sources that are #included from external files
    renderingShaderProgram = createShaderProgram(vertexShaderSource, fragmentShaderSource);
    
    //compile compute shaders
    updateParticleCompute = compileComputeShader(updateParticleComputeShaderSource);
    gravityCompute = compileComputeShader(gravityComputeShaderSource);
    densitiesCompute = compileComputeShader(calculateDensitiesComputeShaderSource);
    pressureCompute = compileComputeShader(pressureForceComputeShaderSource);
    viscosityCompute = compileComputeShader(viscosityForceComputeShaderSource);


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

    initializeParticles(particles);

    // Pre-size grid array to the max size
    grid.resize(gridWidth * gridHeight);

    int cellOffsets[9] = {-1 * gridWidth - 1, -gridWidth, -gridWidth + 1, -1, 0, 1, gridWidth - 1, gridWidth, gridWidth + 1};

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

    glGenBuffers(1, &mouseinfoSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, mouseinfoSSBO);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, mouseinfoSSBO);

    // Assign quad position instructions for vertex shader at layout position 0
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Initialize variables used for deltaTime and average FPS
    double lastTime = glfwGetTime();
    timeCount = 0;
    float deltaTimeSum = 0.0f;

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

        deltaTimeSum += deltaTime;
        timeCount++;

        // Keypress handling
        //if (commaPressed) frameCount++;
        if (spacePressed)
        {
            if (periodHeld)
            {
                runFrame();
            }
            else if (periodPressed)
            {
                runFrame();
                periodPressed = false;
            }
            
            if (commaPressed)
            {
                if (frameCount < 10)
                { 
                    runFrame();
                    frameCount++;
                }
                else
                {
                    frameCount = 0;
                    commaPressed = false;
                }
            }

            if (rPressed)
            {
                particles.clear();
                initializeParticles(particles);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, particleSSBO);
                glBufferData(GL_SHADER_STORAGE_BUFFER, particles.size() * sizeof(Particle), particles.data(), GL_STATIC_DRAW);
                runFrame();
                rPressed = false;
            }
        }
        else
        {
            runFrame();
        }

        glfwSwapBuffers(window);
        glfwPollEvents();

    }

    std::cout << "Average FPS: " << 1 / (deltaTimeSum / timeCount) << " " << timeCount << std::endl;

    // Cleanup
    glfwTerminate();
    return 0;
}