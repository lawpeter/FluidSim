#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/random.hpp>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw_gl3.h>
#include <limits.h>
#include "fragment.glsl"
#include "vertex.glsl"
#include <iostream>
#include "particle.hpp"
#include <vector>

#define gridWidth ceil(SCREEN_WIDTH/smoothingRadius)
#define gridHeight ceil(SCREEN_HEIGHT/smoothingRadius)


const float smoothingRadius = 4.0f * Particle::radius;
float targetDensity = 2.75f;
float pressureMultiplier = 10.0f;
float nearPressureMultiplier = 1.0f;
float gravity = 98.1f;
float Particle::restitution = 0.9f;
float mass = 1.0f;
float viscosityStrength = 100.0f;
float mouseInteractionRadius = 200.0f;
float mouseForceStrength = 500.0f;

bool spacePressed = true;
bool periodPressed = false;
bool periodHeld = false;
bool commaPressed = false;
bool rPressed = false;
int frameCount = 0;

float cursorX;
float cursorY;
bool leftMousePressed = false;
bool rightMousePressed = false;

int numRows = 50;
int numColumns = 50;

GLFWwindow* loadSim()
{
    if (!glfwInit())
    {
        std::cout << "brokey" << std::endl;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window;
    window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "FluidSim", NULL, NULL);

    if (window == NULL)
    {
        std::cout << "brokey 2" << std::endl;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "brokey 3" << std::endl;
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

    // Activate the shader program
    glUseProgram(shaderProgram);

    return shaderProgram;
}

void initializeParticles(std::vector<Particle>& particles)
{
    for (int i = 0; i < numColumns; i++)
    {
        for (int j = 0; j < numRows; j++)
        {
            Particle temp(glm::vec2(50 + (i * 2.5f * Particle::radius), (50 + j * 2.5f * Particle::radius)), glm::vec2(0, 0));
            particles.push_back(temp);
        }
    }
}

float densityKernel(float distance)
{
    if (distance >= smoothingRadius) return 0.0f;

    float scaleFactor = 6.0f / (M_PI * pow(smoothingRadius, 2));
    return (smoothingRadius - distance) * (smoothingRadius - distance) * scaleFactor;
}

float densityKernelDerivative(float distance)
{
    if (distance >= smoothingRadius) return 0.0f;

    float scaleFactor = (pow(smoothingRadius, 2) * M_PI) / 12.0f;
    return -(smoothingRadius - distance) * scaleFactor;
}

float nearDensityKernel(float distance)
{
    if (distance >= smoothingRadius) return 0.0f;

    float scaleFactor = 10.0f / (M_PI * pow(smoothingRadius, 5));
    return (smoothingRadius - distance) * (smoothingRadius - distance) * (smoothingRadius - distance) * scaleFactor;
}

float nearDensityKernelDerivative(float distance)
{
    if (distance >= smoothingRadius) return 0.0f;

    float scaleFactor = 30.0f / (M_PI * pow(smoothingRadius, 5));
    return -(smoothingRadius - distance) * (smoothingRadius - distance) * scaleFactor;
}

int positionToCellArrayIndex(glm::vec2 position)
{
    int cellX = (int) (position.x / smoothingRadius);
    int cellY = (int) (position.y / smoothingRadius);

    return (cellY * gridWidth) + cellX;
}

void checkNearbyParticles(std::vector<int>& results, std::vector<int>& cellOffsets, Particle& particle, std::vector<Particle>& particles, std::vector<std::vector<int>>& grid)
{
    results.clear();

    int cell = positionToCellArrayIndex(particle.position);

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

void calculateDensities(std::vector<Particle>& particles, std::vector<Particle>& predictedParticles, std::vector<float>& densities, std::vector<float>& nearDensities, std::vector<int>& results, std::vector<int>& cellOffsets, std::vector<std::vector<int>>& grid)
{
    for (int currentParticleIndex = 0; currentParticleIndex < particles.size(); currentParticleIndex++)
    {  
        glm::vec2 currentParticlePosition = predictedParticles[currentParticleIndex].position;

        float density = 0.0f;
        float nearDensity = 0.0f;

        checkNearbyParticles(results, cellOffsets, predictedParticles[currentParticleIndex], predictedParticles, grid);

        for (int otherParticleIndex : results)
        {
            glm::vec2 otherParticlePosition = predictedParticles[otherParticleIndex].position;
            glm::vec2 offsetToOtherParticle = otherParticlePosition - currentParticlePosition;
            float sqrDistanceToOtherParticle = glm::dot(offsetToOtherParticle, offsetToOtherParticle);

            if (sqrDistanceToOtherParticle > pow(smoothingRadius, 2)) continue;

            float distance = sqrt(sqrDistanceToOtherParticle);
            float influence = densityKernel(distance);
            float nearInfluence = nearDensityKernel(distance);

            density += mass * influence;
            nearDensity += mass * nearInfluence;
        }

        densities.push_back(density);
        nearDensities.push_back(nearDensity);
    }
}

float convertDensityToPressure(float density)
{
    float densityError = density - targetDensity;
    float pressure = densityError * pressureMultiplier;
    return pressure;
}

float convertNearDensityToPressure(float nearDensity)
{
	return nearPressureMultiplier * nearDensity;
}

glm::vec2 calculatePressureForce(int particleIndex, std::vector<Particle>& particles, std::vector<Particle>& predictedParticles, std::vector<float>& densities, std::vector<float>& nearDensities, std::vector<int>& results)
{   
    glm::vec2 currentParticlePosition = predictedParticles[particleIndex].position;
    float currentDensity = densities[particleIndex];
    float currentNearDensity = nearDensities[particleIndex];
    float currentPressure = convertDensityToPressure(currentDensity);
    float currentNearPressure = convertNearDensityToPressure(currentNearDensity);

    glm::vec2 pressureForce(0, 0);

    int maxIndex = 0;

    for (int otherParticleIndex : results) 
    {
        if (otherParticleIndex == particleIndex) continue;
        glm::vec2 otherParticlePosition = predictedParticles[otherParticleIndex].position;
        
        glm::vec2 offsetToOtherParticle = otherParticlePosition - currentParticlePosition;
        float sqrDistanceToOtherParticle = glm::dot(offsetToOtherParticle, offsetToOtherParticle);
        
        if (sqrDistanceToOtherParticle > pow(smoothingRadius, 2)) continue;

        float distance = sqrt(sqrDistanceToOtherParticle);
        glm::vec2 dirToOtherParticle = distance > 0 ? offsetToOtherParticle / distance : glm::vec2
        (0, 1);

        float otherDensity = densities[otherParticleIndex];
        float otherNearDensity = nearDensities[otherParticleIndex];
        float otherPressure = convertDensityToPressure(otherDensity);
        float otherNearPressure = convertNearDensityToPressure(otherNearDensity);
    
        float sharedPressure = (currentPressure + otherPressure) * 0.5f;
        float sharedNearPressure = (currentNearPressure + otherNearPressure) * 0.5f;

        pressureForce += dirToOtherParticle * densityKernelDerivative(distance) * sharedPressure / otherDensity;
        pressureForce += dirToOtherParticle * nearDensityKernelDerivative(distance) * sharedNearPressure / otherNearDensity;
    }

    return pressureForce;
}

float viscosityKernel(float squareDistance)
{
    if (squareDistance >= smoothingRadius * smoothingRadius) return 0.0f;
    
    float scale = 4 / (M_PI * pow(smoothingRadius, 8));
    return pow(((smoothingRadius * smoothingRadius) - squareDistance), 3) * scale;
}

glm::vec2 calculateViscosity(int particleIndex, std::vector<Particle>& particles, std::vector<Particle>& predictedParticles, std::vector<int>& results)
{
    glm::vec2 viscosityForce(0, 0);
    glm::vec2 currentParticlePosition = predictedParticles[particleIndex].position;
    for (int otherParticleIndex : results)
    {
        if (otherParticleIndex == particleIndex) continue;

        glm::vec2 otherParticlePosition = predictedParticles[otherParticleIndex].position;
        
        glm::vec2 offsetToOtherParticle = otherParticlePosition - currentParticlePosition;
        float sqrDistanceToOtherParticle = glm::dot(offsetToOtherParticle, offsetToOtherParticle);

        if (sqrDistanceToOtherParticle > pow(smoothingRadius, 2)) continue;
        
        glm::vec2 otherVelocity = particles[otherParticleIndex].velocity;
        viscosityForce += (otherVelocity - particles[particleIndex].velocity) * viscosityKernel(sqrDistanceToOtherParticle);
    }

    return viscosityForce * viscosityStrength;
}

void updateParticleCells(std::vector<Particle>& particles, std::vector<std::vector<int>>& grid)
{
    grid.clear();
    grid.resize(gridWidth * gridHeight);

    for (int i = 0; i < particles.size(); i++)
    {
        int cellIndex = positionToCellArrayIndex(particles[i].position);
        if (cellIndex < grid.size() && cellIndex >= 0) 
        {
            grid[cellIndex].push_back(i);
        }
        else
        {
            std::cout << "ROGUE PARTICLE" << std::endl;
        }
    }
}

void constructCellOffsets(std::vector<int>& cellOffsets)
{
    cellOffsets.push_back(-1);
    cellOffsets.push_back(0);
    cellOffsets.push_back(1);
    cellOffsets.push_back(-gridWidth - 1);
    cellOffsets.push_back(-gridWidth);
    cellOffsets.push_back(-gridWidth + 1);
    cellOffsets.push_back(gridWidth - 1);
    cellOffsets.push_back(gridWidth);
    cellOffsets.push_back(gridWidth + 1);
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

int main(int argc, char* argv[])
{
    glm::mat4 projection = glm::ortho(0.0f, SCREEN_WIDTH, 0.0f, SCREEN_HEIGHT, -1.0f, 1.0f);

    GLFWwindow* window = loadSim();

    if (window == NULL)
    {
        return -1;
    }

    ImGui::CreateContext();
    ImGui_ImplGlfwGL3_Init(window, true);
    ImGui::StyleColorsDark();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    GLuint VAO;
    GLuint quadVBO;
    GLuint indexBuffer;
    GLuint particleVBO;
    GLuint instanceVBO;

    GLuint shaderProgram = createShaderProgram(vertexShaderSource, fragmentShaderSource);

    GLint projectionUniformLocation = glGetUniformLocation(shaderProgram, "projection");
    glUniformMatrix4fv(projectionUniformLocation, 1, GL_FALSE, &projection[0][0]);

    GLint radiusUniformLocation = glGetUniformLocation(shaderProgram, "radius");
    glUniform1f(radiusUniformLocation, Particle::radius);

    float quadVertices[8] = {
        -Particle::radius, -Particle::radius, // 0
         Particle::radius, -Particle::radius, // 1
         Particle::radius,  Particle::radius, // 2
        -Particle::radius,  Particle::radius  // 3
    };

    GLuint indices[6] = {
        0, 1, 2,
        2, 3, 0
    };

    std::vector<Particle> particles;
    initializeParticles(particles);

    std::vector<std::vector<int>> grid;
    grid.resize(gridWidth * gridHeight);

    std::vector<int> cellOffsets;
    constructCellOffsets(cellOffsets);

    std::vector<int> results;

    std::vector<float> densities;
    std::vector<float> nearDensities;

    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glGenBuffers(1, &quadVBO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, 8 * sizeof(float), quadVertices, GL_STATIC_DRAW);

    glGenBuffers(1, &indexBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * sizeof(GLuint), indices, GL_STATIC_DRAW);

    // Mesh position
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &particleVBO);
    glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
    glBufferData(GL_ARRAY_BUFFER, particles.size() * sizeof(Particle), particles.data(), GL_DYNAMIC_DRAW);

    GLsizei stride = sizeof(Particle);

    // particlePosition > location 1, vec2
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Particle, position));
    glEnableVertexAttribArray(1);
    glVertexAttribDivisor(1, 1);

    // particleVelocity > location 2, vec2
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Particle, velocity));
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    double lastTime = glfwGetTime();
    int timeCount = 0;
    float deltaTimeSum = 0.0f;

    std::vector<Particle> predictedParticles;
    for (Particle p : particles)
    {
        predictedParticles.push_back(p);
    }

    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    while (!glfwWindowShouldClose(window))
    {
        double nowTime = glfwGetTime();
        float deltaTime = (float) nowTime - lastTime;
        lastTime = nowTime;
        //float deltaTime = 1 / 240.0f;

        deltaTimeSum += deltaTime;
        timeCount++;

        // Keypress handling
        if (commaPressed)
            frameCount++;
        
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

        ImGui_ImplGlfwGL3_NewFrame();

        densities.clear();
        nearDensities.clear();

        
        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].accelerate(glm::vec2(0.0f, -gravity * deltaTime));
            predictedParticles[i].position = particles[i].position + particles[i].velocity * deltaTime;
        }
        
        updateParticleCells(predictedParticles, grid);

        calculateDensities(particles, predictedParticles, densities, nearDensities, results, cellOffsets, grid);

        for (int i = 0; i < particles.size(); i++)
        {
            // Populate results with particles to check
            checkNearbyParticles(results, cellOffsets, particles[i], particles, grid);

            // Apply pressure forces
            glm::vec2 pressureForce = calculatePressureForce(i, particles, predictedParticles, densities, nearDensities, results);
            glm::vec2 pressureAcceleration = pressureForce / glm::max(densities[i], 1e-4f);
            particles[i].accelerate(pressureAcceleration * deltaTime);

            // Apply viscosity forces
            glm::vec2 viscosityForce = calculateViscosity(i, particles, predictedParticles, results);
            glm::vec2 viscosityAcceleration = viscosityForce / glm::max(densities[i], 1e-4f);
            particles[i].accelerate(viscosityAcceleration * deltaTime);

            // Mouse interaction
            if (leftMousePressed || rightMousePressed)
            {
                glm::vec2 cursorPosition(cursorX, SCREEN_HEIGHT - cursorY);
                glm::vec2 directionToCursor = cursorPosition - particles[i].position;
                float distance = glm::length(directionToCursor);
                
                if (distance > 1e-5f && distance < mouseInteractionRadius)
                {
                    glm::vec2 normalizedDireciton = directionToCursor / distance;
                    
                    float distanceFactor = 1.0f - (distance / mouseInteractionRadius);
                    float strength = mouseForceStrength * distanceFactor; // Strength falls off as you get further from the cursor
                    
                    if (leftMousePressed)
                    {
                        particles[i].accelerate(-normalizedDireciton * strength * deltaTime);
                    }
                    if (rightMousePressed)
                    {
                        particles[i].accelerate(normalizedDireciton * strength * deltaTime);
                    }
                }
            }
        }

        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].updatePosition(deltaTime);
        }

        glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
        glBufferData(GL_ARRAY_BUFFER, particles.size()*sizeof(Particle), particles.data(), GL_DYNAMIC_DRAW);
        
        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr, numRows * numColumns);

        ImGui::Begin("Controls");
        ImGui::SliderFloat("PressureMultiplier", &pressureMultiplier, 0.0f, 50.0f);
        ImGui::SliderFloat("nearPressureMultiplier", &nearPressureMultiplier, 0.0f, 100.0f);
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

    ImGui_ImplGlfwGL3_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}