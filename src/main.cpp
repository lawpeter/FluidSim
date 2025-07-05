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
#include "particle.cpp"

#ifndef DIMENSIONS

#define DIMENSIONS
#define SCREEN_WIDTH 800.0f
#define SCREEN_HEIGHT 800.0f

#endif

#ifndef INCLUDE

#define INCLUDE
#include <iostream>
#include <glm/glm.hpp>

#endif

const float smoothingRadius = 4.0f * Particle::radius;
float targetDensity = 2.75f;
float pressureMultiplier = 200.0f;
float nearPressureMultiplier = 1.0f;
float gravity = 1000.0f;
float Particle::restitution = 0.9f;
float mass = 1.0f;
float viscosityStrength = 0.3f;

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

void initializeMesh(std::vector<float>& meshVertices)
{
    int resolution = 40;
    
    meshVertices.push_back(0.0f);
    meshVertices.push_back(0.0f);

    for (int i = 0; i <= resolution; i++) {
        float angle = (2.0f * M_PI * i) / resolution;
        float x = (cos(angle) * Particle::radius);
        float y = (sin(angle) * Particle::radius);
        meshVertices.push_back(x);
        meshVertices.push_back(y);
    }
}

void initializeParticles(std::vector<Particle>& particles)
{
    int numParticles = 40;

    for (int i = 0; i < numParticles; i++)
    {
        for (int j = 0; j < numParticles; j++)
        {
            Particle temp(glm::vec2(50 + (i * 3 * Particle::radius), j * 3 * Particle::radius + 50), glm::vec2(0, 0));
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

glm::vec2 positionToCoordinate(glm::vec2 position)
{
    int cellX = (int) (position.x / smoothingRadius);
    int cellY = (int) (position.y / smoothingRadius);

    return glm::vec2(cellX, cellY);
}

int keyToHash(glm::vec2 cellKey)
{
    return ((cellKey.x * 15823) + (cellKey.y * 9737333));
}

void checkNearbyParticles(std::vector<int>& results, std::vector<glm::vec2>& cellOffsets, Particle& p, std::vector<Particle>& particles, std::vector<std::pair<int, int>>& cellOfParticles, std::vector<int>& startIndexOfCells)
{
    results.clear();

    glm::vec2 cellToCheck = positionToCoordinate(p.position);

    for (glm::vec2 offset : cellOffsets)
    {
        int keyToCheck = keyToHash(cellToCheck + offset) % particles.size();
        int index = startIndexOfCells[keyToCheck];
        for (int i = 0; i < 5; i++)
        {
            std::cout << cellOfParticles[index + i].first << " ";
        }
        std::cout << std::endl;

        while (cellOfParticles[index].first == keyToCheck)
        {
            results.push_back(cellOfParticles[index].second);
            index++;
        }
    }
}

void calculateDensities(std::vector<Particle>& particles, std::vector<Particle>& predictedParticles, std::vector<float>& densities, std::vector<float>& nearDensities, std::vector<int>& results, std::vector<glm::vec2>& cellOffsets, std::vector<std::pair<int, int>>& cellOfParticles, std::vector<int>& startIndexOfCells)
{
    for (int currentParticleIndex = 0; currentParticleIndex < particles.size(); currentParticleIndex++)
    {  
        glm::vec2 currentParticlePosition = predictedParticles[currentParticleIndex].position;

        float density = 0.0f;
        float nearDensity = 0.0f;

        checkNearbyParticles(results, cellOffsets, predictedParticles[currentParticleIndex], predictedParticles, cellOfParticles, startIndexOfCells);

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

void updateParticleCells(std::vector<Particle>& particles, std::vector<std::pair<int, int>>& cellOfParticles, std::vector<int>& startIndexOfCells)
{
    for (int i = 0; i < particles.size(); i++)
    {
        glm::vec2 cellKey = positionToCoordinate(particles[i].position);
        int cellHash = (keyToHash(cellKey) % particles.size());

        cellOfParticles[i] = std::pair<int, int>(cellHash, i);
        startIndexOfCells[i] = UINT_MAX;
    }

    std::sort(cellOfParticles.begin(), cellOfParticles.end());

    for (int i = 0; i < particles.size(); i++)
    {
        int currentKey = cellOfParticles[i].first;
        int previousKey = (i == 0 ? UINT_MAX : cellOfParticles[i - 1].second);
        if (currentKey != previousKey)
        {
            startIndexOfCells[currentKey] = i;
        }
    }
}

void constructCellOffsets(std::vector<glm::vec2>& cellOffsets)
{
    cellOffsets.push_back(glm::vec2(-1, -1));
    cellOffsets.push_back(glm::vec2(0, -1));
    cellOffsets.push_back(glm::vec2(1, -1));
    cellOffsets.push_back(glm::vec2(-1, 0));
    cellOffsets.push_back(glm::vec2(0, 0));
    cellOffsets.push_back(glm::vec2(1, 0));
    cellOffsets.push_back(glm::vec2(-1, 1));
    cellOffsets.push_back(glm::vec2(0, 1));
    cellOffsets.push_back(glm::vec2(1, 1));
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

    GLuint VAO;
    GLuint meshVBO;
    GLuint particleVBO;

    GLuint shaderProgram = createShaderProgram(vertexShaderSource, fragmentShaderSource);

    GLint projectionUniformLocation = glGetUniformLocation(shaderProgram, "projection");
    glUniformMatrix4fv(projectionUniformLocation, 1, GL_FALSE, &projection[0][0]);

    std::vector<float> meshVertices;
    initializeMesh(meshVertices);

    /*float quadPositions[8] = {
        -0.5f, -0.5f, // 0
         0.5f, -0.5f, // 1
         0.5f,  0.5f, // 2
        -0.5f,  0.5f  // 3
    };

    GLuint indices[6] = {
        0, 1, 2,
        2, 3, 0
    };*/

    std::vector<Particle> particles;
    initializeParticles(particles);

    std::vector<std::pair<int, int>> cellOfParticles;
    std::vector<int> startIndexOfCells;
    cellOfParticles.resize(particles.size());
    startIndexOfCells.resize(particles.size());

    std::vector<glm::vec2> cellOffsets;
    constructCellOffsets(cellOffsets);

    std::vector<int> results;

    std::vector<float> densities;
    std::vector<float> nearDensities;

    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glGenBuffers(1, &meshVBO);
    glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
    glBufferData(GL_ARRAY_BUFFER, meshVertices.size() * sizeof(float), meshVertices.data(), GL_STATIC_DRAW);

    /*glGenBuffers(1, &meshVBO);
    glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
    glBufferData(GL_ARRAY_BUFFER, quadPositions.size() * sizeof(float), meshVertices.data(), GL_STATIC_DRAW);*/

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

    // particleVelocity > location 1, vec2
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Particle, velocity));
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    double lastTime = glfwGetTime();
    int timeCount = 0;
    float deltaTimeSum = 0.0f;

    std::vector<Particle> predictedParticles;
    for (Particle p : particles)
    {
        predictedParticles.push_back(p);
    }

    while (!glfwWindowShouldClose(window)/* && particles[20].position.x > 0*/)
    {
        /*double nowTime = glfwGetTime();
        float deltaTime = (float) nowTime - lastTime;
        lastTime = nowTime;*/
        float deltaTime = 1 / 120.0f;

        deltaTimeSum += deltaTime;
        timeCount++;

        // Clear screen each frame
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplGlfwGL3_NewFrame();

        densities.clear();
        nearDensities.clear();

        updateParticleCells(particles, cellOfParticles, startIndexOfCells);

        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].accelerate(glm::vec2(0.0f, -gravity * deltaTime));
            predictedParticles[i].position = particles[i].position + particles[i].velocity * deltaTime;
        }

        calculateDensities(particles, predictedParticles, densities, nearDensities, results, cellOffsets, cellOfParticles, startIndexOfCells);

        for (int i = 0; i < particles.size(); i++)
        {
            /*checkNearbyParticles(results, cellOffsets, particles[i], particles, cellOfParticles, startIndexOfCells);

            glm::vec2 pressureForce = calculatePressureForce(i, particles, predictedParticles, densities, nearDensities, results);
            glm::vec2 pressureAcceleration = pressureForce / glm::max(densities[i], 1e-4f);
            particles[i].accelerate(pressureAcceleration * deltaTime);

            glm::vec2 viscosityForce = calculateViscosity(i, particles, predictedParticles, results);
            glm::vec2 viscosityAcceleration = viscosityForce / glm::max(densities[i], 1e-4f);
            particles[i].accelerate(viscosityAcceleration * deltaTime);*/

            /*std::cout << "Particle: " << i << std::endl;
            std::cout << "Position: " << particles[i].position.x << " " << particles[i].position.y << std::endl;
            std::cout << "Velocity: " << particles[i].velocity.x << " " << particles[i].velocity.y << std::endl;
            std::cout << "Acceleration: " << pressureAcceleration.x << " " << pressureAcceleration.y << std::endl;
            std::cout << "Force: " << pressureForce.x << " " << pressureForce.y << std::endl;
            std::cout << "Density: " << densities[i] << std::endl;
            std::cout << std::endl;*/
            
        }

        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].updatePosition(deltaTime);
        }        

        glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
        glBufferData(GL_ARRAY_BUFFER, particles.size()*sizeof(Particle), particles.data(), GL_DYNAMIC_DRAW);
        

        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, meshVertices.size(), particles.size());


        ImGui::Begin("Controls");
        ImGui::SliderFloat("PressureMultiplier", &pressureMultiplier, 0.0f, 500.0f);
        ImGui::SliderFloat("nearPressureMultiplier", &nearPressureMultiplier, 0.0f, 1000.0f);
        ImGui::SliderFloat("Gravity", &gravity, 0.0f, 1000.0f);
        ImGui::SliderFloat("targetDensity", &targetDensity, 0.0f, 10.0f);
        ImGui::SliderFloat("restitution", &Particle::restitution, 0.0f, 1.0f);
        ImGui::SliderFloat("viscosityStrength", &viscosityStrength, 0.0f, 10.0f);
        ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
        ImGui::End();

        ImGui::Render();
        ImGui_ImplGlfwGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
        glfwPollEvents();
        
    }

    std::cout << "Average FPS: " << 1 / (deltaTimeSum / timeCount) << std::endl;

    ImGui_ImplGlfwGL3_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}