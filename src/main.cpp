#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/random.hpp>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw_gl3.h>
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
float pressureMultiplier = 0.5f;
float nearPressureMultiplier = 10.0f;
float gravity = 0.0f;
float Particle::restitution = 0.9f;

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
    int numParticles = 25;

    for (int i = 0; i < numParticles; i++)
    {
        for (int j = 0; j < numParticles; j++)
        {
            Particle temp(glm::vec2(50 + (i * 3 * Particle::radius), j * 3 * Particle::radius + 50), glm::vec2(0, 0));
            particles.push_back(temp);
        }
    }
}

float smoothingKernel(float smoothingRadius, float distance)
{
    if (distance >= smoothingRadius) return 0.0f;

    float scaleFactor = 6.0f / (M_PI * pow(smoothingRadius, 2));
    return (smoothingRadius - distance) * (smoothingRadius - distance) * scaleFactor;
}

float smoothingKernelDerivative(float smoothingRadius, float distance)
{
    if (distance >= smoothingRadius) return 0.0f;

    // float scaleFactor = 12.0f / (pow(smoothingRadius, 2) * M_PI);
    float scaleFactor = (pow(smoothingRadius, 2) * M_PI) / 12.0f;
    return (distance - smoothingRadius) * scaleFactor;
}

float nearSmoothingKernel(float smoothingRadius, float distance)
{
    if (distance >= smoothingRadius) return 0.0f;

    float scaleFactor = 10 / (M_PI * pow(smoothingRadius, 5));
    return (smoothingRadius - distance) * (smoothingRadius - distance) * (smoothingRadius - distance) * scaleFactor;   
}

float nearSmoothingKernelDerivative(float smoothingradius, float distance)
{
    if (distance >= smoothingRadius) return 0.0f;

    float scaleFactor = 30 / (pow(smoothingRadius, 5) * M_PI);
    return -1.0f * (smoothingRadius - distance) * (smoothingRadius - distance) * scaleFactor;
}

std::pair<float, float> calculateDensities(Particle p, std::vector<Particle>& particles, float smoothingRadius)
{
    float density = 0.0f;
    float nearDensity = 0.0f;
    const float mass = 100.0f;

    for (Particle other : particles)
    {
        float distance = glm::length(glm::vec2(other.position.x - p.position.x, other.position.y - p.position.y));
        float influence = smoothingKernel(smoothingRadius, distance);
        float nearInfluence = nearSmoothingKernel(smoothingRadius, distance);
        density += mass * influence;
        nearDensity += mass * nearInfluence;
    }

    return std::pair<float, float>(density, nearDensity);
}

std::pair<float, float> convertDensityToPressure(float density, float nearDensity, float targetDensity, float pressureMultiplier, float nearPressureMultiplier)
{
    float densityError = -density + targetDensity;
    float pressure = densityError * pressureMultiplier;
    float nearPressure = nearDensity * nearPressureMultiplier;
    return std::pair<float, float>(pressure, nearPressure);
}

std::pair<glm::vec2, glm::vec2> calculatePressureForce(int particleIndex, std::vector<Particle>& particles, std::vector<float>& densities, std::vector<float>& nearDensities, float smoothingRadius, float targetDensity, float pressureMultiplier)
{
    glm::vec2 pressureForce = glm::vec2(0.0f, 0.0f);
    glm::vec2 nearPressureForce = glm::vec2(0.0f, 0.0f);

    for (int i = 0; i < particles.size(); i++)
    {
        if (particleIndex == i) continue;

        float distance = glm::length(glm::vec2(particles[i].position.x - particles[particleIndex].position.x, particles[i].position.y - particles[particleIndex].position.y));
        glm::vec2 direction = glm::vec2(particles[i].position.x - particles[particleIndex].position.x, particles[i].position.y - particles[particleIndex].position.y) / distance;

        if (distance == 0) direction = glm::circularRand(1.0f);
        float slope = smoothingKernelDerivative(smoothingRadius, distance);
        float density = densities[i];
        float nearDensity = nearDensities[i];
        std::pair<float, float> pressurePair = convertDensityToPressure(density, nearDensity, targetDensity, pressureMultiplier, nearPressureMultiplier);
        pressureForce += -1.0f * pressurePair.first * 100.0f * direction * slope / glm::max(density, 1e-4f);

        // nearPressureForce += -1.0f * pressurePair.second * 100.0f * direction * slope / glm::max(density, 1e-4f);
    }

    return std::pair<glm::vec2, glm::vec2>(pressureForce, nearPressureForce);
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

    std::vector<Particle> particles;
    initializeParticles(particles);

    std::vector<float> densities;
    std::vector<float> nearDensities;

    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glGenBuffers(1, &meshVBO);
    glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
    glBufferData(GL_ARRAY_BUFFER, meshVertices.size() * sizeof(float), meshVertices.data(), GL_STATIC_DRAW);

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

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    double lastTime = glfwGetTime();
    int timeCount = 0;
    float deltaTimeSum = 0.0f;

    while (!glfwWindowShouldClose(window))
    {
        double nowTime = glfwGetTime();
        float deltaTime = (float) nowTime - lastTime;
        lastTime = nowTime;

        deltaTimeSum += deltaTime;
        timeCount++;

        // Clear screen each frame
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplGlfwGL3_NewFrame();

        densities.clear();
        nearDensities.clear();

        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].accelerate(glm::vec2(0.0f, -gravity * deltaTime));
            std::pair<float, float> tempDensity = calculateDensities(particles[i], particles, smoothingRadius);
            densities.push_back(tempDensity.first);
            nearDensities.push_back(tempDensity.second);
        }

        for (int i = 0; i < particles.size(); i++)
        {
            std::pair<glm::vec2, glm::vec2> calculatedForce = calculatePressureForce(i, particles, densities, nearDensities, smoothingRadius, targetDensity, pressureMultiplier);
            glm::vec2 pressureForce = calculatedForce.first;
            glm::vec2 nearPressureForce = calculatedForce.second;
            glm::vec2 pressureAcceleration = pressureForce / glm::max(densities[i], 1e-4f);
            glm::vec2 nearPressureAcceleration = nearPressureForce / glm::max(densities[i], 1e-4f);
            particles[i].accelerate(pressureAcceleration * deltaTime);
            //particles[i].accelerate(nearPressureAcceleration * -deltaTime);
            std::cout << pressureAcceleration.x * deltaTime << " " << pressureAcceleration.y * deltaTime << " " << deltaTime << std::endl;

        }
        
        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].updatePosition(deltaTime);
        }
        std::cout << std::endl;

        glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
        glBufferData(GL_ARRAY_BUFFER, particles.size()*sizeof(Particle), particles.data(), GL_DYNAMIC_DRAW);
        

        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, meshVertices.size(), particles.size());

        ImGui::Begin("Controls");
        ImGui::SliderFloat("PressureMultiplier", &pressureMultiplier, 1.0f, 1000.0f);
        ImGui::SliderFloat("Gravity", &gravity, 0.0f, 1000.0f);
        ImGui::SliderFloat("targetDensity", &targetDensity, 0.0f, 20.0f);
        ImGui::SliderFloat("restitution", &Particle::restitution, 0.0f, 1.0f);
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