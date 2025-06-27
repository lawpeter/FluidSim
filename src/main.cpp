#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
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

#endif

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

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
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

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
        float x = (cos(angle) * 10);
        float y = (sin(angle) * 10);
        meshVertices.push_back(x);
        meshVertices.push_back(y);
    }
}

void initializeParticles(std::vector<Particle>& particles)
{
    int numParticles = 10;
    int scale = 10;

    for (int i = 0; i < numParticles; i++)
    {
        Particle temp(glm::vec2(100 + (i * 3 * scale), 100), glm::vec2(0, 0));
        particles.push_back(temp);
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

        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].accelerate(glm::vec2(0.0f, -9.8f * deltaTime));
        }
        
        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].updatePosition(deltaTime);
        }

        glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
        glBufferData(GL_ARRAY_BUFFER, particles.size()*sizeof(Particle), particles.data(), GL_DYNAMIC_DRAW);
        

        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, meshVertices.size(), particles.size());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    std::cout << "Average FPS: " << 1 / (deltaTimeSum / timeCount) << std::endl;

    glfwTerminate();
    return 0;
}