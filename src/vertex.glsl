const char* vertexShaderSource = R"glsl(
    #version 330 core

    layout (location = 0) in vec2 meshPosition;
    layout (location = 1) in vec2 particlePosition;
    layout (location = 2) in vec2 particleVelocity;

    out float particleSpeed;

    uniform mat4 projection;
 
    void main() {
        particleSpeed = (min(sqrt(dot(particleVelocity, particleVelocity)), 800.0f) / 800.0f);
        gl_Position = projection * vec4(meshPosition + particlePosition, 0.0, 1.0);
    }
)glsl";