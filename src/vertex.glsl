const char* vertexShaderSource = R"glsl(
    #version 330 core

    layout (location = 0) in vec2 meshPosition;
    layout (location = 1) in vec2 particlePosition;

    uniform mat4 projection;
 
    void main() {
        gl_Position = projection * vec4(meshPosition + particlePosition, 0.0, 1.0);
    }
)glsl";