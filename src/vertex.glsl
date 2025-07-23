const char* vertexShaderSource = R"glsl(
    #version 460 core

    struct Particle
    {
        vec2 position;
        vec2 velocity;
        float density;
    };

    layout (location = 0) in vec2 meshPosition;
    uint instanceID = uint(gl_InstanceID);
    /*layout (location = 1) in vec2 particlePosition;
    layout (location = 2) in vec2 particleVelocity;*/

    layout(std430, binding = 0) buffer ParticleBuffer
    {
        Particle particles[];
    };

    out float particleSpeed;
    out vec2 localPosition;

    uniform mat4 projection;
    uniform float radius;
 
    void main() {
        particleSpeed = (min(sqrt(dot(particles[gl_InstanceID].velocity, particles[gl_InstanceID].velocity)), 300.0f) / 300.0f);
        localPosition = meshPosition / radius;
        gl_Position = projection * vec4(meshPosition + particles[gl_InstanceID].position, 0.0, 1.0);
    }

    
)glsl";