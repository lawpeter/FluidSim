const char* computeShaderSource = R"glsl(
    #version 460 core

    layout(local_size_x = 128) in;

    struct Particle
    {
        vec2 position;
        vec2 velocity;
        float density;
    };

    layout (std430, binding = 0) buffer ParticleBuffer
    {
        Particle particles[];
    };


    const vec2 gravity = vec2(0.0, -98.1);
    const float deltaTime = 0.008;

    void main()
    {
        uint gID = gl_GlobalInvocationID.x;

        if (gID >= 2500) return; 

        particles[gID].velocity = vec2(10, 10);
        particles[gID].position = vec2(10, 10);
    }
)glsl";