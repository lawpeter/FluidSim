const char* gravityComputeShaderSource = R"glsl(
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

    vec2 gravity = vec2(0.0, -9.81);

    void main()
    {
        uint gID = gl_GlobalInvocationID.x;

        if (gID >= particles.length()) return; 

        //particles[gID].velocity += gravity;
    }
)glsl";