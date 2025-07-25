const char* updateParticleComputeShaderSource = R"glsl(
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
    const float deltaTime = 0.0008;
    const float radius = 2.0;
    const float restitution = 0.9;
    const float SCREEN_WIDTH = 1000.0;
    const float SCREEN_HEIGHT = 800.0;

    void main()
    {
        uint gID = gl_GlobalInvocationID.x;
        
        if (gID >= particles.length()) return; 

        particles[gID].position += particles[gID].velocity * deltaTime;

        if (particles[gID].position.x - radius < 0.0)
        {
            particles[gID].position.x = radius;
            particles[gID].velocity.x *= -restitution;
        }
        else if (particles[gID].position.x + radius > SCREEN_WIDTH)
        {
            particles[gID].position.x = SCREEN_WIDTH - radius;
            particles[gID].velocity.x *= -restitution;
        }
        if  (particles[gID].position.y - radius < 0.0)
        {
            particles[gID].position.y = radius;
            particles[gID].velocity.y *= -restitution;           
        }
        else if (particles[gID].position.y + radius > SCREEN_HEIGHT)
        {
            particles[gID].position.y = SCREEN_HEIGHT - radius;
            particles[gID].velocity.y *= -restitution;
        }

    }
)glsl";