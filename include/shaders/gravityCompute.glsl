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

    layout (std430, binding = 4) buffer MouseBuffer
    {
        float mouseinfo[4];
    };


    vec2 gravity = vec2(0.0, 0.0);
    int SCREEN_WIDTH = 1800;
    int SCREEN_HEIGHT = 1000;
    float mouseInteractionRadius = 60;
    float mouseForceStrength = 400;

    void main()
    {
        uint gID = gl_GlobalInvocationID.x;

        if (gID >= particles.length()) return; 

        particles[gID].velocity += gravity;

        bool leftMousePressed = mouseinfo[0] == 1;
        bool rightMousePressed = mouseinfo[1] == 1;
        float cursorX = mouseinfo[2];
        float cursorY = mouseinfo[3];
        const float deltaTime = 0.0022222222222;

        if (leftMousePressed || rightMousePressed)
        {
            vec2 cursorPosition = vec2(cursorX, SCREEN_HEIGHT - cursorY);
            vec2 directionToCursor = cursorPosition - (particles[gID].position + particles[gID].velocity * deltaTime);
            float distance = length(directionToCursor);
            
            if (distance > 1e-5f && distance < mouseInteractionRadius)
            {
                vec2 normalizedDireciton = directionToCursor / distance;
                
                float distanceFactor = 1.0f - (distance / mouseInteractionRadius);
                float strength = mouseForceStrength * distanceFactor; // Strength falls off as you get further from the cursor
                
                if (leftMousePressed)
                {
                    // Pull toward from mouse position (due to positive)
                    particles[gID].velocity += vec2(normalizedDireciton * strength);
                }
                if (rightMousePressed)
                {
                    // Push away from mouse position (due to negative)
                    particles[gID].velocity += vec2(-normalizedDireciton * strength);
                }
            }
        }
    }
)glsl";