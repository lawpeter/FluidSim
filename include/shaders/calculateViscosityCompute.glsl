const char* viscosityForceComputeShaderSource = R"glsl(
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

    layout (std430, binding = 1) buffer OneDBuffer
    {
        int oneDimensionalGrid[];
    };

    layout (std430, binding = 2) buffer StartIndicesBuffer
    {
        int startIndices[];
    };

    layout (std430, binding = 3) buffer EndIndicesBuffer
    {
        int endIndices[];
    };
    
    const float smoothingRadius = 8.0;
    const float M_PI = 3.14159265358979323846264338327950288419716939937510;
    float viscosityStrength = 100.0;
    float deltaTime = 0.0008;
    int SCREEN_WIDTH = 1000;
    int SCREEN_HEIGHT = 800;
    int gridWidth = int(ceil(SCREEN_WIDTH/smoothingRadius));
    int gridHeight = int(ceil(SCREEN_HEIGHT/smoothingRadius));
    int cellOffsets[9] =  int[9](-gridWidth - 1, -gridWidth, -gridWidth + 1, -1, 0, 1, gridWidth - 1, gridWidth, gridWidth + 1);

    int positionToCellArrayIndex(vec2 position)
    {
        // makes a grid based on smoothing radius
        int cellX = int(position.x / smoothingRadius);
        int cellY = int(position.y / smoothingRadius);

        return (cellY * gridWidth) + cellX;
    }

    float viscosityKernel(float squareDistance)
    {
        if (squareDistance >= smoothingRadius * smoothingRadius) return 0.0f;
        
        float scale = 4 / (M_PI * pow(smoothingRadius, 8));
        return pow(((smoothingRadius * smoothingRadius) - squareDistance), 3) * scale;
    }

    void main()
    {
        uint gID = gl_GlobalInvocationID.x;

        if (gID >= particles.length()) return; 

        vec2 viscosityForce = vec2(0.0, 0.0);
        vec2 currentParticlePosition = particles[gID].position + particles[gID].velocity * deltaTime;

        int gridIndex = positionToCellArrayIndex(currentParticlePosition);
        for (int i = 0; i < 9; i++)
        {
            int gridToCheck = gridIndex + cellOffsets[i];

            if (gridToCheck < 0 || gridToCheck >= gridWidth * gridHeight) continue;

            for (int oneDimensionalIndex = startIndices[gridToCheck]; oneDimensionalIndex < endIndices[gridToCheck]; oneDimensionalIndex++)
            {
                int otherParticleIndex = oneDimensionalGrid[oneDimensionalIndex];

                // Skip if checking itself
                if (otherParticleIndex == gID) continue;

                vec2 otherParticlePosition = particles[otherParticleIndex].position + particles[otherParticleIndex].velocity * deltaTime;
                vec2 offsetToOtherParticle = otherParticlePosition - currentParticlePosition;
                float sqrDistanceToOtherParticle = dot(offsetToOtherParticle, offsetToOtherParticle);

                // Skip if other particle is outside smoothing radius
                if (sqrDistanceToOtherParticle > pow(smoothingRadius, 2)) continue;
                
                vec2 otherVelocity = particles[otherParticleIndex].velocity;
                viscosityForce += (otherVelocity - particles[gID].velocity) * viscosityKernel(sqrDistanceToOtherParticle);  
            }
        }

        if (particles[gID].density != 0.0)
        particles[gID].velocity += viscosityForce * viscosityStrength * deltaTime / particles[gID].density;
                
    }
)glsl";