const char* pressureForceComputeShaderSource = R"glsl(
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

    float targetDensity = 0.005;
    const float smoothingRadius = 4.0;
    const float M_PI = 3.14159265358979323846264338327950288419716939937510;
    float deltaTime = 0.0022222222222;
    float pressureMultiplier = 0.16 / deltaTime;
    int SCREEN_WIDTH = 1800;
    int SCREEN_HEIGHT = 1000;
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

    float convertDensityToPressure(float density)
    {
        float densityError = density - targetDensity;
        float pressure = densityError * pressureMultiplier;
        return pressure;
    } 
        
    float densityKernelDerivative(float distance)
    {
        if (distance >= smoothingRadius) return 0.0f;

        float scaleFactor = (pow(smoothingRadius, 2) * M_PI) / 12.0f;
        return -(smoothingRadius - distance) * scaleFactor;
    }

    void main()
    {
        uint gID = gl_GlobalInvocationID.x;

        if (gID >= particles.length()) return; 

        vec2 currentParticlePosition = particles[gID].position + particles[gID].velocity * deltaTime;
        float currentDensity = particles[gID].density;
        float currentPressure = convertDensityToPressure(currentDensity);

        vec2 pressureForce = vec2(0.0, 0.0);

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
                
                // Skip if outside smoothing radius
                if (sqrDistanceToOtherParticle > pow(smoothingRadius, 2)) continue;

                float distance = sqrt(sqrDistanceToOtherParticle);
                vec2 dirToOtherParticle = distance > 0 ? offsetToOtherParticle / distance : vec2(0, 1);

                float otherDensity = particles[otherParticleIndex].density;
                float otherPressure = convertDensityToPressure(otherDensity);
            
                // Apply the average of the two pressures to each particle
                float sharedPressure = (currentPressure + otherPressure) * 0.5;

                if (otherDensity != 0.0)
                {
                    pressureForce += dirToOtherParticle * densityKernelDerivative(distance) * sharedPressure / otherDensity;
                }

            }
        }

        if (particles[gID].density != 0.0)
        particles[gID].velocity += pressureForce / particles[gID].density;

    }
)glsl";