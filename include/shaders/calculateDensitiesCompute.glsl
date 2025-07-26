const char* calculateDensitiesComputeShaderSource = R"glsl(
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
        int oneD[];
    };

    layout (std430, binding = 2) buffer StartIndicesBuffer
    {
        int startIndices[];
    };
    
    layout (std430, binding = 3) buffer EndIndicesBuffer
    {
        int endIndices[];
    };

    const float mass = 1.0;
    const float smoothingRadius = 80.0;
    const float M_PI = 3.14159265358979323846264338327950288419716939937510;
    const float deltaTime = 0.0008;
    int SCREEN_WIDTH = 1000;
    int SCREEN_HEIGHT = 800;
    int gridWidth = int(ceil(SCREEN_WIDTH/smoothingRadius));
    int gridHeight = int(ceil(SCREEN_HEIGHT/smoothingRadius));
    int cellOffsets[9] =  int[9](-gridWidth - 1, -gridWidth, -gridWidth + 1, -1, 0, 1, gridWidth - 1, gridWidth, gridWidth + 1);

    int positionToCellArrayIndex(vec2 position)
    {
        // makes a grid based on smoothing radius
        int cellX = int((position.x / smoothingRadius));
        int cellY = int((position.y / smoothingRadius));

        return (cellY * gridWidth) + cellX;
    }

    float densityKernel(float distance)
    {
        if (distance >= smoothingRadius) return 0.0f;

        float scaleFactor = 6.0f / (M_PI * pow(smoothingRadius, 2));
        return (smoothingRadius - distance) * (smoothingRadius - distance) * scaleFactor;
    }

    void main()
    {
        uint gID = gl_GlobalInvocationID.x;
        //gID = 15960;

        if (gID >= particles.length()) return; 

        vec2 currentParticlePosition = particles[gID].position + particles[gID].velocity * deltaTime;
        float density = 0.0;

        //int cell = positionToCellArrayIndex(currentParticlePosition);
        //cell = 1598;

        for (int i = 0; i < startIndices.length(); i++)
        {
            for (int j = startIndices[i]; j < endIndices[i]; j++)
            {
                particles[j].velocity = vec2(i, 0);
            }
        }


        /*int cellX = int((currentParticlePosition.x / smoothingRadius));
        int cellY = int((currentParticlePosition.y / smoothingRadius));
        
        for (int vertical = -1; vertical <= 1; vertical++)
        {
            for (int horizontal = -1; horizontal <= 1; horizontal++)
            {
                int neighborX = cellX + horizontal;
                int neighborY = cellY + vertical;
                if (neighborX < 0 || neighborY < 0 || neighborX >= gridWidth || neighborY >= gridHeight) continue;

                int cellToCheck = neighborY * gridWidth + neighborX;
                
                if (cellToCheck < gridWidth * gridHeight && cellToCheck >= 0)
                {
                    //int endIndex = cellToCheck + 1 < gridWidth * gridHeight ? startIndices[cellToCheck + 1] : particles.length();
                    for (int otherParticleIndex = startIndices[cellToCheck]; otherParticleIndex < endIndices[cellToCheck]; otherParticleIndex++)
                    {
                        if (oneD[otherParticleIndex] == -1) break;
                        vec2 otherParticlePosition = particles[oneD[otherParticleIndex]].position + particles[oneD[otherParticleIndex]].velocity * deltaTime;
                        vec2 offsetToOtherParticle = otherParticlePosition - currentParticlePosition;
                        float sqrDistanceToOtherParticle = dot(offsetToOtherParticle, offsetToOtherParticle);
                        
                        // Skip particle if its outside the smoothing radius
                        if (sqrDistanceToOtherParticle > pow(smoothingRadius, 2)) continue;
                        
                        float distance = sqrt(sqrDistanceToOtherParticle);
                        float influence = densityKernel(distance);
                        
                        density += mass * influence; 
                        
                        if (int(particles[otherParticleIndex].position.x) % int(smoothingRadius) == 0 || int(particles[otherParticleIndex].position.y) % int(smoothingRadius) == 0)
                        {
                            particles[otherParticleIndex].velocity = vec2(800, 800);
                        }   
                        //particles[gID].velocity = vec2(400, 400);
                        //particles[gID].velocity = vec2(float(cell % gridWidth), float(cell / gridWidth)) * 7.0;
                    }
                }
            }
        }   */

        /*
        for (int i = 0; i < 9; i++)
        {
            int offset = cellOffsets[i];
            int cellToCheck = cell + offset;
            
            //cellToCheck = 1598 + (gridWidth + 1);
            if (cellToCheck < gridWidth * gridHeight && cellToCheck >= 0 && startIndices[cellToCheck] >= 0)
            {
                int endIndex = cellToCheck + 1 < gridWidth * gridHeight ? startIndices[cellToCheck + 1] : particles.length();
                for (int otherParticleIndex = startIndices[cellToCheck]; otherParticleIndex < endIndex; otherParticleIndex++)
                {
                    if (oneD[otherParticleIndex] == -1) break;
                    vec2 otherParticlePosition = particles[oneD[otherParticleIndex]].position + particles[oneD[otherParticleIndex]].velocity * deltaTime;
                    vec2 offsetToOtherParticle = otherParticlePosition - currentParticlePosition;
                    float sqrDistanceToOtherParticle = dot(offsetToOtherParticle, offsetToOtherParticle);
                    
                    // Skip particle if its outside the smoothing radius
                    if (sqrDistanceToOtherParticle > pow(smoothingRadius, 2)) continue;
                    
                    float distance = sqrt(sqrDistanceToOtherParticle);
                    float influence = densityKernel(distance);
                    
                    density += mass * influence; 
                    
                    //particles[otherParticleIndex].velocity = vec2(100 + i * 77.77, 100 + i * 77.77);
                    particles[gID].velocity = vec2(float(cell % gridWidth), float(cell / gridWidth)) * 7.0;
                }
            }
        }
        */
        // Add calculated density for the particle to it's density
        particles[gID].density = density;
    }
)glsl";