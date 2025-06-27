#include <utility>

#ifndef DIMENSIONS

#define DIMENSIONS
#define SCREEN_WIDTH 800.0f
#define SCREEN_HEIGHT 800.0f

#endif
#ifndef INCLUDE

#define INCLUDE
#include <iostream>
#include <glm/glm.hpp>

#endif

class Particle 
{
private:


public:
    glm::vec2 position;
    glm::vec2 velocity;
    Particle(glm::vec2 position, glm::vec2 velocity) : 
        position(position), velocity(velocity) {}
    void accelerate(glm::vec2 acceleration)
    {
        velocity.x += acceleration.x;
        velocity.y += acceleration.y;
    }
    void updatePosition(float deltaTime)
    {
        position.x += velocity.x * deltaTime;
        position.y += velocity.y * deltaTime;

        if (position.x - 10.0f < 0.0f)
        {
            position.x = 10.0f;
            velocity.x *= -1.0f;
        }
        else if (position.x + 10.0f > SCREEN_WIDTH)
        {
            position.x = SCREEN_WIDTH - 10.0f;
            velocity.x *= -1.0f;
        }
        if (position.y - 10.0f < 0.0f)
        {
            position.y = 10.0f;
            velocity.y *= -1.0f;           
        }
        else if (position.y + 10.0f > SCREEN_HEIGHT)
        {
            position.y = SCREEN_HEIGHT - 10.0f;
            velocity.y *= -1.0f;
        }


    }
};