#include <utility>
#include <iostream>
#include <glm/glm.hpp>

#define SCREEN_WIDTH 1000.0f
#define SCREEN_HEIGHT 800.0f

class Particle 
{
private:


public:
    glm::vec2 position;
    glm::vec2 velocity;
    static constexpr float radius = 5.0f;
    static float restitution;
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

        if (position.x - radius < 0.0f)
        {
            position.x = radius;
            velocity.x *= -restitution;
        }
        else if (position.x + radius > SCREEN_WIDTH)
        {
            position.x = SCREEN_WIDTH - radius;
            velocity.x *= -restitution;
        }
        if (position.y - radius < 0.0f)
        {
            position.y = radius;
            velocity.y *= -restitution;           
        }
        else if (position.y + radius > SCREEN_HEIGHT)
        {
            position.y = SCREEN_HEIGHT - radius;
            velocity.y *= -restitution;
        }
    }
};