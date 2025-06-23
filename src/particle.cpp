#include <utility>

class Particle
{
private:


public:
    std::pair<float, float> position;
    std::pair<float, float> velocity;
    Particle(std::pair<float, float> position, std::pair<float, float> velocity) : 
        position(position), velocity(velocity) {}
};