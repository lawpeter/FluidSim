const char* fragmentShaderSource = R"glsl(
    #version 330 core

    in float particleSpeed;
    in vec2 localPosition;

    out vec4 FragColor;

    void main() {
        float distance = dot(localPosition, localPosition);
        float alpha = smoothstep(1.0, 0.9, distance);

        if (alpha < 0.005) discard;

        vec4 deepBlue = vec4(0.14, 0.14, 0.74, 1.0);
        vec4 blue = vec4(0.35, 0.56, 0.90, 1.0);
        vec4 yellow = vec4(1.0, 0.91, 0.14, 1.0);
        vec4 red = vec4(0.87, 0.24, 0.24, 1.0);
        vec4 white = vec4(1.0, 1.0, 1.0, 1.0);

        if (particleSpeed < 0.5) FragColor = mix(blue, yellow, 2.0 * particleSpeed);
        else FragColor = mix(yellow, red, 2.0 * (particleSpeed - 0.5f));

        vec4 baseColor = mix(deepBlue, white, particleSpeed);
        //FragColor = vec4(baseColor.rgb, alpha);
    }
)glsl";