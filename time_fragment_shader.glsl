#version 330

uniform sampler2D tex;
in vec2 TexCoord;

out vec4 fragColor;

void main() {
    fragColor = vec4(texture(tex, TexCoord).rgb, 1.0);
}