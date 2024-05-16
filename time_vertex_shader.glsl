#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec2 texCoord; 


out vec2 TexCoord;

uniform mat4 mvpMatrix;
uniform int selectedFace;


void main()
{
    gl_Position = mvpMatrix * vec4(position, 1.0);

    TexCoord = texCoord;
}