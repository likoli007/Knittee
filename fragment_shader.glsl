#version 330 core

in vec4 fragColor;
in float shading;

out vec4 resultColor;

uniform float shadingValue;

void main()
{
    resultColor = fragColor*(shadingValue+0.4);
    

}