#version 330 core

in vec4 fragColor;
in float shading;

out vec4 resultColor;

uniform float shadingValue;
uniform int selectedFace;

void main()
{
    

    if (selectedFace == 0){
        resultColor = fragColor*(shadingValue+0.5);
    }
    else{
        resultColor = fragColor;
    }
    
    

}