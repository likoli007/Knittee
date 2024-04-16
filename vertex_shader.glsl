#version 330 core

layout(location = 0) in vec3 position;
//layout(location = 1) in float shadingValue; // Assuming shadingValue is the attribute name

out vec4 fragColor;

uniform mat4 mvpMatrix;
uniform int selectedFace;

void main()
{
    gl_Position = mvpMatrix * vec4(position, 1.0);
    
    if(selectedFace == 1){
        fragColor = vec4(0.0, 1.0, 0, 1.0);
    }
    else{
        fragColor = vec4(1.0, 0, 0, 1.0);
    }



}