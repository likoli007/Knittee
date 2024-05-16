#version 330 core

layout(location = 0) in vec3 position;

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
        fragColor = vec4(0.5450, 0.5450, 0.5450, 1.0);
    }



}