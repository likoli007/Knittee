#version 330 core

layout (triangles) in;
layout (triangle_strip, max_vertices = 9) out;

out vec3 Color; // Output variable for passing color to fragment shader

void main()
{
    for (int i = 0; i < gl_in.length(); i++) {
        gl_Position = gl_in[i].gl_Position;
        Color = vec3(1.0, 0.0, 0.0); // Set color to red for the outline
        EmitVertex();
    }
    EndPrimitive();

    // Emit vertices for the interior of the triangle
    for (int i = 0; i < 3; i++) {
        gl_Position = gl_in[i].gl_Position;
        Color = vec3(0.0, 1.0, 0.0); // Set color to green for the interior
        EmitVertex();
    }
    EndPrimitive();
}