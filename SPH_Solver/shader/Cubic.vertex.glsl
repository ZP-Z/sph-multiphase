#version 400

layout(location=0) in vec3 in_Position;
layout(location=1) in vec4 in_Color;
layout(location=2) in mat4 in_Rotation;

out vec4 ex_Color;
out mat4 ex_Rotation;

uniform mat4 ModelMatrix;
uniform mat4 ViewMatrix;

void main(void)
{
  gl_Position = ViewMatrix * vec4(in_Position.xyz,1.0);
  ex_Color = in_Color;
  ex_Rotation = in_Rotation;
}