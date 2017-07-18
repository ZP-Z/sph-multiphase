#version 400

layout (points) in;
layout (triangle_strip,max_vertices=4) out;

uniform mat4 ProjectionMatrix;

uniform float particle_size;

in vec4 ex_Color[];
in mat4 LocalRotation[];

out vec2 Vertex_UV;
out vec4 int_Color;

void generateCube(vec4 P){
  //up
  vec3 displacement =  vec3(-0.5,0.5,0.5);
  vec4 rotatedDx = LocalRotation[0] * vec4(displacement, 1); 
  gl_Position = P+rotatedDx;
  int_Color=vec4(1.0,0,0,1);
  EmitVertex();

  displacement = vec3(-0.5,0.5,-0.5);
  rotatedDx=LocalRotation[0]*vec4(displacement,1);
  gl_Position = P+rotatedDx;
  int_Color=vec4(1.0,0,0,1);
  EmitVertex();

  displacement = vec3(0.5,0.5,0.5);
  rotatedDx=LocalRotation[0]*vec4(displacement,1);
  gl_Position = P+rotatedDx;
  int_Color=vec4(1.0,0,0,1);
  EmitVertex();

  displacement = vec3(0.5,0.5,-0.5);
  rotatedDx=LocalRotation[0]*vec4(displacement,1);
  gl_Position = P+rotatedDx;
  int_Color=vec4(1.0,0,0,1);
  EmitVertex();

  EndPrimitive();
  //down

  //left

  //right

  //front

  //back
}

void main(void)
{

  vec4 P = gl_in[0].gl_Position;

  generateCube(P);
  

//  // a: left-bottom
//
//  vec2 va = P.xy + vec2(-0.5,-0.5)*particle_size;
//  gl_Position = ProjectionMatrix * vec4(va,P.zw);
//  Vertex_UV = vec2(0.0, 0.0);
//  int_Color = ex_Color[0];
//  EmitVertex();
//
//  // b: left-top
//  vec2 vb = P.xy + vec2(-0.5,0.5)*particle_size;
//  gl_Position = ProjectionMatrix * vec4(vb,P.zw);
//  Vertex_UV = vec2(0.0, 1.0);
//  int_Color = ex_Color[0];
//  EmitVertex();
//
//  // c: right-bottom
//  vec2 vd = P.xy + vec2(0.5,-0.5)*particle_size;
//  gl_Position = ProjectionMatrix * vec4(vd,P.zw);
//  Vertex_UV = vec2(1.0, 0.0);
//  int_Color = ex_Color[0];
//  EmitVertex();
//
//  // d: right-top
//  vec2 vc = P.xy + vec2(0.5,0.5)*particle_size;
//  gl_Position = ProjectionMatrix * vec4(vc,P.zw);
//  Vertex_UV = vec2(1.0, 1.0);
//  int_Color = ex_Color[0];
//  EmitVertex();
  

  EndPrimitive();
}