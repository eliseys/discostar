#version 330 core

precision highp float;
precision highp vec2;
precision highp vec3;
precision highp vec4;

in highp vec2 VertexUV;

// Input vertex data, different for all executions of this shader.
layout(location = 0) in highp vec3 vertexPosition_modelspace;

// Output data ; will be interpolated for each fragment.
out highp vec2 UV;

void main(){
	gl_Position =  vec4(vertexPosition_modelspace,1);
//   	UV = (vertexPosition_modelspace.xy+vec2(1,1)) / 2.0;
    UV = VertexUV;
}