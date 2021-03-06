#version 330

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec3 color;
layout(location = 3) in vec3 texCoord;

// varying variables
out vec3 varColor;
out vec3 varTexCoord;
out vec3 varNormal;
out vec3 varPosition;

// uniform variables
uniform int mode;	// 1 -- color / 2 -- texture
uniform mat4 mvpMatrix;
uniform mat4 mvMatrix;

uniform int shadowState;	// 1 -- normal / 2 -- shadow
uniform mat4 light_mvpMatrix;
uniform vec3 lightDir;

void main(){
	varColor = color;
	varTexCoord = texCoord;
	varPosition = position;

	// SHADOW: From light
	if (shadowState == 2) {
		gl_Position = light_mvpMatrix * vec4(varPosition, 1.0);
		return;
	}
	
	varNormal = normalize(normal).xyz;

	gl_Position = mvpMatrix * vec4(position, 1.0);
}