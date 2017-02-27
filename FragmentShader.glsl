#version 330 core
// Interpolated values from the vertex shaders
in vec2 UV;
in vec3 Position_worldspace;
in vec3 Normal_cameraspace;
in vec3 LightDirection_cameraspace;
in vec3 EyeDirection_cameraspace;

// Ouput data
//out vec3 color;
layout(location = 0) out vec3 color;

// Values that stay constant for the whole mesh.
//uniform sampler2D myTextureSampler;
uniform vec3 LightPosition_worldspace;
uniform vec3 LightColor;

void main(){
    // Material properties
//    vec3 MaterialDiffuseColor = texture( myTextureSampler, UV ).rgb;
    vec3 MaterialDiffuseColor = vec3(1.0, 0.0, 0.0);
    vec3 MaterialAmbientColor = vec3(0.5, 0.5, 0.5) * MaterialDiffuseColor;
//    vec3 MaterialSpecularColor = vec3(0.0,0.0,0.0);

    // Distance to the light
    float distance = length( LightPosition_worldspace - Position_worldspace );

	// Normal of the computed fragment, in camera space
    vec3 n = normalize( Normal_cameraspace );
    // Direction of the light (from the fragment to the light)
    vec3 l = normalize( LightDirection_cameraspace );

    // Cosine of the angle between the normal and the light direction,
    // clamped above 0
    //  - light is at the vertical of the triangle -> 1
    //  - light is perpendicular to the triangle -> 0
    //  - light is behind the triangle -> 0
    float cosTheta = clamp( dot( n,l ), 0, 1 );

    // Eye vector (towards the camera)
//    vec3 E = normalize(EyeDirection_cameraspace);
    // Direction in which the triangle reflects the light
//    vec3 R = reflect(-l,n);
    // Cosine of the angle between the Eye vector and the Reflect vector,
    // clamped to 0
    //  - Looking into the reflection -> 1
    //  - Looking elsewhere -> < 1
//    float cosAlpha = clamp( dot( E,R ), 0, 1 );

//    color = MaterialDiffuseColor * LightColor * LightPower * cosTheta / (distance*distance);
    color =
    	// Ambient : simulates indirect lighting
    	MaterialAmbientColor +
    	// Diffuse : "color" of the object
    	MaterialDiffuseColor * LightColor * cosTheta / (distance*distance);
    	// Specular : reflective highlight, like a mirror
//    	MaterialSpecularColor * LightColor * LightPower * pow(cosAlpha,5) / (distance*distance);

	// Output color = color of the texture at the specified UV
//    	color = texture( myTextureSampler, UV ).rgb;
//    	color = vec3(0.8,0.4,0.2);
}
