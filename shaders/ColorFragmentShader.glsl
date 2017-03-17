#version 330 core
// Interpolated values from the vertex shaders
in vec2 UV;
in vec3 Position_worldspace;
in vec3 Normal_cameraspace;
in vec3 LightDirection_cameraspace;
//in vec3 EyeDirection_cameraspace;
in vec4 ShadowCoord;

// Ouput data
//out vec3 color;
layout(location = 0) out vec3 color;

// Values that stay constant for the whole mesh.
//uniform sampler2D myTextureSampler;
uniform vec3 LightPosition_worldspace;
uniform vec3 LightColor;
uniform sampler2DShadow shadowMap;
uniform sampler2D textureSampler;
uniform vec4 limbDarkingCoeffs;


void main(){
    vec3 MaterialDiffuseColor = vec3(1.0, 1.0, 1.0);
    vec3 MaterialAmbientColor = texture( textureSampler, UV ).rgb;;

    float distance = length( LightPosition_worldspace - Position_worldspace );

    float bias = 0.1;
    float visibility = 0;
    visibility = texture( shadowMap,  vec3(ShadowCoord.xy, ShadowCoord.z-bias)/ShadowCoord.w );
    visibility = float(distance < length(LightPosition_worldspace) || visibility > 0);

    vec3 n = normalize( Normal_cameraspace );
    vec3 l = normalize( LightDirection_cameraspace );

    float cosTheta = clamp( dot( n,l ), 0, 1 );
    float cosAlpha = clamp( n.z, 0, 1 );

    float limbDarking = 1 +
                        limbDarkingCoeffs[0] * (1 - cosAlpha) +
                        limbDarkingCoeffs[1] * (1 - cosAlpha*cosAlpha) +
                        limbDarkingCoeffs[2] * (1 - pow(cosAlpha, 3)) +
                        limbDarkingCoeffs[3] * (1 - pow(cosAlpha, 4));

    color = (
    	    MaterialAmbientColor +
    	    MaterialDiffuseColor * visibility * LightColor * cosTheta / (distance*distance)
    	) * limbDarking;
//    color = vec3(0,0,limbDarking);
}
