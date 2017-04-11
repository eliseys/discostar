#version 330 core

precision highp float;
precision highp vec2;
precision highp vec3;
precision highp vec4;
precision highp sampler2D;
precision highp sampler2DShadow;

// Interpolated values from the vertex shaders
in highp vec2 UV;
in highp vec3 Position_worldspace;
in highp vec3 Normal_cameraspace;
in highp vec3 LightDirection_cameraspace;
//in vec3 EyeDirection_cameraspace;
in highp vec4 ShadowCoord;

// Ouput data
//out vec3 color;
layout(location = 0) out highp vec3 color;

// Values that stay constant for the whole mesh.
//uniform sampler2D myTextureSampler;
uniform highp vec3 LightPosition_worldspace;
uniform highp vec3 LightColor;
uniform highp sampler2DShadow shadowMap;
uniform highp sampler2D textureSampler;
uniform highp vec4 limbDarkingCoeffs;


void main(){
    highp vec3 MaterialDiffuseColor = vec3(1.0, 1.0, 1.0);
    highp vec3 MaterialAmbientColor = texture( textureSampler, UV ).rgb;;

    highp float distance = length( LightPosition_worldspace - Position_worldspace );

    highp float bias = 0.1;
    highp float visibility = 0;
    visibility = texture( shadowMap,  vec3(ShadowCoord.xy, ShadowCoord.z-bias)/ShadowCoord.w );
    visibility = float(distance < length(LightPosition_worldspace) || visibility > 0);

    highp vec3 n = normalize( Normal_cameraspace );
    highp vec3 l = normalize( LightDirection_cameraspace );

    highp float cosTheta = clamp( dot( n,l ), 0, 1 );
    highp float cosAlpha = clamp( n.z, 0, 1 );

    highp float limbDarking = 1 +
                        limbDarkingCoeffs[0] * (1 - cosAlpha) +
                        limbDarkingCoeffs[1] * (1 - cosAlpha*cosAlpha) +
                        limbDarkingCoeffs[2] * (1 - pow(cosAlpha, 3)) +
                        limbDarkingCoeffs[3] * (1 - pow(cosAlpha, 4));

    color = (
    	    MaterialAmbientColor +
    	    MaterialDiffuseColor * visibility * LightColor * cosTheta / (12.566370614359172 * distance*distance)
    	) * limbDarking;
    color = sqrt(sqrt(color));
}
