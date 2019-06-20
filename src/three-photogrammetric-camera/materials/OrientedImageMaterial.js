import { Uniform, ShaderMaterial, ShaderLib, ShaderChunk, Matrix4, Vector3, Vector4, Color } from 'three';
import { default as RadialDistortion } from '../cameras/distortions/RadialDistortion';

function pop(options, property, defaultValue) {
    if (options[property] === undefined) return defaultValue;
    const value = options[property];
    delete options[property];
    return value;
}

function popUniform(options, property, defaultValue) {
    const value = pop(options, property, defaultValue);
    if (options.uniforms[property])
        return options.uniforms[property];
    return new Uniform(value);
}

ShaderChunk.common = `${ShaderChunk.common}
#ifdef USE_WORLDPOS
varying vec3 vEyePosition;
varying float vDepth;
#endif
#ifdef USE_MAP4
#undef USE_MAP
  uniform vec3 uvwPosition;
#endif
`;

ShaderChunk.worldpos_vertex = `
#if defined( USE_WORLDPOS ) || defined( USE_ENVMAP ) || defined( DISTANCE ) || defined ( USE_SHADOWMAP )
  vec4 worldPosition = modelMatrix * vec4( transformed, 1.0 );
#endif
#ifdef USE_WORLDPOS
  mat4 m = modelMatrix;
  m[3].xyz -= uvwPosition;
  vEyePosition = (m * vec4( transformed, 1.0 )).xyz;
  //vDepth = distance(cameraPosition, worldPosition.xyz);
#endif

`;

ShaderChunk.color_pars_fragment = `${ShaderChunk.color_pars_fragment}
${RadialDistortion.chunks.radial_pars_fragment}
uniform bool diffuseColorGrey;
uniform sampler2D tDepth;
uniform float foggy;
uniform float near;
uniform float far;
uniform float parallax;

#ifdef USE_MAP4
  uniform mat4 uvwPreTransform;
  uniform mat4 uvwPostTransform;
  uniform RadialDistortion uvDistortion;
  uniform sampler2D map;
  uniform float borderSharpness;

#endif

const float GOLDEN_ANGLE = 2.39996323;
const float MAX_BLUR_SIZE = 20.0;
const float RAD_SCALE = 0.5; // Smaller = nicer blur, larger = faster
vec2 uPixelSize = vec2(1.0/1920., 1.0/1080.);
float uFar = far; 

/*
float getBlurSize(float depth, float focusPoint, float focusScale)
{
 float coc = clamp((1.0 / focusPoint - 1.0 / depth)*focusScale, -1.0, 1.0);
 return abs(coc) * MAX_BLUR_SIZE;
}


vec3 depthOfField(vec2 texCoord, float focusPoint, float focusScale)
{
 float centerDepth = texture2D(tDepth, texCoord).r * uFar;
 float centerSize = getBlurSize(centerDepth, focusPoint, focusScale);
 vec3 color = texture2D(map, vTexCoord).rgb;
 float tot = 1.0;

 
 float radius = RAD_SCALE;
 for (float ang = 0.0; radius<MAX_BLUR_SIZE; ang += GOLDEN_ANGLE)
 {
  vec2 tc = texCoord + vec2(cos(ang), sin(ang)) * uPixelSize * radius;

  vec3 sampleColor = texture2D(map, tc).rgb;
  float sampleDepth = texture2D(tDepth, tc).r * uFar;
  float sampleSize = getBlurSize(sampleDepth, focusPoint, focusScale);
  if (sampleDepth > centerDepth)
   sampleSize = clamp(sampleSize, 0.0, centerSize*2.0);

  float m = smoothstep(radius-0.5, radius+0.5, sampleSize);
  color += mix(color/tot, sampleColor, m);
  tot += 1.0;
  radius += RAD_SCALE/radius;
 }
 
 return color /= tot;
}
*/


vec3 averageColor(sampler2D tex, vec2 coord, float dist){

  return (texture2D(tex, vec2(coord.x + dist, coord.y)) +
          texture2D(tex, vec2(coord.x, coord.y + dist)) + 
          texture2D(tex, vec2(coord.x - dist, coord.y)) + 
          texture2D(tex, vec2(coord.x, coord.y - dist))).rgb / 4.;
}

`;

ShaderChunk.color_fragment = `${ShaderChunk.color_fragment}


  if (diffuseColorGrey) {
    diffuseColor.rgb = vec3(dot(diffuseColor.rgb, vec3(0.333333)));
  }
#ifdef USE_MAP4
    vec4 uvw = uvwPreTransform * vec4(vEyePosition, 1.);
  distort_radial(uvw, uvDistortion);
  uvw = uvwPostTransform * uvw;
  uvw.xyz /= 2. * uvw.w;
  uvw.xyz += vec3(0.5);
  // diffuseColor.rgb = vec3(fract(uvw.xyz));
  vec3 border = min(uvw.xyz, 1. - uvw.xyz);
//  if (all(greaterThan(border,vec3(0.))))
///  {
    vec4 color = texture2D(map, uvw.xy);
    color.a *= min(1., borderSharpness*min(border.x, border.y));
    diffuseColor.rgb = mix(diffuseColor.rgb, color.rgb, color.a);
//  }
 // diffuseColor.rgb = vec3(vDepth / 1000.);


 
        // Alex modif
				float fragCoordZ = texture2D( tDepth, uvw.xy ).x;
        float viewZ = ( near * far ) / ( ( far - near ) * fragCoordZ - far );
				float d = ( viewZ + near ) / ( near - far);
        d = 1. - d;

        if(parallax > 0.){
          diffuseColor.rgb =  texture2D(map, vec2(uvw.x + d * parallax / 10., uvw.y ) ).rgb;
        }
      
        if(foggy > 0.){
          float fogFactor = d; // - foggy;
          vec3 blurry = (averageColor(map, uvw.xy, foggy / 50.) + averageColor(map, uvw.xy, foggy / 65.) + averageColor(map, uvw.xy, foggy / 80.) + averageColor(map, uvw.xy, foggy / 100.)) / 4.; 
         // diffuseColor.rgb = mix(diffuseColor.rgb, blurry, 1. - fogFactor);
          float fogFactorX = 1. / exp(pow(d * 5., 2.));
          diffuseColor.rgb = mix(diffuseColor.rgb, vec3( 1.) , foggy * ( fogFactorX )); // mix(diffuseColor.rgb, vec3( 1.) , foggy * (1. - fogFactor)); 
        }

        



#endif
 
`;

function definePropertyUniform(object, property, defaultValue) {
    object.uniforms[property] = new Uniform(object[property] || defaultValue);
    Object.defineProperty(object, property, {
        get: () => object.uniforms[property].value,
        set: (value) => {
            if (object.uniforms[property].value != value) {
                object.uniformsNeedUpdate = true;
                object.uniforms[property].value = value;
            }
        }
    });
}

class OrientedImageMaterial extends ShaderMaterial {
    constructor(options = {}) {
        const size = pop(options, 'size', 1);
        const diffuse = pop(options, 'diffuse', new Color(0xeeeeee));
        const uvwPosition = pop(options, 'uvwPosition', new Vector3());
        const uvwPreTransform = pop(options, 'uvwPreTransform', new Matrix4());
        const uvwPostTransform = pop(options, 'uvwPostTransform', new Matrix4());
        const uvDistortion = pop(options, 'uvDistortion', {R: new Vector4(), C: new Vector3()});
        const map = pop(options, 'map', null);
        const tDepth = pop(options, 'tDepth', new Color(0xeeeeee));
        const alphaMap = pop(options, 'alphaMap', null);
        const scale = pop(options, 'scale', 1);
        const borderSharpness = pop(options, 'borderSharpness', 100);
        const diffuseColorGrey = pop(options, 'diffuseColorGrey', true);
        options.vertexShader = options.vertexShader || ShaderLib.points.vertexShader;
        options.fragmentShader = options.fragmentShader || ShaderLib.points.fragmentShader;
        options.defines = options.defines || {};
        if (map) {
            options.defines.USE_MAP4 = '';
            options.defines.USE_WORLDPOS = '';
        }
        if (alphaMap) options.defines.USE_ALPHAMAP = '';
        if (options.vertexColors) options.defines.USE_COLOR = '';
        if (options.logarithmicDepthBuffer) options.defines.USE_LOGDEPTHBUF = '';
        if (pop(options, 'sizeAttenuation')) options.defines.USE_SIZEATTENUATION = '';
        super(options);
        definePropertyUniform(this, 'size', size);
        definePropertyUniform(this, 'diffuse', diffuse);
        definePropertyUniform(this, 'tDepth', tDepth);
        definePropertyUniform(this, 'uvwPosition', uvwPosition);
        definePropertyUniform(this, 'uvwPreTransform', uvwPreTransform);
        definePropertyUniform(this, 'uvwPostTransform', uvwPostTransform);
        definePropertyUniform(this, 'uvDistortion', uvDistortion);
        definePropertyUniform(this, 'opacity', this.opacity);
        definePropertyUniform(this, 'map', map);
        definePropertyUniform(this, 'alphaMap', alphaMap);
        definePropertyUniform(this, 'scale', scale);
        definePropertyUniform(this, 'foggy', 0.);
        definePropertyUniform(this, 'near', 10.);
        definePropertyUniform(this, 'far', 500.);
        definePropertyUniform(this, 'parallax', 0.);
        definePropertyUniform(this, 'borderSharpness', borderSharpness);
        definePropertyUniform(this, 'diffuseColorGrey', diffuseColorGrey);
    }
}

export default OrientedImageMaterial;