// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'

Shader "Unlit/SteerablePerlinNoise"
{
    Properties
    {
        _UseNoise3d("Use 3D Noise",  Int) = 0
        _EigenValueSum("Eigenvalue Sum", float) = 4        
        _NoiseFrequency("Noise Frequency", Vector) = (1.,1.,1.)
        _NoiseOffset("Noise Offset", Vector) = (0,0,0)
        _NumOctaves("Number of Octaves", Int) = 5 
        _OctaveBias("Noise Bias", float) = .65      
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            // make fog work
            #pragma multi_compile_fog

            #include "UnityCG.cginc"
            #define PI 3.1415926535
            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
                float3 normal : NORMAL;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float3 coords : TEXCOORD1;
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
                float3 normal : TEXCOORD2;
            };

            int _UseNoise3d;
            float3 _NoiseFrequency;
            float3 _NoiseOffset;
            int _NumOctaves;
            float _OctaveBias;
            float _EigenValueSum;
            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.coords = mul (unity_ObjectToWorld, v.vertex);// v.vertex;
                UNITY_TRANSFER_FOG(o,o.vertex);
                o.normal =  normalize( mul( float4( v.normal, 0.0 ), unity_WorldToObject ).xyz );;
                return o;
            }
            float random3(float3 pos){
                return frac(sin(dot(pos, float3(64.25375463f, 23.27536534f, 86.29678483f))) * 59482.7542f);
            }

            float3 random33(float3 pos){
                return float3(random3(pos + .01f), random3(pos + .02f), random3(pos + .03f));
            }

            float3 rsphere(float3 p) {
                float3 vals = random33(p);
                
                float u = vals.x;
                float v = vals.y;
                float theta = u * 2.0f * PI;
                float phi = acos(2.0f * v - 1.0f);
                float r = pow(vals.z, .33333);
                float sinTheta = sin(theta);
                float cosTheta = cos(theta);
                float sinPhi = sin(phi);
                float cosPhi = cos(phi);
                float x = r * sinPhi * cosTheta;
                float y = r * sinPhi * sinTheta;
                float z = r * cosPhi;
                return float3(x,y,z);
            }

            float smootherstep(float x){
                return clamp(6 * (x * x * x * x * x) - 15 * (x * x * x * x) + 10 * (x * x * x), 0, 1);
            }
            
            float interp(float u){
                return 1.0 - smootherstep(abs(u));//smoothstep(0.0, 1.0, abs(u));
            }

            float3 interp3(float3 u){
                return float3(interp(u.x), interp(u.y), interp(u.z));
            }

            float fitrange(float x, float in_low, float in_high, float out_low, float out_high){
                float u = clamp(x, in_low, in_high);
                return ((out_high - out_low) * (u - in_low)) / (in_high - in_low) + out_low;
            }
            float2 fitrange_2(float2 x, float in_low, float in_high, float out_low, float out_high){
                return float2(fitrange(x.x, in_low, in_high, out_low, out_high),
                            fitrange(x.y, in_low, in_high, out_low, out_high));
            }
            float3 fitrange_3(float3 x, float in_low, float in_high, float out_low, float out_high){
                return float3(fitrange(x.x, in_low, in_high, out_low, out_high),
                            fitrange(x.y, in_low, in_high, out_low, out_high),
                            fitrange(x.z, in_low, in_high, out_low, out_high));
            }
            float3x3 outerprod3(float3 x, float3 y){
                return float3x3(x.x * y.x, x.y * y.x, x.z * y.x, 
                            x.x * y.y, x.y * y.y, x.z * y.y,
                            x.x * y.z, x.y * y.z, x.z * y.z);
            }
            float2x2 outerprod2(float2 x, float2 y){
                return float2x2(x.x * y.x, x.y * y.x,
                                x.x * y.y, x.y * y.y);

            }
            float3x3 generate_metric(float3 _x){
                float3 x = normalize(_x);
                
                //hack to avoid division by zero when extracting the eigenvectors :)
                //without this i think there'd need to be a bunch of if statements for all the different cases
                //like if x = (.5,.5,0) or x = (0,.5,.5) or x = (0,0,1), etc.
                x = max(abs(x),1e-5) * (x >= 0.0 ? 1 : -1);
                
                float3 evals = float3(0,0,1.0);
                float3 evec0, evec1, evec2;
                //eigenvectors of outerproduct(x, x) have a closed form solution:

                evec0 = normalize(float3(-x.z / x.x,0,1));
                evec1 = normalize(float3(-x.y / x.x,1.0,0.0));
                evec2 = normalize(float3(x.x/x.z, x.y/x.z, 1.0));    
                
                float k = .51;
                int dimensions = 3;
                float denom = 1. / float(dimensions - 1);
                float3 mapped_evals = fitrange_3(evals, 0.0, 1.0,  (_EigenValueSum - k) * denom - 1e-4, k);
                //float3 mapped_evals = float3(1.5, 1.5, 1.);
                return  mapped_evals.x * outerprod3(evec0, evec0) +
                        mapped_evals.y * outerprod3(evec1, evec1) + 
                        mapped_evals.z * outerprod3(evec2, evec2);
            }

            float2x2 generate_metric(float2 _x){
                float2 x = normalize(_x);
                
                //hack to avoid division by zero when extracting the eigenvectors :)
                //without this i think there'd need to be a bunch of if statements for all the different cases
                //like if x = (.5,.5,0) or x = (0,.5,.5) or x = (0,0,1), etc.
                x = max(abs(x),1e-5) * (x >= 0.0 ? 1 : -1);
                
                float2 evals = float2(0,1.0);
                float2 evec0, evec1;
                //eigenvectors of outerproduct(x, x) have a closed form solution:

                evec0 = normalize(float2(-x.y / x.x,1.0));
                evec1 = normalize(float2(x.x / x.y,1.0));
                
                float k = .51;
                int dimensions = 2;
                float denom = 1. / float(dimensions - 1);
                float2 mapped_evals = fitrange_2(evals, 0.0, 1.0,  (_EigenValueSum - k) * denom - 1e-4, k);
                //float3 mapped_evals = float3(1.5, 1.5, 1.);
                return  mapped_evals.x * outerprod2(evec0, evec0) +
                        mapped_evals.y * outerprod2(evec1, evec1);
            }            
            float3x3 vec_projector(float3 N){
                return float3x3(1,0,0,0,1,0,0,0,1) - outerprod3(N,N);
            }
            float3x3 dihedral(float3 a, float3 b){
                float3 v = cross(a,b);
                float c = dot(a,b);
                float3x3 I = float3x3(1,0,0,0,1,0,0,0,1);
                float3x3 v_prime;
                v_prime[0] = float3(0, -v.z, v.y);
                v_prime[1] = float3(v.z, 0, -v.x);
                v_prime[2] = float3(-v.y, v.x, 0);
                float scale = 1.0 / (1.f + c);
                float3x3 v_prime_2 = mul(v_prime, v_prime);
                return I + v_prime + v_prime_2 * scale;
            }

            float steerable_perlin(float3 pos, float3x3  metric){
                float3 noise_p = floor(pos);
                float3 noise_f = pos - noise_p;//frac(pos);
                
                float out_val = 0.0;
                
                //perlin weights
                float3 blend = interp3(noise_f);
                
                //we opt to use a weighted average instead of lerps.
                //if we remove the anisotropy from the dot product and the interpolation
                //the weighted average results in the same output as standard perlin noise.
                for(int i = 0; i <= 1; i++)
                for(int j = 0; j <= 1; j++)
                for(int k = 0; k <= 1; k++){

                    float3 o = float3(i,j,k);
                    
                    float3 g = noise_p  + o;
                    float3 r = rsphere(g);
                    float3 v = (o - noise_f);
                    
                    float3 metric_v = mul(metric, v);
                    
                    //regular perlin is dot(r, v)
                    //applying the metric to v, adds one level of anisotropy
                    float d = dot(r, metric_v);
                    
                    float3 wv = abs(o-blend);
                    
                    //we get another level of anisotropy by introducing anisotropic weights to the algorithm       
                    float w =  wv.x * wv.y * wv.z * interp(dot(v, metric_v));

                    out_val += d * w;
        
                } 

                return out_val;
            }
            float steerable_perlin_projected(float3 pos, float2x2  metric, float3x3 projection){
                float3 noise_p = floor(pos);
                float3 noise_f = pos - noise_p;//frac(pos);
                
                float out_val = 0.0;
                
                //perlin weights
                float3 blend = interp3(noise_f);

                //we opt to use a weighted average instead of lerps.
                //if we remove the anisotropy from the dot product and the interpolation
                //the weighted average results in the same output as standard perlin noise.
                for(int i = 0; i <= 1; i++)
                for(int j = 0; j <= 1; j++)
                for(int k = 0; k <= 1; k++){

                    float3 o = float3(i,j,k);
                    
                    float3 g = noise_p  + o;
                    float2 r = mul(projection, rsphere(g)).xy;
                    float2 v = mul(projection, o - noise_f).xy;
                    
                    float2 metric_v = mul(metric, v);
                    
                    //regular perlin is dot(r, v)
                    //applying the metric to v, adds one level of anisotropy
                    float d = dot(r, metric_v);
                    
                    float3 wv = abs(o-blend);
                    
                    //we get another level of anisotropy by introducing anisotropic weights to the algorithm       
                    float w =  wv.x * wv.y * wv.z * interp(dot(v, metric_v));
                    out_val += d * w;
        
                } 

                return out_val;
            }

           
            float fbm_3d(float3 p, float3x3 metric)
            {
                float bias = _OctaveBias;
                float freq = _NoiseFrequency;
                int octaves = max(1,_NumOctaves);
                float3 offset = _NoiseOffset;
                float out_val  = 0;

                for(int i = 0; i < octaves; i++){
                    out_val += pow(bias, float(i)) * steerable_perlin(pow(2.0, float(i)) * (p + offset) * freq , metric);
                }
                        
                return out_val;
            }
            float fbm_projected(float3 p, float2x2 metric, float3x3 projection)
            {
                float bias = _OctaveBias;
                float freq = _NoiseFrequency;
                int octaves = max(1,_NumOctaves);
                float3 offset = _NoiseOffset;
                float out_val  = 0;

                for(int i = 0; i < octaves; i++){
                    out_val += pow(bias, float(i)) * steerable_perlin_projected(pow(2.0, float(i)) * (p + offset) * freq , metric, projection);
                }
                        
                return out_val;
            }
            //this is really only useful if you want a mesh to move through the noise field.
            float fbm_3d_artifact_free(float3 p, float3x3 metric)
            {
                float bias = _OctaveBias;
                float freq = _NoiseFrequency;
                int octaves = max(1,_NumOctaves);
                float3 offset = _NoiseOffset;
                float out_val  = 0;

                for(int i = 0; i < octaves; i++){
                    //since the weights are always heighest at the .5 position, combine two noises at the same octave to remove artifacts.
                    out_val += pow(bias, float(i)) * (steerable_perlin(pow(2.0, float(i)) * (p + offset) * freq , metric) + steerable_perlin(pow(2.0, float(i)) * ((p+offset) * freq + float3(.5,.5,.5)), metric)) * .5;
                }
                        
                return out_val;
            }


            fixed4 frag (v2f i) : SV_Target
            {
                float col;
                //this is the anisotropy direction to steer the noise!
                float3 anisotropy_dir =i.coords.xyz * float3(1,-1,-1) + float3(0,0,sin(i.coords.x * i.coords.y));

                //3d version
                if(_UseNoise3d){
                    float3x3 x = generate_metric(anisotropy_dir);
                    col = fbm_3d(i.coords.xyz, x);
                }else{

                    float3x3 projector = dihedral(normalize(i.normal), float3(0,0,1));
                    float3x3 to_normal = vec_projector(i.normal);
                    projector = mul(projector, to_normal);

                    float2x2 x = generate_metric(mul(projector, anisotropy_dir).xy);
                    col = fbm_projected(i.coords.xyz, x, projector);                    
                }
                return col;
            }
            ENDCG
        }
    }
}
