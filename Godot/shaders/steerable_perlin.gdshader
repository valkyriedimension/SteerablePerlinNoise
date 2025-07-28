shader_type spatial;

varying vec3 world_position;
varying vec3 normal;
uniform int Full_3D_Noise = 0;

uniform vec3 NOISE_FREQUENCY = vec3(1,1,1);
uniform vec3 NOISE_OFFSET = vec3(0,0,0);

uniform float OCTAVE_BIAS = .67;
uniform int NUMBER_OF_OCTAVES = 6;
uniform float EIGENVALUE_SUM = 4;



void vertex()
{
    world_position = (MODEL_MATRIX * vec4(VERTEX, 1.0)).xyz;
	normal = NORMAL;
	
}

#define PI 3.1415926535

float random3(vec3 pos){
	return fract(sin(dot(pos, vec3(64.25375463f, 23.27536534f, 86.29678483f))) * 59482.7542f);
}

vec3 random33(vec3 pos){
	return vec3(random3(pos + .01f), random3(pos + .02f), random3(pos + .03f));
}

vec3 rsphere(vec3 p) {
	vec3 vals = random33(p);
	
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
	return vec3(x,y,z);
}

float smootherstep(float x){
	return clamp(6.0f * (x * x * x * x * x) - 15.0f * (x * x * x * x) + 10.0f * (x * x * x), 0.0f, 1.0f);
}

float interp(float u){
	return 1.0 - smootherstep(abs(u));//smoothstep(0.0, 1.0, abs(u));
}

vec3 interp3(vec3 u){
	return vec3(interp(u.x), interp(u.y), interp(u.z));
}

float fitrange(float x, float in_low, float in_high, float out_low, float out_high){
	float u = clamp(x, in_low, in_high);
	return ((out_high - out_low) * (u - in_low)) / (in_high - in_low) + out_low;
}
vec2 fitrange_2(vec2 x, float in_low, float in_high, float out_low, float out_high){
	return vec2(fitrange(x.x, in_low, in_high, out_low, out_high),
				fitrange(x.y, in_low, in_high, out_low, out_high));
}

vec3 fitrange_3(vec3 x, float in_low, float in_high, float out_low, float out_high){
	return vec3(fitrange(x.x, in_low, in_high, out_low, out_high),
				fitrange(x.y, in_low, in_high, out_low, out_high),
				fitrange(x.z, in_low, in_high, out_low, out_high));
}
mat2 outerprod2(vec2 x, vec2 y){
	return mat2(vec2(x.x * y.x, x.y * y.x), 
				vec2(x.x * y.y, x.y * y.y));
}
mat3 outerprod3(vec3 x, vec3 y){
	return mat3(vec3(x.x * y.x, x.y * y.x, x.z * y.x), 
				vec3(x.x * y.y, x.y * y.y, x.z * y.y),
				vec3(x.x * y.z, x.y * y.z, x.z * y.z));
}
mat3 vec_projector(vec3 N){
    return mat3(1.0) - outerprod3(N, N);
    
}

mat3 dihedral(vec3 a, vec3 b){

    vec3 v = cross(a, b);
    float c = dot(a, b);
    mat3 v_prime;
    mat3 I = mat3(1.0);
    v_prime[0] = vec3(0.f,-v.z, v.y);
    v_prime[1] = vec3(v.z, 0, -v.x);
    v_prime[2] = vec3(-v.y, v.x, 0);
    
    float scale= 1.f / (1.f + c);
    mat3 v_prime_2 = v_prime * v_prime;

    return I + v_prime + v_prime_2 * scale;
   
}
mat3 make_projection(vec3 N){
	return vec_projector(N) *  (dihedral(N, vec3(0,0,1))) ;
}

mat3 generate_metric(vec3 _x){
	vec3 x = normalize(_x);
	
	//hack to avoid division by zero when extracting the eigenvectors :)
	//without this i think there'd need to be a bunch of if statements for all the different cases
	//like if x = (.5,.5,0) or x = (0,.5,.5) or x = (0,0,1), etc.
	x = max(abs(x),1e-5) * (mix(vec3(1.0f), vec3(-1.0f), greaterThan(x, vec3(0.0))));
	
	vec3 evals = vec3(0,0,1.0);
	vec3 evec0, evec1, evec2;
	//eigenvectors of outerproduct(x, x) have a closed form solution:

	evec0 = normalize(vec3(-x.z / x.x,0,1));
	evec1 = normalize(vec3(-x.y / x.x,1.0,0.0));
	evec2 = normalize(vec3(x.x/x.z, x.y/x.z, 1.0));    
	
	float k = .51;
	int dimensions = 3;
	float denom = 1. / float(dimensions - 1);
	vec3 mapped_evals = fitrange_3(evals, 0.0, 1.0,  (EIGENVALUE_SUM - k) * denom - 1e-4, k);
	//vec3 mapped_evals = vec3(1.5, 1.5, 1.);
	return  mapped_evals.x * outerprod3(evec0, evec0) +
			mapped_evals.y * outerprod3(evec1, evec1) + 
			mapped_evals.z * outerprod3(evec2, evec2);
}
mat2  generate_metric_2d(vec2 _x){
	vec2 x = normalize(_x);
	
	//hack to avoid division by zero when extracting the eigenvectors :)
	//without this i think there'd need to be a bunch of if statements for all the different cases
	//like if x = (.5,.5,0) or x = (0,.5,.5) or x = (0,0,1), etc.
	x = max(abs(x),1e-5) * (mix(vec2(1.0f), vec2(-1.0f), greaterThan(x, vec2(0.0))));
	
	vec2 evals = vec2(0,1.0);
	vec2 evec0, evec1;
	//eigenvectors of outerproduct(x, x) have a closed form solution:

	evec0 = normalize(vec2(-x.y / x.x,1));
	evec1 = normalize(vec2(x.x / x.y,1));
	
	float k = .51;
	int dimensions = 2;
	float denom = 1. / float(dimensions - 1);
	vec2 mapped_evals = fitrange_2(evals, 0.0, 1.0,  (EIGENVALUE_SUM - k) * denom - 1e-4, k);
	//vec3 mapped_evals = vec3(1.5, 1.5, 1.);
	return  mapped_evals.x * outerprod2(evec0, evec0) +
			mapped_evals.y * outerprod2(evec1, evec1); 
}

float steerable_perlin(vec3 pos, mat3  metric){
	vec3 noise_p = floor(pos);
	vec3 noise_f = pos - noise_p;//frac(pos);
	
	float out_val = 0.0;
	
	//perlin weights
	vec3 blend = interp3(noise_f);
	
	//we opt to use a weighted average instead of lerps.
	//if we remove the anisotropy from the dot product and the interpolation
	//the weighted average results in the same output as standard perlin noise.
	for(int i = 0; i <= 1; i++)
	for(int j = 0; j <= 1; j++)
	for(int k = 0; k <= 1; k++){

		vec3 o = vec3(float(i),float(j),float(k));
		
		vec3 g = noise_p  + o;
		vec3 r = rsphere(g);
		vec3 v = (o - noise_f);
		
		vec3 metric_v = metric * v;
		
		//regular perlin is dot(r, v)
		//applying the metric to v, adds one level of anisotropy
		float d = dot(r, metric_v);
		
		vec3 wv = abs(o-blend);
		
		//we get another level of anisotropy by introducing anisotropic weights to the algorithm       
		float w =  wv.x * wv.y * wv.z * interp(dot(v, metric_v));
		out_val += d * w;

	} 

	return out_val;
}


float steerable_perlin_projected(vec3 pos, mat2  metric, mat3 projection){
	vec3 noise_p = floor(pos);
	vec3 noise_f = pos - noise_p;//frac(pos);
	
	float out_val = 0.0;
	
	//perlin weights
	vec3 blend = interp3(noise_f);
	

	//we opt to use a weighted average instead of lerps.
	//if we remove the anisotropy from the dot product and the interpolation
	//the weighted average results in the same output as standard perlin noise.
	for(int i = 0; i <= 1; i++)
	for(int j = 0; j <= 1; j++)
	for(int k = 0; k <= 1; k++){

		vec3 o = vec3(float(i),float(j),float(k));
		
		vec3 g = noise_p  + o;
		vec2 r = ( rsphere(g) * projection).xy;
		vec2 v = ((o - noise_f) * projection).xy;
		
		vec2 metric_v = metric * v;
		
		//regular perlin is dot(r, v)
		//applying the metric to v, adds one level of anisotropy
		float d = dot(r, metric_v);
		
		vec3 wv = abs(o-blend);
		
		//we get another level of anisotropy by introducing anisotropic weights to the algorithm       
		float w =  wv.x * wv.y * wv.z * interp(dot(v, metric_v));

		out_val += d * w;

	} 

	return out_val;
}


float fbm(vec3 p, mat3 metric)
{

	float out_val  = 0.0;

	for(int i = 0; i < NUMBER_OF_OCTAVES; i++){
		out_val += pow(OCTAVE_BIAS, float(i)) * steerable_perlin(pow(2.0, float(i)) * (p + NOISE_OFFSET) * NOISE_FREQUENCY , metric);
	}
			
	return out_val;
}

//this is really only useful if you want a mesh to move through the noise field.
float fbm_artifact_free(vec3 p, mat3 metric)
{


	float out_val  = 0.0;

	for(int i = 0; i < NUMBER_OF_OCTAVES; i++){
		//since the weights are always heighest at the .5 position, combine two noises at the same octave to remove artifacts.
		out_val += pow(OCTAVE_BIAS, float(i)) * (steerable_perlin(pow(2.0, float(i)) * (p + NOISE_OFFSET) * NOISE_FREQUENCY , metric) + steerable_perlin(pow(2.0, float(i)) * ((p + NOISE_OFFSET + vec3(.5,.5,.5)) * NOISE_FREQUENCY ), metric)) * .5;
	}
			
	return out_val;
}
float fbm_projected(vec3 p, mat2 metric, mat3 projection)
{

	float out_val  = 0.0;

	for(int i = 0; i < NUMBER_OF_OCTAVES; i++){
		out_val += pow(OCTAVE_BIAS, float(i)) * steerable_perlin_projected(pow(2.0, float(i)) * (p + NOISE_OFFSET) * NOISE_FREQUENCY , metric, projection);
	}
			
	return out_val;
}



void fragment() {
	float c;
	vec3 anisotropy_dir = world_position.zyx * vec3(1,0,-1);
	if(Full_3D_Noise >= 1){
		mat3 x = generate_metric(anisotropy_dir);
		c = fbm_artifact_free(world_position, x);
	}else{
		mat3 projection = make_projection(normal);
		mat2 x = generate_metric_2d((anisotropy_dir * projection).xy);
		c = fbm_projected(world_position, x, projection);		
	}
	ALBEDO = vec3(c);
}

