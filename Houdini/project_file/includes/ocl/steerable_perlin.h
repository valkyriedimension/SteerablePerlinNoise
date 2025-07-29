#define PI 3.1415926535

float interp(float x){
    return 1.0f - smoothstep(0.0f, 1.0f, fabs(x));
}
float3 interp3(float3 x){
    return (float3)(interp(x.x), interp(x.y), interp(x.z));
}

float2 interp2(float2 x){
    return (float2)(interp(x.x), interp(x.y));
}

float random3(float3 pos){

    float _x;
    
    float x = sin((float)dot(pos, (float3)(64.25375463f, 23.27536534f, 86.29678483f))) * 59482.7542f;
    return x - floor(x);
}


float random2(float2 pos){

    float _x;
    
    float x = sin(dot(pos, (float2)(12.9898f,78.233f))) * 43758.5453123f;
    return (x - floor(x));
}
float3 random33(float3 pos){
    return (float3)(random3(pos + .01f), random3(pos + .02f), random3(pos + .03f));
}

float3 rsphere(float3 p) {
    float3 vals = random33(p);
    
    float u = vals.x;
    float v = vals.y;
    float theta = u * 2.0f * PI;
    float phi = acos(2.0f * v - 1.0f);
    float r = pow(vals.z, .33333f);
    float sinTheta = sin(theta);
    float cosTheta = cos(theta);
    float sinPhi = sin(phi);
    float cosPhi = cos(phi);
    float x = r * sinPhi * cosTheta;
    float y = r * sinPhi * sinTheta;
    float z = r * cosPhi;
    return (float3)(x,y,z);
}


float fitrange(float x, float in_low, float in_high, float out_low, float out_high){
    float u = clamp(x, in_low, in_high);
    return ((out_high - out_low) * (u - in_low)) / (in_high - in_low) + out_low;
}
float2 fitrange_2(float2 x, float in_low, float in_high, float out_low, float out_high){
    return (float2)(fitrange(x.x, in_low, in_high, out_low, out_high),
                fitrange(x.y, in_low, in_high, out_low, out_high));
}
float3 fitrange_3(float3 x, float in_low, float in_high, float out_low, float out_high){
    return (float3)(fitrange(x.x, in_low, in_high, out_low, out_high),
                fitrange(x.y, in_low, in_high, out_low, out_high),
                fitrange(x.z, in_low, in_high, out_low, out_high));
}


void vec_projector(float3 N, mat3 out){
    //I - outerprod(N,N);
    out[0] = (float3)(1,0,0);
    out[1] = (float3)(0,1,0);
    out[2] = (float3)(0,0,1);
    out[0] -= N * N.x;
    out[1] -= N * N.y;
    out[2] -= N * N.z;
}

void dihedral(float3 a, float3 b, mat3 out){
    out[0] = (float3)0.0f;
    out[1] = (float3)0.0f;
    out[2] = (float3)0.0f;

    float3 v = cross(a,b);
    float c = dot(a,b);
    mat3 I;
    I[0] = (float3)(1,0,0);
    I[1] = (float3)(0,1,0);
    I[2] = (float3)(0,0,1);

    mat3 v_prime;
    v_prime[0] = (float3)(0, -v.z, v.y);
    v_prime[1] = (float3)(v.z, 0, -v.x);
    v_prime[2] = (float3)(-v.y, v.x, 0);
    float scale = 1.0 / (1.f + c);
    mat3 v_prime_2;
    mat3mul(v_prime, v_prime, v_prime_2);
    out[0] = I[0] + v_prime[0] + v_prime_2[0] * scale;
    out[1] = I[1] + v_prime[1] + v_prime_2[1] * scale;    
    out[2] = I[2] + v_prime[2] + v_prime_2[2] * scale;    

}


void make_projection(float3 N, mat3 out){
    mat3 rot;
    mat3 projector;
    mat3 x;
    dihedral(N, (float3)(0,0,1), rot);
    vec_projector(N, projector);
    transpose3(rot, x);
    mat3mul(projector, x, out);
}


void generate_metric_3d(float3 aniso_dir, float eigenvalue_sum, float k, mat3 out){
    float3 x = normalize(aniso_dir);
    
    //hack to avoid division by zero when extracting the eigenvectors :)
    //without this i think there'd need to be a bunch of if statements for all the different cases
    //like if x = (.5,.5,0) or x = (0,.5,.5) or x = (0,0,1), etc.
    x = max(fabs(x),1e-5);
    
    x.x *= (aniso_dir.x >= 0.0f ? 1.0f : -1.0f);
    x.y *= (aniso_dir.y >= 0.0f ? 1.0f : -1.0f);   
    x.z *= (aniso_dir.z >= 0.0f ? 1.0f : -1.0f);   

    float3 evals =(float3)(0,0,1.0f);
    float3 evec0, evec1, evec2;
    //eigenvectors of outerproduct(x, x) have a closed form solution:

    evec0 = normalize((float3)(-x.z / x.x,0,1));
    evec1 = normalize((float3)(-x.y / x.x,1.0,0.0));
    evec2 = normalize((float3)(x.x/x.z, x.y/x.z, 1.0));    
    

    int dimensions = 3;
    float denom = 1. / (float)(dimensions - 1);
    float3 mapped_evals = fitrange_3(evals, 0.0, 1.0,  (eigenvalue_sum - k) * denom - 1e-4, k);

    mat3 e00;
    outerprod3(evec0 * mapped_evals.x, evec0, e00);
    mat3 e11;
    outerprod3(evec1 * mapped_evals.y, evec1, e11);
    mat3 e22;
    outerprod3(evec2 * mapped_evals.z, evec2, e22);    
    out[0] = e00[0] + e11[0] + e22[0];
    out[1] = e00[1] + e11[1] + e22[1];
    out[2] = e00[2] + e11[2] + e22[2];      
}



void generate_metric_2d(float2 aniso_dir,  float eigenvalue_sum, float k, mat2 *out){
    float2 x = normalize(aniso_dir);
    
    //hack to avoid division by zero when extracting the eigenvectors :)
    //without this i think there'd need to be a bunch of if statements for all the different cases
    //like if x = (.5,.5,0) or x = (0,.5,.5) or x = (0,0,1), etc.
    x = max(fabs(x),1e-5);
    x.x *= (aniso_dir.x >= 0.0f ? 1.0f : -1.0f);
    x.y *= (aniso_dir.y >= 0.0f ? 1.0f : -1.0f);  


    float2 evals =(float2)(0, 1.0);
    float2 evec0, evec1;
    //eigenvectors of outerproduct(x, x) have a closed form solution:

    evec0 = normalize((float2)(-x.y / x.x,1.0));
    evec1 = normalize((float2)(x.x / x.y,1.0));  
    
    int dimensions = 2;
    float denom = 1. / (float)(dimensions - 1);
    float2 mapped_evals = fitrange_2(evals, 0.0, 1.0,  (eigenvalue_sum - k) * denom - 1e-4, k);
    mat2 e00 = (float4)(evec0 * evec0.x, evec0 * evec0.y) * mapped_evals.x;
    mat2 e11 = (float4)(evec1 * evec1.x, evec1 * evec1.y) * mapped_evals.y;
    *out = 0;
    *out += e00;
    *out += e11;
}

float steerable_perlin_3d(float3 pos, mat3  metric){
    float3 noise_p = floor(pos);
    float3 noise_f = pos - noise_p;
    
    float out_val = 0.0;
    
    //perlin weights
    float3 blend = interp3(noise_f);
    
    //we opt to use a weighted average instead of lerps.
    //if we remove the anisotropy from the dot product and the interpolation
    //the weighted average results in the same output as standard perlin noise.
    for(int i = 0; i <= 1; i++)
    for(int j = 0; j <= 1; j++)
    for(int k = 0; k <= 1; k++){

        float3 o = (float3)(i,j,k);
        
        float3 g = noise_p  + o;
        float3 r = rsphere(g);
        float3 v = (o - noise_f);
        
        float3 metric_v = mat3vecmul(metric, v);
        
        //regular perlin is dot(r, v)
        //applying the metric to v, adds one level of anisotropy
        float d = dot(r, metric_v);
        
        float3 wv = fabs(o-blend);
        
        //we get another level of anisotropy by introducing anisotropic weights to the algorithm       
        float w =  wv.x * wv.y * wv.z * interp(dot(v, metric_v));

        out_val += d * w;

    } 

    return out_val;
}


float steerable_perlin_projected(float3 pos, mat2 metric, mat3 projection){
    float3 noise_p = floor(pos);
    float3 noise_f = pos - noise_p;
    
    float out_val = 0.0;
    
    //perlin weights
    float3 blend = interp3(noise_f);
    
    //we opt to use a weighted average instead of lerps.
    //if we remove the anisotropy from the dot product and the interpolation
    //the weighted average results in the same output as standard perlin noise.
    for(int i = 0; i <= 1; i++)
    for(int j = 0; j <= 1; j++)
    for(int k = 0; k <= 1; k++){

        float3 o = (float3)(i,j,k);
        
        float3 g = noise_p  + o;
        float2 r = mat3vecmul(projection, rsphere(g)).xy;
        float2 v = mat3vecmul(projection, (o - noise_f)).xy;
        
        float2 metric_v = mat2vecmul(metric, v);
        
        //regular perlin is dot(r, v)
        //applying the metric to v, adds one level of anisotropy
        float d = dot(r, metric_v);
        
        float3 wv = fabs(o-blend);
        
        //we get another level of anisotropy by introducing anisotropic weights to the algorithm       
        float w =  wv.x * wv.y * wv.z * interp(dot(v, metric_v));

        out_val += d * w;

    } 

    return out_val;
}

float steerable_perlin_2d(float2 pos, mat2 metric){
    float2 noise_p = floor(pos);
    float2 noise_f = pos - noise_p;
    
    float out_val = 0.0;
    
    //perlin weights
    float2 blend = interp2(noise_f);
    
    //we opt to use a weighted average instead of lerps.
    //if we remove the anisotropy from the dot product and the interpolation
    //the weighted average results in the same output as standard perlin noise.
    for(int i = 0; i <= 1; i++)
    for(int j = 0; j <= 1; j++){

        float2 o = (float2)(i,j);
        
        float2 g = noise_p + o;
        float2 r = (float2)(cos(random2(g) * 2.0f * PI), sin(random2(g) * 2.0f * PI));
        float2 v = (o - noise_f);
        
        float2 metric_v = mat2vecmul(metric, v);
        
        //regular perlin is dot(r, v)
        //applying the metric to v, adds one level of anisotropy
        float d = dot(r, metric_v);
        
        float2 wv = fabs(o-blend);
        
        //we get another level of anisotropy by introducing anisotropic weights to the algorithm       
        float w =  wv.x * wv.y *  interp(dot(v, metric_v));

        out_val += d * w;

    } 

    return out_val;
}