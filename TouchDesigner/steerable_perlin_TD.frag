/* License:
Copyright 2025 Jacob Rice

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
layout(location = 0) out vec4 fragColor;

#define tolerance 0.1e-20
#define PI 3.1415926535

uniform float iTime;

uniform vec3 iResolution;

uniform int ANISOTROPY_MAP;//0 is an irrotational field, 1 is the gradient of the image in channel 0, scaled up so that it's not super noisey....

uniform float ANISOTROPY_VECTOR_SCALE; //multiplier on the anisotropy vector, to see how vector length affects anisotropy

uniform int ANIMATED_NOISE;

uniform float ANISOTROPY_STRENGTH; //setting this to zero results in effectively perlin noise (it's technically different though)

uniform int CROSS_FIELD;

uniform int DRAW_VECTOR_FIELD_VISUALIZER;

uniform float EIGENVALUE_SUM;//a value of 4.0 means the noise can be normalized, any value higher results in more anisotropy, but allows for zero valued weights.

uniform int MAKE_NOISE_POSITIVE_ONLY;

uniform float NOISE_FREQUENCY;

uniform int NOISE_ORDER; //how many extra neighbors to look up.

uniform int NUMBER_OF_OCTAVES;

uniform float OCTAVE_BIAS;

uniform int ROTATE_FEATURES;

float square(float a){
    return a*a;
}
float metric_gamma(float x, float n){
    return pow(n + abs(x),-1.0);
}

float fitrange(float x, float in_low, float in_high, float out_low, float out_high){
    float u = clamp(x, in_low, in_high);
    return ((out_high - out_low) * (u - in_low)) / (in_high - in_low) + out_low;
}

vec2 fitrange_2(vec2 x, float in_low, float in_high, float out_low, float out_high){
    return vec2(fitrange(x.x, in_low, in_high, out_low, out_high),
                  fitrange(x.y, in_low, in_high, out_low, out_high));
}

void solve_eig(mat2 M, inout vec2 evals, inout vec2 evec0, inout vec2 evec1) 
{

    //references:
    //https://github.com/AndreaCensi/2x2_matrix_eigenvalues
    //Andrea Censi's 2x2 matrix eigenvalues code
    
    float A = M[0][0], B = M[0][1], C = M[1][0], D = M[1][1];
    if(B * C <= tolerance  ) {
        evals.x = A; evec0 = vec2(1.,0.);
        evals.y = D; evec1 = vec2(0.,1.);
        return;
    }

    float tr = A + D;
    float det = A * D - B * C;
    float S = sqrt( square(tr/2.0) - det );
    evals.x = tr/2.0+ S;
    evals.y = tr/2.0 - S;

    float SS = sqrt( max(square((A-D)/2.0) + B * C, 0.0) );
    if( A - D < 0.0 ) {
        evec0.x = C;
        evec0.y = - (A-D)/2.0 + SS;
        evec1.x = + (A-D)/2.0 - SS;
        evec1.y = B;
    } else {
        evec1.x = C;
        evec1.y = - (A-D)/2.0 - SS;
        evec0.x = + (A-D)/2.0 + SS;
        evec0.y = B;
    }
    evec0 /= sqrt(dot(evec0, evec0));
    evec1 /= sqrt(dot(evec1, evec1));
}


mat2 outerprod(vec2 x, vec2 y){
    return mat2(x.x * y.x, x.y * y.x, 
                x.x * y.y, x.y * y.y);
}

float random (vec2 st) {
    return fract(sin(dot(st.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}

vec2 rand_dir(vec2 p){
    float r = random(p);
    r *= PI * 2.0;
    if (ANIMATED_NOISE == 0) {
        return vec2(cos(r), sin(r));   
    } else {
        return vec2(cos(r + iTime), sin(r + iTime));
    }
}

float interp(float u){
    return 1.0 - smoothstep(0.0, 1.0, abs(u));
}

vec2 interp2(vec2 u){
    return vec2(interp(u.x), interp(u.y));
}

mat2 generate_metric(vec2 aniso_direction, int quality){

    mat2 m = outerprod(aniso_direction, aniso_direction);

    vec2 evals, evec0, evec1;
    solve_eig(m, evals, evec0, evec1);
    float k = .5 / (float(quality) + 1.0);
    int dimensions = 2;
    float denom = 1. / float(dimensions - 1);
    vec2 mapped_evals = fitrange_2(evals, 0.0, dot(aniso_direction, aniso_direction),  (EIGENVALUE_SUM - k) * denom - 1e-4, k);
    
    return  mapped_evals.x * outerprod(evec0, evec0)  + mapped_evals.y * outerprod(evec1, evec1);

}

float aniso_perlin(vec2 p, mat2 metric){
    vec2 noise_p = floor(p);
    vec2 noise_f = fract(p);
    float out_val = 0.0;
    int start = -(NOISE_ORDER);
    int end = NOISE_ORDER + 1;
    float scale = 2. / float(abs(start) + end + 1);
    
    for(int i = start; i <= end; i++)
    for(int j = start; j <= end; j++){
        vec2 o = vec2(i,j);
        vec2 g = noise_p  + o;
        vec2 r = rand_dir(g); // random vector
        vec2 v = o - noise_f; //dir to corner
        vec2 metric_v = v * metric;
        float d = dot(r, metric_v); //inner product
        float w = interp(v.x * scale) * interp(v.y * scale);
        w *= interp(dot(v, metric_v)); //aniso weights
        out_val += d * w;
    }        
    return out_val;
    
}

vec2 image_grad(sampler2D img, vec2 p, float h){
    float grad_x = length(texture(img, p - vec2(h, 0))) - length(texture(img, p + vec2(h, 0)));
    float grad_y = length(texture(img, p - vec2(0, h))) - length(texture(img, p + vec2(0, h)));
    return vec2(grad_x, grad_y) / (2. * h);
}

void main()
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = gl_FragCoord.xy/iResolution.xy;
    vec2 p = uv;
    
    p -= .5;
    p *= 2.;
    p.x *= iResolution.x / iResolution.y;
    vec2 aniso_dir;
    float theta = PI / 2.0;
    
    if (ANISOTROPY_MAP == 0) {
        aniso_dir = vec2(p.y, -p.x);
    } else {
        float h = .05;
        vec2 sample_pos = uv * .2;
        aniso_dir = image_grad(sTD2DInputs[0], sample_pos, h);
    }
    
    if(ROTATE_FEATURES == 1) {
        aniso_dir = vec2(aniso_dir.y, -aniso_dir.x);
    }

    mat2 metric = generate_metric(aniso_dir * ANISOTROPY_VECTOR_SCALE, NOISE_ORDER);
    
    metric = ANISOTROPY_STRENGTH * metric + mat2(1) * (1.0 - ANISOTROPY_STRENGTH);

    float bias = OCTAVE_BIAS;
    float freq = NOISE_FREQUENCY;
    int octaves = NUMBER_OF_OCTAVES;

    float out_val = 0.0;

    for(int i = 0; i < octaves; i++){
        out_val += pow(bias, float(i)) * aniso_perlin(pow(2.0, float(i)) * p * freq, metric);
    }

    if( CROSS_FIELD == 1) {
        float out_val2 = 0.0;
        mat2 r = mat2(cos(theta), -sin(theta), sin(theta),cos(theta));
        mat2 metric2 = generate_metric(aniso_dir * ANISOTROPY_VECTOR_SCALE * r, NOISE_ORDER);
        metric2 = ANISOTROPY_STRENGTH * metric2 + mat2(1) * (1.0 - ANISOTROPY_STRENGTH);
        for(int i = 0; i < octaves; i++){
            out_val2 += pow(bias, float(i)) * aniso_perlin(pow(2.0, float(i)) * p * freq, metric2);
        }
        out_val = max(out_val, out_val2);
    }
    
    if(MAKE_NOISE_POSITIVE_ONLY == 1) {
        out_val = (out_val + 1.0) * .5;
    }

    vec3 col = vec3(out_val);

    if (DRAW_VECTOR_FIELD_VISUALIZER == 1) { 
        vec2 vis_pos = fract(p * 10.) - .5;
        float line_dist = length(vis_pos - aniso_dir * clamp(dot(vis_pos, aniso_dir) / dot(aniso_dir, aniso_dir), 0., 1.));
        line_dist = 1. - smoothstep(.01, .04, line_dist);
        col = max(col, vec3(line_dist, line_dist, 0.));
    }

    fragColor = vec4(col * 1.0,1.0);
}