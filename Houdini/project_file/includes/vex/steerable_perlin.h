matrix3 vec_projector(vector N){
    return matrix3(1) - outerproduct(N,N);
}

matrix3 make_projection(vector N){
    return vec_projector(N) * matrix3(dihedral(N, set(0,0,1)));
}

vector rsphere(vector p) {
        vector vals = random(p);
        
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
        return set(x,y,z);
}

float interp(float u){
        return 1.0 - smooth(0.0, 1.0, abs(u));
}

vector interp3(vector u){
        return set(interp(u.x), interp(u.y), interp(u.z));
}


matrix3 generate_metric(vector _x; float eigenvalue_sum, k){
        vector x = normalize(_x);
        
        //hack to avoid division by zero when extracting the eigenvectors :)
        //without this i think there'd need to be a bunch of if statements for all the different cases
        //like if x = (.5,.5,0) or x = (0,.5,.5) or x = (0,0,1), etc.
        x = max(abs(x),1e-5);
        x.x *= _x.x > 0 ? 1 : -1;
        x.y *= _x.y > 0 ? 1 : -1;
        x.z *= _x.z > 0 ? 1 : -1;

        
        vector evals = set(0,0,1.0);
        vector evec0, evec1, evec2;
        //eigenvectors of outerproduct(x, x) have a closed form solution:

        evec0 = normalize(set(-x.z / x.x,0,1));
        evec1 = normalize(set(-x.y / x.x,1.0,0.0));
        evec2 = normalize(set(x.x/x.z, x.y/x.z, 1.0));    
        
        int dimensions = 3;
        float denom = 1. / float(dimensions - 1);
        vector mapped_evals = fit(evals, 0.0, 1.0,  (eigenvalue_sum - k) * denom - 1e-4, k);
        //vector mapped_evals = vector(1.5, 1.5, 1.);
        return  mapped_evals.x * outerproduct(evec0, evec0) +
                        mapped_evals.y * outerproduct(evec1, evec1) + 
                        mapped_evals.z * outerproduct(evec2, evec2);
}


matrix2 generate_metric(vector2 _x; float eigenvalue_sum, k){
        vector2 x = normalize(_x);
        
        //hack to avoid division by zero when extracting the eigenvectors :)
        //without this i think there'd need to be a bunch of if statements for all the different cases
        //like if x = (.5,.5,0) or x = (0,.5,.5) or x = (0,0,1), etc.
        x = max(abs(x),1e-5) * lerp(vector2(1.0f), vector2(-1.0f), set(x.x > 0, x.y > 0));
        
        vector2 evals = set(0,1.0);
        vector2 evec0, evec1;
        //eigenvectors of outerproduct(x, x) have a closed form solution:

        evec0 = normalize(set(-x.y / x.x,1.0));
        evec1 = normalize(set(x.x / x.y,1.0));  
        
        int dimensions = 2;
        float denom = 1. / float(dimensions - 1);
        vector2 mapped_evals = fit(evals, 0.0, 1.0,  (eigenvalue_sum - k) * denom - 1e-4, k);

        return  mapped_evals.x * outerproduct(evec0, evec0) +
                mapped_evals.y * outerproduct(evec1, evec1);
}


float steerable_perlin(vector pos; matrix3  metric){
        vector noise_p = floor(pos);
        vector noise_f = pos - noise_p;//frac(pos);
        
        float out_val = 0.0;
        
        //perlin weights
        vector blend = interp3(noise_f);
        
        float w_sum = 0.0;
        //we opt to use a weighted average instead of lerps.
        //if we remove the anisotropy from the dot product and the interpolation
        //the weighted average results in the same output as standard perlin noise.
        for(int i = 0; i <= 1; i++)
        for(int j = 0; j <= 1; j++)
        for(int k = 0; k <= 1; k++){

                vector o = set(float(i),float(j),float(k));
                
                vector g = noise_p  + o;
                vector r = rsphere(g);
                vector v = (o - noise_f);
                
                vector metric_v = metric * v;
                
                //regular perlin is dot(r, v)
                //applying the metric to v, adds one level of anisotropy
                float d = dot(r, metric_v);
                
                vector wv = abs(o-blend);
                
                //we get another level of anisotropy by introducing anisotropic weights to the algorithm       
                float w =  wv.x * wv.y * wv.z * interp(dot(v, metric_v));
                w_sum += w;
                out_val += d * w;

        } 

        return out_val;// / w_sum;
}

float steerable_perlin(vector pos; matrix2 metric; matrix3 projection){
        vector noise_p = floor(pos);
        vector noise_f = pos - noise_p;//frac(pos);
        
        float out_val = 0.0;
        
        //perlin weights
        vector blend = interp3(noise_f);
        
        //we opt to use a weighted average instead of lerps.
        //if we remove the anisotropy from the dot product and the interpolation
        //the weighted average results in the same output as standard perlin noise.
        for(int i = 0; i <= 1; i++)
        for(int j = 0; j <= 1; j++)
        for(int k = 0; k <= 1; k++){

                vector o = set(float(i),float(j),float(k));
                
                vector g = noise_p  + o;
                vector r = rsphere(g) * projection;
                vector2 r2d = set(r.x,r.y);
                vector v = (o - noise_f) * projection;
                vector2 v2d = set(v.x, v.y);
                vector2 metric_v = metric * v2d;
                
                //regular perlin is dot(r, v)
                //applying the metric to v, adds one level of anisotropy
                float d = dot(r2d, metric_v);
                
                vector wv = abs(o-blend);
                
                //we get another level of anisotropy by introducing anisotropic weights to the algorithm       
                float w =  wv.x * wv.y * wv.z * interp(dot(v2d, metric_v));
                out_val += d * w;

        } 

        return out_val;
}

/*  

example useage:

matrix3 m = generate_metric(curlnoise(@P), 4.0, .5);
noise = fbm_3d(@P, chv("frequency"), chv("offset"), ch("bias"), chi("octaves"), m);

*/

float steerable_perlin_fbm_3d(vector p, freq, offset; float bias; int _octaves; matrix3 metric)
{

    int octaves = max(1,_octaves);
    float out_val  = 0;
    
    for(int i = 0; i < octaves; i++){
        out_val += pow(bias, float(i)) * steerable_perlin(pow(2.0, float(i)) * (p + offset) * freq , metric);
    }
            
    return out_val;
}

/*  

example useage:
matrix3 projection =make_projection(v@N);

vector aniso_dir = curlnoise(@P);
aniso_dir = aniso_dir * projection;

vector2 aniso_dir_projected = set(aniso_dir.x, aniso_dir.y)
matrix2 m = generate_metric(aniso_dir_projected, 4.0, .5);
noise = fbm_projected(@P, chv("frequency"), chv("offset"), ch("bias"), chi("octaves"), m, projection);

*/


float steerable_perlin_fbm_projected(vector p, freq, offset; float bias; int _octaves; matrix2 metric; matrix3 projection){
    int octaves = max(1,_octaves);
    float out_val  = 0;
    
    for(int i = 0; i < octaves; i++){
        out_val += pow(bias, float(i)) * steerable_perlin(pow(2.0, float(i)) * (p + offset) * freq , metric, projection);
    }
            
    return out_val;
}