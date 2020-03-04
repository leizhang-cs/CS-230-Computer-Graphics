#include "driver_state.h"
#include <cstring>
#include <vector>
#include <limits>


static float small_t = 1e-5;

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color = new pixel [width*height]();
    state.image_depth = new float [width*height];
    for(int i=0; i<width*height; i++) state.image_depth[i] = std::numeric_limits<float>::infinity();
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    std::cout<<"rendering..."<<std::endl;
    // calculate indices, three geometries per triangle.
    if(type==render_type::triangle && state.num_triangles==0 && state.num_vertices!=0){
        state.num_triangles = state.num_vertices/3;
        state.index_data = new int[state.num_vertices];
        for(int i=0; i<state.num_vertices; i++) state.index_data[i] = i;
    }
    else if(type==render_type::fan && state.num_triangles==0 && state.num_vertices!=0){
        state.num_triangles = state.num_vertices - 2;
        state.index_data = new int[state.num_triangles*3];
        for(int i=0, j=0; i<state.num_triangles; i++, j+=3){
            state.index_data[j] = 0;
            state.index_data[j+1] = i + 1;
            state.index_data[j+2] = i + 2;
        }
    }
    else if(type==render_type::strip && state.num_triangles==0 && state.num_vertices!=0){
        state.num_triangles = state.num_vertices - 2;
        state.index_data = new int[state.num_triangles*3];
        for(int i=0, j=0; i<state.num_triangles; i++, j+=3){
            state.index_data[j] = i;
            state.index_data[j+1] = i + 1;
            state.index_data[j+2] = i + 2;
        }
    }
    const data_geometry **triangle_geometry = new const data_geometry *[state.num_triangles];

    for(int i=0; i<state.num_triangles; i++){
        // three vertices per triangle geometry.
        data_geometry *out = new data_geometry[3];
        std::vector<data_vertex> in(3);

        for(int j=0, k=i*3; j<3; j++){
            out[j].data = new float[state.floats_per_vertex];
            int v_index = state.index_data[k+j];
            in[j].data = state.vertex_data + v_index*state.floats_per_vertex;
            state.vertex_shader(in[j], out[j], state.uniform_data);
        }
        triangle_geometry[i] = out;
    }
    // ???
    //delete [] state.index_data;

    std::cout<<"clipping..."<<std::endl;
    clip_triangle(state, triangle_geometry, 0);
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3], int face)
{
    if(face==6)
    {
        std::cout<<"rasterizing..."<<std::endl;
        rasterize_triangle(state, in);
        return;
    }

    // output triangles after clipping
    std::vector<data_geometry*> vec_triangle;
    // clip plane x y or z, xyz = face/2
    int xyz = face / 2;
    // plane: x, normal. x: left btm most or right top most
    vec3 x, normal;
    // Set clip plane of 3D world. So, convert vertex (divide w) into 3D world 
    // before judge
    switch (face)
    {
        case 0: xyz = 0; x = {-1,-1,-1}; normal = {1,0,0};  break;
        case 1: xyz = 0; x = {1,1,1};    normal = {-1,0,0}; break;
        case 2: xyz = 1; x = {-1,-1,-1}; normal = {0,1,0};  break;
        case 3: xyz = 1; x = {1,1,1};    normal = {0,-1,0}; break;
        case 4: xyz = 2; x = {-1,-1,-1}; normal = {0,0,1};  break;
        case 5: xyz = 2; x = {1,1,1};    normal = {0,0,-1}; break;
        default: break;
    }
    
    for(int i=0; i<state.num_triangles; i++){
        // vertices inside or outside judge: n*(v-x) > 0 ? inside: outside
        // equivalent to |d| < |w|. Specifically, clip 6 plane in turn. 
        // So, if w<0, d>w; if w>0, d<w.
        std::vector<int> index_in, index_out;
        
        for(int j=0; j<3; j++){
            
            float d = in[i][j].gl_Position[xyz]/in[i][j].gl_Position[3];
            float w = 1;
            if(face%2==0) w = -1;
            
            // save inside and outside vertex_index separately
            if((w<0 && d>w-small_t) || (w>0 && d<w+small_t)){
                index_in.push_back(j);
            }
            else{
                index_out.push_back(j);
            }
        }

        // triangle interpolate
        if(index_in.size()==0) continue;
        else if(index_in.size()==3){
            data_geometry *out = new data_geometry[3];
            for(int j=0; j<3; j++) out[j] = in[i][j];
            vec_triangle.push_back(out);
        }
        else{
            // interpolate coefficient. r: perspective incorrect, k: correct
            // k is for position. k or r is for other attributes
            float r1, r2, k1, k2;
            // k = v2 / (v1 - v2), P = k*v1 + (1-k)*v2, new vertices
            bool in_out = true;
            // two triangles in "in"
            if(index_in.size()==1){
                swap(index_in, index_out);
                in_out = false;
            }
            // n*(x-P) = 0 => nx = nP, P: intersect A: outside, BC: inside
            // k = n*(B-P)/n*(B-A) = (nB - nx) / (nB - nA), 0<k<1
            vec4 A = in[i][index_out[0]].gl_Position, B = in[i][index_in[0]].gl_Position,\
                C = in[i][index_in[1]].gl_Position;
            // convert value into 3D world (divide w)
            float nA = dot(normal, vec3(A))/in[i][index_out[0]].gl_Position[3];
            float nB = dot(normal, vec3(B))/in[i][index_in[0]].gl_Position[3];
            float nC = dot(normal, vec3(C))/in[i][index_in[1]].gl_Position[3];
            float nx = dot(normal, x);
            
            r1 = (nB - nx) / (nB - nA);
            r2 = (nC - nx) / (nC - nA);
            
            k1 = (r1*B[3])/((r1*B[3])+(1-r1)*A[3]);
            k2 = (r2*C[3])/((r2*C[3])+(1-r2)*A[3]);
            // invalid k
            if(k1<0 || k1>1 || k2<0 || k2>1) continue;

            // ck: color interpolate coefficient
            float ck1, ck2;
            data_geometry p1, p2;
            p1.data = new float[state.floats_per_vertex];
            p2.data = new float[state.floats_per_vertex];

            switch (state.interp_rules[state.floats_per_vertex-1])
            {
                case interp_type::flat: ck1=-1,ck2=-1; break;
                case interp_type::noperspective: ck1=r1,ck2=r2; break;
                case interp_type::smooth: ck1=k1, ck2=k2; break;
                default: assert("invalid interpolation type" && 0); break;
            }
            
            // position interpolate
            p1.gl_Position = (1-k1)*in[i][index_in[0]].gl_Position + \
                            k1*in[i][index_out[0]].gl_Position;
            p2.gl_Position = (1-k2)*in[i][index_in[1]].gl_Position + \
                            k2*in[i][index_out[0]].gl_Position;
            
            // attributes (color) interpolate
            if(ck1<0 && ck2<0){
                p1.data = in[i][0].data;
                p2.data = in[i][0].data;
            }
            else{
                for(int j=0; j<state.floats_per_vertex; j++){
                    p1.data[j] = (1-ck1)*in[i][index_in[0]].data[j] + \
                                ck1*in[i][index_out[0]].data[j];
                    p2.data[j] = (1-ck2)*in[i][index_in[1]].data[j] + \
                            ck2*in[i][index_out[0]].data[j];
                }
            }
            
            // add triangle. two inside or one inside
            if(in_out){
                data_geometry *out1 = new data_geometry[3], *out2 = new data_geometry[3];

                out1[0] = in[i][index_in[0]];
                out1[1] = in[i][index_in[1]];
                out1[2] = p1;
                out2[0] = in[i][index_in[1]];
                out2[1] = p1;
                out2[2] = p2;
                vec_triangle.push_back(out1);
                vec_triangle.push_back(out2);
            }
            else{
                data_geometry *out = new data_geometry[3];
                
                out[0] = in[i][index_out[0]];
                out[1] = p1;
                out[2] = p2;
                vec_triangle.push_back(out);
            }
        }
    }

    for(int i=0; i<state.num_triangles; i++) delete [] in[i];
    delete [] in;
    
    // transfer format
    const data_geometry** out_triangle = new const data_geometry *[vec_triangle.size()];
    state.num_triangles = vec_triangle.size();
    for(int i=0; i<state.num_triangles; i++)
        out_triangle[i] = vec_triangle[i];
    clip_triangle(state,out_triangle,face+1);
}

// optimized calculation for triangle area
template <class T >
T area(const vec<T,2> & u, const vec<T,2> & v)
{
    return u[0] * v[1] - u[1] * v[0];
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    std::cout<<"num_tri: "<<state.num_triangles<<std::endl;
    for(int i=0; i<state.num_triangles; i++){ // int &n = i;
        std::vector<vec4> v = {in[i][0].gl_Position, \
            in[i][1].gl_Position, in[i][2].gl_Position};
        
        std::vector<data_fragment> in_frag(3);
        for(int j=0; j<3; j++){
            in_frag[j].data = in[i][j].data;
        }
        
        // perspective correct, (x, y, z) divided by w
        for(int i=0; i<3; i++){
            v[i][0] /= v[i][3];
            v[i][1] /= v[i][3];
            v[i][2] /= v[i][3];
        }
        
        vec2 A = vec2(v[0]), B = vec2(v[1]), C = vec2(v[2]), \
            O = {0,0}, e1 = {1,0}, e2 = {0,1};
        // use barycentric coordinate to judge point inside or outside triangle
        // coefficient for (x,y): a = area(PBC)/area(ABC) = a0 + a1*x + a2*y
        float a, b, c, a0, a1, a2, b0, b1, b2;
        // optimize: a0 = area(OBC)/area(ABC), a0 + a1 = area(e1BC)/area(ABC), e1 = (1,0)
        float S_A, S_B, S_a0, S_a1, S_a2, S_b0, S_b1, S_b2;
        S_A = 1.0f / area(A-B, C-B); S_B = 1.0f / area(B-A, C-A);
        S_a0 = area(O-B, C-B); S_a1 = area(e1-B, C-B); S_a2 = area(e2-B, C-B);
        S_b0 = area(O-A, C-A); S_b1 = area(e1-A, C-A); S_b2 = area(e2-A, C-A);
        a0 = S_a0 * S_A; a1 = S_a1 * S_A - a0; a2 = S_a2 * S_A - a0;
        b0 = S_b0 * S_B; b1 = S_b1 * S_B - b0; b2 = S_b2 * S_B - b0;
        // perspective correct interpolate: k_pc = a/wa + b/wb + c/wc, w = 1/w
        float wa = 1/v[0][3], wb = 1/v[1][3], wc = 1/v[2][3], k_pc = 1;

        // border of triangle
        float xmin = std::min(std::min(A[0], B[0]), C[0]),\
            xmax = std::max(std::max(A[0], B[0]), C[0]),\
            ymin = std::min(std::min(A[1], B[1]), C[1]),\
            ymax = std::max(std::max(A[1], B[1]), C[1]);
        int left = (xmin+1)/2*state.image_width, right = (xmax+1)/2*state.image_width,\
            btm = (ymin+1)/2*state.image_height, top = (ymax+1)/2*state.image_height;
        
        // a = a + delta a, delta a = ka1 * delta x + ka2 * delta y; 
        float y = 2.0f/state.image_height * (btm-1) - 1 + 1.0f/state.image_height;
        float ay = a2*y, by = b2*y;
        float ka1 = a1*2.0f/state.image_width, ka2 = a2*2.0f/state.image_height, \
            kb1 = b1*2.0f/state.image_width, kb2 = b2*2.0f/state.image_height;
            
        for(int i=btm; i<=top; i++){
            ay += ka2;
            by += kb2;
            
            float x = 2.0f/state.image_width * (left-1) - 1 + 1.0f/state.image_width;
            // Initialize a b c again
            a = a0 + a1*x + ay; b = b0 + b1*x + by;
            for(int j=left; j<=right; j++){
                a += ka1;
                b += kb1;
                c = 1 - a - b;
                
                // Interpolate if (x,y) in triangle
                if(a>-small_t && b>-small_t && c>-small_t){
                    // z-buffer coefficient : k_zbuf, perspective-correct
                    vec3 k_interp, k_zbuf = {a,b,c};
                    data_output out;
                    float depth = 0;

                    // attributes (color) interpolate
                    switch (state.interp_rules[state.floats_per_vertex-1])
                    {
                        case interp_type::flat: k_interp = {1,0,0}; break;
                        case interp_type::noperspective: k_interp = {a,b,c}; break;
                        case interp_type::smooth: k_pc = 1 / (a*wa + b*wb + c*wc);
                                k_interp = {a*wa*k_pc,b*wb*k_pc,c*wc*k_pc}; break;
                        default: assert("invalid interpolation type" && 0); break;
                    }
                    
                    // z-buffer judge
                    for(int i=0; i<3; i++) depth += k_zbuf[i]*v[i][2];
                    if(depth<state.image_depth[i*state.image_width+j]){
                        state.image_depth[i*state.image_width+j] = depth;
                    }
                    else continue;
                    
                    // output image, call fragment shader
                    vec3 rgb;
                    for(int i=0; i<3; i++){
                        state.fragment_shader(in_frag[i], out, state.uniform_data);
                        for(int j=0; j<3; j++) rgb[j] += k_interp[i]*out.output_color[j]*255;
                    }
                    state.image_color[i*state.image_width+j] = make_pixel(rgb[0], rgb[1], rgb[2]);
                }
            }
        }
    }
    for(int i=0; i<state.num_triangles; i++) delete [] in[i];
    delete [] in;
}