#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"
#include "sphere.h"
#include <time.h>
#include "inline_func.cpp"

extern bool disable_hierarchy;

Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3), anti_aliasing_samples(1)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find and return the Hit structure for the closest intersection.  Be careful
// to ensure that hit.dist>=small_t.
Hit Render_World::Closest_Intersection(const Ray& ray)
{
    Hit hit = {NULL, 0, 0}, temp;
    std::vector<int> candidates;
    //s1 = clock();
    hierarchy.Intersection_Candidates(ray, candidates);
    //e1 = clock();
    //if(debug_pixel) std::cout<<"entries.size: "<<hierarchy.entries.size()<<std::endl;
    //if(debug_pixel) std::cout<<"candidates.size: "<<candidates.size()<<std::endl;
    for(auto i: candidates){
        temp = hierarchy.entries[i].obj->Intersection(ray, hierarchy.entries[i].part);
        if(temp.dist>=small_t && (!hit.object || temp.dist<hit.dist)){
            hit = temp;
        }
    }
    for(auto en: hierarchy.entries_inf){
        temp = en.obj->Intersection(ray, en.part);
        if(temp.dist>=small_t && (!hit.object || temp.dist<hit.dist)){
            hit = temp;
        }
    }
    //if(debug_pixel) std::cout<<"dist: "<<hit.dist<<" part: "<<hit.part<<std::endl;
    return hit;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    // set up the initial view ray here
    std::vector<Ray> ray;
    vec3 color;
    camera.Set_Ray(pixel_index, ray, anti_aliasing_samples);
    for(auto r: ray){
        color += Cast_Ray(r, 1);
    }
    color /= ray.size();
    //Ray ray = camera.Set_Ray(pixel_index);
    //vec3 color = Cast_Ray(ray,1);
    camera.Set_Pixel(pixel_index,Pixel_Color(color));
}

void Render_World::Render()
{
    double start, end;
    if(!disable_hierarchy)
        Initialize_Hierarchy();

    start = clock();
    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));

    end = clock();
    std::cout<<"render time: "<<(end-start)/CLOCKS_PER_SEC<<std::endl;
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
    // determine the color here
    vec3 intersection_point, normal;
    Hit hit = Closest_Intersection(ray);
    if(!hit.object){
        return background_shader->Shade_Surface(ray, \
            intersection_point, normal, recursion_depth);
    }
    else{
        intersection_point = ray.Point(hit.dist);
        normal = hit.object->Normal(intersection_point, hit.part);
        //if(debug_pixel) std::cout<<normal[0]<<std::endl;
        return hit.object->material_shader->Shade_Surface(ray, \
            intersection_point, normal, recursion_depth);
    }
}

void Render_World::Initialize_Hierarchy()
{
    for(auto obj: objects){
        for(int i=0; i<obj->number_parts; i++){
            if(std::isinf(obj->Bounding_Box(i).hi[0])){
                hierarchy.entries_inf.push_back({obj, i, obj->Bounding_Box(i)});
            }
            else hierarchy.entries.push_back({obj, i, obj->Bounding_Box(i)});
        }
    }
    // Fill in hierarchy.entries; there should be one entry for
    // each part of each object.
    hierarchy.Reorder_Entries();
    hierarchy.Build_Tree();
    //std::cout<<hierarchy.entries.size()<<" "<<hierarchy.tree.size()<<std::endl;
}

/*for(auto en: hierarchy.entries){
        temp = en.obj->Intersection(ray, en.part);
        if(temp.dist>=small_t && (!hit.object || temp.dist<hit.dist)){
            hit = temp;
        }
    }
*/