#include "globe.h"
#include "ray.h"

// Determine if the ray intersects with the globe
Hit Globe::Intersection(const Ray& ray, int part) const
{
    vec3 w = ray.endpoint - center;
    // t = t1 - sqrt(t2);
    // w * u > 0?
    double t1 = - dot(w, ray.direction), t2;
    t2 = t1 * t1 - dot(w, w) + radius * radius;
    if(t2<0) return {NULL, 0, 0};
    else{
        if(t1-sqrt(t2)>small_t)
            return {this, t1 - sqrt(t2), part};
        else if(t1+sqrt(t2)>small_t)
            return {this, t1 + sqrt(t2), part};
        else return {NULL, 0, 0};
    }
}

vec3 Globe::Normal(const vec3& point, int part) const
{
    vec3 normal;
    // compute the normal direction
    normal = (point - center).normalized();
    return normal;
}

Box Globe::Bounding_Box(int part) const
{
    Box box;
    box.lo = center;
    box.hi = center;
    for(int i=0; i<3; i++){
        box.lo[i] -= radius;
        box.hi[i] += radius;
    }
    return box;
}
