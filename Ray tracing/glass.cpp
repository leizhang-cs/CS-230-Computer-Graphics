#include "glass.h"
#include "ray.h"
#include <cfloat>
#include <limits>
#include "vec.h"

// Intersect with the half space defined by the glass.  The glass's normal
// points outside.
Hit Glass::Intersection(const Ray& ray, int part) const
{
    vec3 w = ray.endpoint - x1;
    double wn, t;
    // t = d * n. if t==0, ray is parallel with glass
    t = dot(ray.direction, normal);
    // wn = max(w*n, 0). if wn<0, ray is inside
    wn = dot(w, normal);
    if(abso(t)<small_t || abso(wn)<small_t) return {NULL,0,0};
    else{
        // (w+td) * n = 0, dist: t = -w*n/d*n
        t = - wn / t;
        return {this, t, part};
    }
}

vec3 Glass::Normal(const vec3& point, int part) const
{
    return normal;
}

// There is not a good answer for the bounding box of an infinite object.
// The safe thing to do is to return a box that contains everything.
Box Glass::Bounding_Box(int part) const
{
    Box b;
    b.hi.fill(std::numeric_limits<double>::infinity());
    b.lo=-b.hi;
    return b;
}
