#include <limits>
#include "box.h"
#include "object.h"


// Compute the smallest box that contains both *this and bb.
Box Box::Union(const Box& bb) const
{
    Box box;
    box.lo = componentwise_min(bb.lo, lo);
    box.hi = componentwise_max(bb.hi, hi);
    return box;
}

// Enlarge this box (if necessary) so that pt also lies inside it.
void Box::Include_Point(const vec3& pt)
{
    lo = componentwise_min(pt, lo);
    hi = componentwise_max(pt, hi);
}

// Create a box to which points can be correctly added using Include_Point.
void Box::Make_Empty()
{
    lo.fill(std::numeric_limits<double>::infinity());
    hi=-lo;
}

/*
bool Box::Intersection(const Ray& ray) const
{
    // t = (lo-ray.e)*n/(ray.d*n), normal = (0,0,1) or (0,1,0) or (0,0,1)
    double t, dn;
    vec3 x;
    for(int i=0; i<3; i++){
        dn = ray.direction[i];
        if(std::abs(dn)>1e-4){
            t = (lo-ray.endpoint)[i]/dn;
            if(t>=small_t){
                x = ray.Point(t);
                if((x[0]<=hi[0]+small_t && x[0]>=lo[0]-small_t) &&\
                    (x[1]<=hi[1]+small_t && x[1]>=lo[1]-small_t) &&\
                    (x[2]<=hi[2]+small_t && x[2]>=lo[2]-small_t))
                    return true;
            }
        }
    }
    for(int i=0; i<3; i++){
        dn = ray.direction[i];
        if(std::abs(dn)>1e-4){
            t = (hi-ray.endpoint)[i]/dn;
            if(t>=small_t){
                x = ray.Point(t);
                if(t>small_t && (x[0]<=hi[0]+small_t && x[0]>=lo[0]-small_t) &&\
                    (x[1]<=hi[1]+small_t && x[1]>=lo[1]-small_t) &&\
                    (x[2]<=hi[2]+small_t && x[2]>=lo[2]-small_t))
                    return true;
            }
        }
    }
    return false;
}*/