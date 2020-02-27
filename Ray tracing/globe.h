#ifndef __GLOBE_H__
#define __GLOBE_H__

#include "object.h"
class Globe : public Object
{
    vec3 center;
    double radius, refractive;

public:
    Globe(const vec3& center_input,double radius_input)
        :center(center_input),radius(radius_input),\
        refractive(0)
    {}

    virtual Hit Intersection(const Ray& ray, int part) const override;
    virtual vec3 Normal(const vec3& point, int part) const override;
    virtual Box Bounding_Box(int part) const override;
};
#endif
