#include "ray_tracing.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <limits>

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    float A, a1, a2, a3;
    float a, b, c;
    A = glm::length(glm::cross(v1 - v0, v2 - v0)) / 2.0f;
    a1 = glm::length(glm::cross(v2 - v1, p - v1)) / 2.0f;
    a2 = glm::length(glm::cross(v0 - v2, p - v2)) / 2.0f;
    a3 = glm::length(glm::cross(v1 - v0, p - v0)) / 2.0f;

    a = a1 / A;
    b = a2 / A;
    c = a3 / A;
    if (a < 0.0f || b < 0.0f || c < 0.0f)
        return false;


    float epsilon = 10e-3f;
    if (((a + b + c) - 1.0f) < epsilon)
        return true;

    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    // no intersection
    if (glm::dot(ray.direction, plane.normal) == 0.0f)
        return false;

    float temp_t = (plane.D - glm::dot(ray.origin, plane.normal)) / glm::dot(ray.direction, plane.normal); // find new ray.t

    if (ray.t < temp_t)
        return false;


    if (temp_t > 0.0f) {
        ray.t = temp_t;
        return true;
    }

    return false;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    plane.normal = glm::cross((v0 - v2), (v1 - v2));
    plane.normal = glm::normalize(plane.normal);

    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    Plane plane = trianglePlane(v0, v1, v2);
    float curr_t = ray.t;
    if (intersectRayWithPlane(plane, ray) == false) {

        return false;
    }

    glm::vec3 p = ray.origin + ray.t * ray.direction;
    if (pointInTriangle(v0, v1, v2, plane.normal, p)) {
        float bigTriangleArea = glm::dot((v1 - v0), (v2 - v0)) / 2.0f;

        //P = w*v0 + u*v1 + v*v2

        float w = (glm::length(glm::cross((p - v1), (v2 - v1))) / 2.0f) / bigTriangleArea;
        float u = (glm::length(glm::cross((p - v0), (v2 - v0))) / 2.0f) / bigTriangleArea;
        float v = (glm::length(glm::cross((v0 - v1), (p - v1))) / 2.0f) / bigTriangleArea;

        glm::vec3 normalV0 = glm::cross(v2 - v0, v1 - v0);
        glm::vec3 normalV1 = glm::cross(v0 - v1, v2 - v1);
        glm::vec3 normalV2 = glm::cross(v0 - v2, v1 - v2);

        hitInfo.normal = glm::normalize(plane.normal);
        return true;
    }

    ray.t = curr_t;
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    float A = pow(ray.direction.x, 2.0) + pow(ray.direction.y, 2.0) + pow(ray.direction.z, 2.0);
    float B = 2.0f * (ray.direction.x * (ray.origin.x - sphere.center.x) + ray.direction.y * (ray.origin.y - sphere.center.y) + ray.direction.z * (ray.origin.z - sphere.center.z));
    float C = pow(ray.origin.x - sphere.center.x, 2.0) + pow(ray.origin.y - sphere.center.y, 2.0) + pow(ray.origin.z - sphere.center.z, 2.0) - pow(sphere.radius, 2.0);
    float t_0, t_1;



    float discriminant = pow(B, 2.0) - (4.0f * A * C);

    if (discriminant < 0.0f) {
        return false;
    }


    if (discriminant == 0.0f) {
        t_0 = (-B) / (2.0f * A);
        if (t_0 > 0.0f) {

            if (ray.t < t_0)
                return false;

            ray.t = t_0;
            hitInfo.normal = glm::normalize((ray.origin + ray.t * ray.direction) - sphere.center);
            return true;
        }
        return false;
    }

    t_0 = ((-B) + sqrt(discriminant)) / (2.0f * A);
    t_1 = ((-B) - sqrt(discriminant)) / (2.0f * A);


    if (t_0 > 0.0f && t_0 < t_1) {

        if (ray.t < t_0)
            return false;

        ray.t = t_0;
        hitInfo.normal = glm::normalize((ray.origin + ray.t * ray.direction) - sphere.center);
        return true;
    }

    if (t_1 > 0.0f && t_1 < t_0) {

        if (ray.t < t_1)
            return false;

        ray.t = t_1;
        hitInfo.normal = glm::normalize((ray.origin + ray.t * ray.direction) - sphere.center);
        return true;
    }

    // ray from the inside
    if (t_1 > 0.0f && t_0 < 0.0f) {
        if (ray.t < t_1)
            return false;

        ray.t = t_1;
        hitInfo.normal = glm::normalize((ray.origin + ray.t * ray.direction) - sphere.center);
        return true;
    }

    if (t_0 > 0.0f && t_1 < 0.0f) {
        if (ray.t < t_0)
            return false;

        ray.t = t_0;
        hitInfo.normal = glm::normalize((ray.origin + ray.t * ray.direction) - sphere.center);
        return true;
    }

    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float tx_min = (box.lower.x - ray.origin.x) / ray.direction.x;
    float tx_max = (box.upper.x - ray.origin.x) / ray.direction.x;
    float ty_min = (box.lower.y - ray.origin.y) / ray.direction.y;
    float ty_max = (box.upper.y - ray.origin.y) / ray.direction.y;
    float tz_min = (box.lower.z - ray.origin.z) / ray.direction.z;
    float tz_max = (box.upper.z - ray.origin.z) / ray.direction.z;


    float t_inx, t_outx;
    float t_iny, t_outy;
    float t_inz, t_outz;
    float t_in, t_out;

    if (tx_min < tx_max) {
        t_inx = tx_min;
        t_outx = tx_max;
    }
    else {
        t_inx = tx_max;
        t_outx = tx_min;
    }

    if (ty_min < ty_max) {
        t_iny = ty_min;
        t_outy = ty_max;
    }
    else {
        t_iny = ty_max;
        t_outy = ty_min;
    }

    if (tz_min < tz_max) {
        t_inz = tz_min;
        t_outz = tz_max;
    }
    else {
        t_inz = tz_max;
        t_outz = tz_min;
    }



    if (t_inx > t_iny && t_inx > t_inz)
        t_in = t_inx;
    else if (t_iny > t_inx && t_iny > t_inz)
        t_in = t_iny;
    else
        t_in = t_inz;


    if (t_outx < t_outy && t_outx < t_outz)
        t_out = t_outx;
    else if (t_outy < t_outx && t_outy < t_outz)
        t_out = t_outy;
    else
        t_out = t_outz;


    if (t_in > t_out || t_out < 0)
        return false;



    if (t_in > 0.0f) {

        if (ray.t < t_in)
            return false;

        ray.t = t_in;
        return true;
    }

    if (t_in < 0.0f && t_out > 0.0f) {
        if (ray.t < t_out)
            return false;

        ray.t = t_out;
        return true;
    }

    return false;
}
