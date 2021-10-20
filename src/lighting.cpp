#include "lighting.h"
// Suppress warnings in third-party code.
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
#include <glm/gtc/constants.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <limits>

const float RAY_STEP = 1.0e-3f;     // Defines the 'step' to take in a ray's direction to prevent erroneous intersections. 

bool hardShadow(const Ray& ray, const PointLight& light, BoundingVolumeHierarchy &bvh)
{
    // Calculate the point of intersection and the shadow ray direction.
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 ray_direction = glm::normalize(light.position - p);

    // Construct the shadow ray.
    Ray shadowRay = { p + (RAY_STEP * ray_direction), ray_direction, glm::distance(p, light.position) };

    // Draw a red debug ray if the shadow ray hits another source.
    HitInfo hitInfo;
    if (bvh.intersect(shadowRay, hitInfo))
    {
        drawRay(shadowRay, glm::vec3{ 1.0f, 0.0f, 0.0f });
        return true;
    }

    // Draw a yellow debug ray if there is no intersection with another object
    drawRay(shadowRay, glm::vec3{ 1.0f, 1.0f, 0.0f });
    return false;
}

glm::vec3 lightRay(const Ray& ray, const HitInfo& hitInfo, const Scene& scene, BoundingVolumeHierarchy& bvh)
{
    // Calculate the point of intersection.
    glm::vec3 p = ray.origin + ray.t * ray.direction;

    // Draw a blue debug ray if the ray hit.
    drawRay(ray, glm::vec3{ 0.0f, 0.0f, 1.0f });

    //// Draw a green debug ray for the normal if the ray hit.
    //debugRay(p, hitInfo.normal, 1.0f, glm::vec3{ 0.0f, 1.0f, 0.0f });

    // Calculate the lighting considering each point light source.
    glm::vec3 lighting = glm::vec3{ 0.0f };
    for (const auto& light : scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            const PointLight pointLight = std::get<PointLight>(light);
            if (hardShadow(ray, std::get<PointLight>(light), bvh))
            {
                return lighting;
            }
        }
    }
    return glm::vec3{ 1.0f };
}
