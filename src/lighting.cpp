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

const int RECURSION_LIMIT = 3;
const float RAY_STEP = 1.0e-3f;     //The 'step' taken in a ray's direction to prevent wrong intersections. 
const float LINE_SAMPLE_LIMIT = 100;
const float LINE_PARALLELOGRAM_SAMPLE_LIMIT = 100;

static // Standard lambertian shading: Kd * dot(N,L), clamped to zero when negative. Where L is the light vector.
glm::vec3 diffuseOnly(HitInfo hitInfo, glm::vec3& vertexPos, glm::vec3 lightPos)
{
    glm::vec3 normal = hitInfo.normal;
    glm::vec3 lightVec = lightPos - vertexPos;
    glm::vec3 lambertian = hitInfo.material.kd * glm::dot(normal, glm::normalize(lightVec));

    if (glm::dot(normal, glm::normalize(lightVec)) < -10e-3f) return glm::vec3(0.0f);

    return lambertian;
}


glm::vec3 phongSpecularOnly(HitInfo hitInfo, glm::vec3& vertexPos, glm::vec3 lightPos, glm::vec3& cameraPos)
{
    glm::vec3 normal = hitInfo.normal; 
    glm::vec3 viewVector = cameraPos - vertexPos;
    viewVector = glm::normalize(viewVector);
    glm::vec3 lightVec = lightPos - vertexPos;
    lightVec = glm::normalize(lightVec);
    glm::vec3 reflectionV = 2 * glm::dot(lightVec, normal) * normal - lightVec;
    if (glm::dot(lightVec, normal) < -10e-3f) {
        return glm::vec3(0, 0, 0);
    }
    else {
        return glm::vec3{ hitInfo.material.ks * glm::pow(glm::dot(viewVector, reflectionV), std::round(hitInfo.material.shininess)) };
    }
}

//debug purposes
//std::vector<glm::vec3> randomPointsOnLine(glm::vec3 ep0, glm::vec3 ep1, const int& amount) {
std::vector<glm::vec3> randomPointsOnLine(const SegmentLight& segmentLine, const int& amount) {
    glm::vec3 ep0 = segmentLine.endpoint0;
    glm::vec3 ep1 = segmentLine.endpoint1; 
    float xInterval = glm::length(ep1.x - ep0.x) / amount;
    float yInterval = glm::length(ep1.y - ep0.y) / amount;
    float zInterval = glm::length(ep1.z - ep0.z) / amount;

    float xValue = std::min(ep0.x, ep1.x); 
    float yValue = std::min(ep0.y, ep1.y);
    float zValue = std::min(ep0.z, ep1.z);
    std::vector<glm::vec3> randomPoints;
    for (int splits = 0; splits < amount; splits++) {
        xValue = xValue + xInterval; 
        yValue = yValue + yInterval; 
        zValue = zValue + zInterval; 
        glm::vec3 pointOnLine = { xValue, yValue, zValue };
        randomPoints.push_back(pointOnLine); 
    }
    //std::generate(ep0, ep1,  )
    //    std::generate(v.begin(), v.end(), [n = 0, &a]() mutable { return a * n++; });
    //std::uniform_int_distribution<> x_val(ep0.x, ep1.x); 
    //std::uniform_int_distribution<> y_val(ep0.y, ep1.y); 
    //std::uniform_int_distribution<> z_val(ep0.z, ep1.z); 
    return randomPoints; 
}

//NEED TO CHECK ITS ON THE RIGHT SIDE
//NEED TO GET THE COLOR
std::vector<glm::vec3> randomPointsOnParallelogram(const ParallelogramLight& parallelogram, const int& amount) {
    std::vector<glm::vec3> randomPoints; 
    for (int i = 0; i < LINE_PARALLELOGRAM_SAMPLE_LIMIT; i++) {
        float uniformVariate1 = (rand()) / static_cast <float> (RAND_MAX);
        float uniformVariate2 = (rand()) / static_cast <float> (RAND_MAX);
        glm::vec3 randomPoint = parallelogram.v0 + (uniformVariate1 * parallelogram.edge01) + (uniformVariate2 * parallelogram.edge02);
        randomPoints.push_back(randomPoint);
    }
    return randomPoints;
}

bool hardShadow (const Ray& ray, const PointLight& light, BoundingVolumeHierarchy& bvh)
{
    glm::vec3 point_intersection = ray.origin + ray.t * ray.direction;
    glm::vec3 shadowRay_direction = glm::normalize(light.position - point_intersection);
    float lightDistance = glm::distance(point_intersection, light.position);
    Ray shadowRay = { point_intersection + (RAY_STEP * shadowRay_direction), shadowRay_direction, glm::distance(point_intersection, light.position) };
  
    //Draw a red debug ray if the shadow ray hits another object
    HitInfo hitInfo;
    if (bvh.intersect(shadowRay, hitInfo))
    {
        drawRay(shadowRay, glm::vec3{ 1.0f, 0.0f, 0.0f });
        return true;
    }
    //Draw a yellow debug ray if there is no intersection with another object
    drawRay(shadowRay, glm::vec3{ 1.0f, 1.0f, 0.0f });
    return false;
}

glm::vec3 segmentShadow (const HitInfo& hitInfo, Ray& ray, const SegmentLight& light, BoundingVolumeHierarchy& bvh)
{
    glm::vec3 final_lighting = glm::vec3{ 0.0f };
    std::vector<glm::vec3> sample_points = randomPointsOnLine(light, LINE_SAMPLE_LIMIT);

    glm::vec3 ep0 = light.endpoint0;
    glm::vec3 ep1 = light.endpoint1;
    float lightDistance = glm::distance(ep0, ep1); 

    for (glm::vec3 point_position : sample_points)
    {
        glm::vec3 intersection_point = ray.origin + (ray.t * ray.direction);
        float intersection_to_point_dist = glm::distance(point_position, intersection_point);

        float pDistance = glm::distance(point_position, ep1); 
        float alpha = pDistance / lightDistance; 
        glm::vec3 color = ((1 - alpha) * light.color0) + (alpha * light.color1); 
 
        //need to do something with the color here. 
        PointLight point_sample = { point_position, color};
        if (!hardShadow(ray, point_sample, bvh))
        {
            //need to fix this stil. 
            //glm::vec3 diffuseOnly(HitInfo hitInfo, glm::vec3& vertexPos, glm::vec3& normal, glm::vec3 lightPos)
            // //                glm::vec3 diffuseTerm = diffuseOnly(hitInfo, point_intersection, hitInfo.normal, pointLight);
            //glm::vec3 phongSpecularOnly(HitInfo hitInfo, glm::vec3& vertexPos, glm::vec3& normal, PointLight light, glm::vec3& cameraPos)
            final_lighting += color * (diffuseOnly(hitInfo, intersection_point, point_position) + phongSpecularOnly(hitInfo, intersection_point, point_position, ray.origin));
        };
    }

    return final_lighting * (1.0f / float(LINE_SAMPLE_LIMIT));
}

bool checkSide(glm::vec3 v0, glm::vec3 v1) {
    return v0.x <= v1.x && v0.y <= v1.y && v0.z <= v1.z;
}

float closestPoint(glm::vec3 p, glm::vec3 endEdge, glm::vec3 startEdge) {
    glm::vec3 normalisedV = glm::normalize(endEdge - startEdge);
    glm::vec3 closestPoint1 = normalisedV * glm::dot(p, normalisedV);
    return glm::distance(closestPoint1, p);
}

glm::vec3 bilinearInterpolation(glm::vec3 p, float lightDistancev0, float lightDistancev1, const ParallelogramLight& light) {

    float pdistance1 = closestPoint(p, light.edge01, light.v0);
    float alpha = pdistance1 / lightDistancev0;
    glm::vec3 color1 = ((1 - alpha) * light.color0) + (alpha * light.color1);

    float pdistance2 = closestPoint(p, light.edge02, light.v0); 
    float beta = pdistance2 / lightDistancev1;
    glm::vec3 color2 = ((1 - beta) * light.color2) + (beta * light.color3);
    return color1 + color2;
}



glm::vec3 parallelogramShadow (const HitInfo& hitInfo, Ray& ray, const ParallelogramLight& light, BoundingVolumeHierarchy& bvh)
{
    glm::vec3 final_lighting = glm::vec3{ 0.0f };
    std::vector<glm::vec3> sample_points = randomPointsOnParallelogram(light, LINE_SAMPLE_LIMIT);

    float lightDistancev0 = glm::distance(light.edge01, light.v0); 
    float lightDistancev1 = glm::distance(light.edge02, light.v0);

    for (glm::vec3 point_position : sample_points)
    {
        glm::vec3 intersection_point = ray.origin + (ray.t * ray.direction);
        float intersection_to_point_dist = glm::distance(point_position, intersection_point);
        glm::vec3 color = bilinearInterpolation(point_position, lightDistancev0, lightDistancev1, light);

        //need to do something with the color here. 
        PointLight point_sample = { point_position, color };
        if (!hardShadow(ray, point_sample, bvh))
        {
            //need to fix this stil. 
            //glm::vec3 diffuseOnly(HitInfo hitInfo, glm::vec3& vertexPos, glm::vec3& normal, glm::vec3 lightPos)
            // //                glm::vec3 diffuseTerm = diffuseOnly(hitInfo, point_intersection, hitInfo.normal, pointLight);
            //glm::vec3 phongSpecularOnly(HitInfo hitInfo, glm::vec3& vertexPos, glm::vec3& normal, PointLight light, glm::vec3& cameraPos)
            final_lighting += color * (diffuseOnly(hitInfo, intersection_point, point_position) + phongSpecularOnly(hitInfo, intersection_point, point_position, ray.origin));
        };
    }

    return final_lighting * (1.0f / float(LINE_SAMPLE_LIMIT));
}

glm::vec3 lightRay(Ray& ray, HitInfo& hitInfo, const Scene& scene, BoundingVolumeHierarchy& bvh)
{
    glm::vec3 point_intersection = ray.origin + ray.t * ray.direction;
    // Draw a blue debug ray if the ray hit.
    drawRay(ray, glm::vec3{ 0.0f, 0.0f, 1.0f });

    // Calculates the lighting for a point light source
    glm::vec3 lighting = glm::vec3{ 0.0f };
    for (const auto& light : scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            const PointLight pointLight = std::get<PointLight>(light);
            if (!hardShadow(ray, pointLight, bvh))
            {
                //return glm::vec3{ 1.0f }; 
                glm::vec3 diffuseTerm = diffuseOnly(hitInfo, point_intersection, pointLight.position);
                glm::vec3 specularTerm = phongSpecularOnly(hitInfo, point_intersection, pointLight.position, ray.origin);
                lighting = lighting + pointLight.color * (diffuseTerm + specularTerm); 
            }
        }

        if (std::holds_alternative<SegmentLight>(light)) {
            const SegmentLight segmentLight = std::get<SegmentLight>(light);
            lighting = lighting + segmentShadow(hitInfo, ray, segmentLight, bvh); 
        }

        if (std::holds_alternative<ParallelogramLight>(light)) {
            const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
            lighting = lighting + parallelogramShadow(hitInfo, ray, parallelogramLight, bvh);
        }

    }
    return lighting;
}

glm::vec3 recursiveRayTrace(Ray& intersectionRay, HitInfo& hitInfo, const Scene& scene,
    BoundingVolumeHierarchy& bvh, int rayLevel)
{
    if (rayLevel < RECURSION_LIMIT) // Control the recursion level.
    {
        if (hitInfo.material.ks != glm::vec3{ 0.0f }) // Check if current intersection surface is specular.
        {
            // Construct the mirror reflection ray.
            glm::vec3 reverse_incidence_direction = glm::normalize(glm::vec3{ -intersectionRay.direction.x, -intersectionRay.direction.y, -intersectionRay.direction.z }); // TODO: Normalisation could potentially be removed for optimisation.
            glm::vec3 normalised_normal = glm::normalize(hitInfo.normal);
            glm::vec3 reflection_ray_direction = (2.0f * glm::dot(reverse_incidence_direction, normalised_normal) * normalised_normal) - reverse_incidence_direction;
            Ray reflection_ray = {
                glm::vec3{intersectionRay.origin + (intersectionRay.t * intersectionRay.direction) + (RAY_STEP * reflection_ray_direction)}, // Slight offset in the reflection ray direction to prevent intersection with intersectionRay's intersection point.
                reflection_ray_direction,
                std::numeric_limits<float>::max() };

            // Compute new ray lighting data and recursively add lighting data of subsequent mirror reflection rays.
            HitInfo newRayInfo;
            if (bvh.intersect(reflection_ray, newRayInfo))
            {
                return lightRay(intersectionRay, hitInfo, scene, bvh) + (hitInfo.material.ks * recursiveRayTrace(reflection_ray, newRayInfo, scene, bvh, ++rayLevel));
            }
        }
        return lightRay(intersectionRay, hitInfo, scene, bvh);
    }
    return glm::vec3{ 0.0f };
}

