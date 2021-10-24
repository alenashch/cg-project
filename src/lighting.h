#pragma once
#include "ray_tracing.h"
#include "draw.h"
#include "bounding_volume_hierarchy.h"
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include "lighting.cpp"	

glm::vec3 diffuseOnly(HitInfo hitInfo, glm::vec3& vertexPos, glm::vec3 lightPos);

glm::vec3 phongSpecularOnly(HitInfo hitInfo, glm::vec3& vertexPos, glm::vec3 lightPos, glm::vec3& cameraPos);

std::vector<glm::vec3> randomPointsOnLine(const SegmentLight& segmentLine, const int& amount); 

std::vector<glm::vec3> randomPointsOnParallelogram(const ParallelogramLight& parallelogram, const int& amount);

glm::vec3 bilinearInterpolation(glm::vec3 p, float lightDistancev0, float lightDistancev1, const ParallelogramLight& light);

float closestPoint(glm::vec3 p, glm::vec3 endEdge, glm::vec3 startEdge);

glm::vec3 segmentShadow(const HitInfo& hitInfo, Ray& ray, const SegmentLight& light, BoundingVolumeHierarchy& bvh);

glm::vec3 parallelogramShadow(const HitInfo& hitInfo, Ray& ray, const ParallelogramLight& light, BoundingVolumeHierarchy& bvh);

bool hardShadow(const Ray& ray, const PointLight& light, BoundingVolumeHierarchy &bvh);

glm::vec3 lightRay(Ray& ray, HitInfo& hitInfo, const Scene& scene, BoundingVolumeHierarchy& bvh);

glm::vec3 recursiveRayTrace(Ray& intersectionRay, HitInfo& hitInfo, const Scene& scene,
    BoundingVolumeHierarchy& bvh, int rayLevel);