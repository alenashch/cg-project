#pragma once
#include "ray_tracing.h"
#include "draw.h"
#include "bounding_volume_hierarchy.h"
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include "lighting.cpp"	
// Returns a boolean indicating if the ray's intersection point is occluded from the light, and therefore is in shadow
// with regards to that light.
bool hardShadow(const Ray& ray, const PointLight& light, BoundingVolumeHierarchy &bvh);

glm::vec3 lightRay(const Ray& ray, const HitInfo& hitInfo, const Scene& scene, BoundingVolumeHierarchy& bvh);