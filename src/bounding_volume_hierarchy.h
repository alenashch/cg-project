#pragma once
#include "ray_tracing.h"
#include "scene.h"
#include <array>
#include <span>

struct Triangle
{
    glm::vec3 center;
    glm::vec3 min;
    glm::vec3 max;

};

class BoundingVolumeHierarchy {

    struct Node {
        int index;
        int level;
        bool isLeaf = false;
        AxisAlignedBox aabb;
        std::vector<Triangle> triangles;
    };

public:
    int MAX_LEVEL = 20;
    int levels = 0.0;
    std::vector<Node> nodes;

    BoundingVolumeHierarchy(Scene* pScene);

    void countLevels(std::vector<Triangle> triangles, int level);

    AxisAlignedBox createAABB(std::vector<Triangle> triangles);

    void constructBVHTree(std::vector<Triangle> triangles, int index, int level, int axis);
    

    // Implement these two functions for the Visual Debug.
    // The first function should return how many levels there are in the tree that you have constructed.
    // The second function should draw the bounding boxes of the nodes at the selected level.
    int numLevels();
    void debugDraw(int level);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo) const;

private:
    Scene* m_pScene;
};
