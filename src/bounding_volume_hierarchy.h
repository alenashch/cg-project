#pragma once
#include "ray_tracing.h"
#include "scene.h"
#include <array>
#include <span>
#include <queue>

struct Triangle
{
    glm::vec3 center;
    glm::vec3 min;
    glm::vec3 max;
    Vertex V0;
    Vertex V1;
    Vertex V2;

    Material material;

};

struct rayAABB {
    Ray ray;
    AxisAlignedBox AABB;
    float hitT{ std::numeric_limits<float>::max() };
};

class BoundingVolumeHierarchy {

    struct Node {
        int index;
        int level;
        bool isLeaf = false;
        AxisAlignedBox aabb;
        std::vector<Triangle> triangles;
        float hitT{ std::numeric_limits<float>::max() };
    };

    struct compare {
        bool operator() (const Node& a, const Node& b) {
            return a.hitT < b.hitT;
        }
    };

public:
    int MAX_LEVEL = 20;
    int numberOfLevels = 0.0;
    std::vector<Node> nodes;

    BoundingVolumeHierarchy(Scene* pScene);

    void setHitT(Node node, float t) const;

    void levelCount(std::vector<Triangle> triangles, int level);

    AxisAlignedBox createAABB(std::vector<Triangle> triangles);

    void constructBVHTree(std::vector<Triangle> triangles, int index, int level, int axis);

    void nodeIntersection(std::vector<Node> nodes, Ray& ray, HitInfo& hitInfo, std::priority_queue<Node, std::vector<Node>, compare> rayAABBintersections) const;
    
   // bool BoundingVolumeHierarchy::compare(Node a, Node b);
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
