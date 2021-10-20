#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include <queue>
#include <limits>
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>



BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{

    std::vector<Triangle> triangles;

    // Traverse trhough all the triangles in the scene and add them to the triangles vector
    for (int i = 0; i < m_pScene->meshes.size(); i++) {
        Mesh currentMesh = m_pScene->meshes[i];

        for (int j = 0; j < currentMesh.triangles.size(); j++) {

            glm::vec3 V0 = currentMesh.vertices[currentMesh.triangles[j][0]].position;
            glm::vec3 V1 = currentMesh.vertices[currentMesh.triangles[j][1]].position;
            glm::vec3 V2 = currentMesh.vertices[currentMesh.triangles[j][2]].position;

            float x_min = std::min(V0.x, V1.x);
            x_min = std::min(x_min, V2.x);
            float y_min = std::min(V0.y, V1.y);
            y_min = std::min(y_min, V2.y);
            float z_min = std::min(V0.z, V1.z);
            z_min = std::min(z_min, V2.z);

            float x_max = std::max(V0.x, V1.x);
            x_max = std::max(x_max, V2.x);
            float y_max = std::max(V0.y, V1.y);
            y_max = std::max(y_max, V2.y);
            float z_max = std::max(V0.z, V1.z);
            z_max = std::max(z_max, V2.z);

            Triangle triangle;
            triangle.center = { (V0.x + V1.x + V2.x) / 3, (V0.y + V1.y + V2.y) / 3, (V0.z + V1.z + V2.z) / 3 };     
            triangle.min = { x_min, y_min, z_min };
            triangle.max = { x_max, y_max, z_max };
            triangles.push_back(triangle);
        }

    }
    // Count the levels
    levelCount(triangles, 0);

    Node n;
    n.level = -1;
    nodes.resize(std::pow(2, numberOfLevels + 1) - 1, n);
    constructBVHTree(triangles, 0, 0, 0);
}

/*
* Counts the levels of the BoundingVolumeHierarchy
*/
void BoundingVolumeHierarchy::levelCount(std::vector<Triangle> triangles, int level)
{
    if (numberOfLevels < level + 1) {
        numberOfLevels = level + 1;
    }
    if (triangles.size() != 1 && level != MAX_LEVEL) {

        std::vector<Triangle> childrenLeft(triangles.begin(), triangles.begin() + triangles.size() / 2);
        std::vector<Triangle> childrenRight(triangles.begin() + triangles.size() / 2, triangles.end());

        if(childrenLeft.size() > childrenRight.size())
            levelCount(childrenLeft, level + 1);
        else
            levelCount(childrenRight, level + 1);

    }
}

/*
* Compare the coordinates of centroids of traingles (X,Y,Z)
*/
bool compareCentroidsX(Triangle t1, Triangle t2)
{
    return (t1.center.x < t2.center.x);
}

bool compareCentroidsY(Triangle t1, Triangle t2)
{
    return (t1.center.y < t2.center.y);
}

bool compareCentroidsZ(Triangle t1, Triangle t2)
{
    return (t1.center.z < t2.center.z);
}


/*
* Creates an AABB around the given triangles
*/
AxisAlignedBox BoundingVolumeHierarchy::createAABB(std::vector<Triangle> triangles) {

    float x_min = triangles[0].min[0];
    float y_min = triangles[0].min[1];
    float z_min = triangles[0].min[2];

    float x_max = triangles[0].max[0];
    float y_max = triangles[0].max[1];
    float z_max = triangles[0].max[2];


    for (int i = 1; i < triangles.size(); i++) {
        if (triangles[i].min[0] < x_min)
            x_min = triangles[i].min[0];

        if (triangles[i].min[1] < y_min)
            y_min = triangles[i].min[1];

        if (triangles[i].min[2] < z_min)
            z_min = triangles[i].min[2];


        if (triangles[i].max[0] > x_max)
            x_max = triangles[i].max[0];

        if (triangles[i].max[1] > y_max)
            y_max = triangles[i].max[1];

        if (triangles[i].max[2] > z_max)
            z_max = triangles[i].max[2];
    }

    return { glm::vec3(x_min, y_min, z_min), glm::vec3(x_max, y_max, z_max) };
}


/*
* Constructs the BVH Tree
*/
void BoundingVolumeHierarchy::constructBVHTree(std::vector<Triangle> triangles, int index, int level, int axis) {

    if (triangles.size() == 1 || level == MAX_LEVEL) {

        Node node;
        node.index = index;
        node.level = level;
        node.isLeaf = true;
        node.aabb = createAABB(triangles);
        node.triangles = triangles;
        nodes[index] = node;
    }
    else {

        Node node;
        node.index = index;
        node.level = level;
        node.aabb = createAABB(triangles);
        nodes[index] = node;

        // Sort triangles based on the axis
        if (axis == 0)
            sort(triangles.begin(), triangles.end(), compareCentroidsX);
        if (axis == 1)
            sort(triangles.begin(), triangles.end(), compareCentroidsY);
        if (axis == 2)
            sort(triangles.begin(), triangles.end(), compareCentroidsZ);

        // Change the next axis of division
        int nextAxis = 0;
        if (axis == 0)
            nextAxis = 1;
        if (axis == 1)
            nextAxis = 2;
        if (axis == 2)
            nextAxis = 0;

        // Median division
        std::vector<Triangle> childrenLeft = std::vector<Triangle>(triangles.begin(), triangles.begin() + triangles.size() / 2);
        std::vector<Triangle> childrenRight = std::vector<Triangle>(triangles.begin() + triangles.size() / 2, triangles.end());

        if (childrenLeft.size() != 0) {
            constructBVHTree(childrenLeft, index * 2 + 1, level + 1, nextAxis);
        }
        if (childrenRight.size() != 0) {
            constructBVHTree(childrenRight, index * 2 + 2, level + 1, nextAxis);
        }

    }

}


// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display.
int BoundingVolumeHierarchy::numLevels()
{
    return numberOfLevels;
}

// Use this function to visualize your BVH. This can be useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDraw(int level)
{
    for (int i = 0; i < nodes.size(); i++) {
        if (nodes[i].level == level) {
            //drawAABB(nodes[i].aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
            drawAABB(nodes[i].aabb, DrawMode::Wireframe, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
        }
    }
}

void normalInterpolation(const auto& v0, const auto& v1, const auto& v2, Ray& ray, HitInfo& hitInfo) 
{
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    float bigTriangleArea = glm::length(glm::cross((v1.position - v0.position), (v2.position - v0.position))) / 2.0f;

    //P = w*v0 + u*v1 + v*v2
    float w = (glm::length(glm::cross((p - v1.position), (v2.position - v1.position))) / 2.0f) / bigTriangleArea;
    float u = (glm::length(glm::cross((p - v0.position), (v2.position - v0.position))) / 2.0f) / bigTriangleArea;
    float v = (glm::length(glm::cross((v0.position - v1.position), (p - v1.position))) / 2.0f) / bigTriangleArea;

    hitInfo.normal = glm::normalize(v0.normal * w + v1.normal * u + v2.normal * v);

    // Normal for vector v0
    Ray ray0;
    ray0.origin = v0.position;
    ray0.direction = v0.normal * w;
    ray0.t = 3.0;
    glm::vec3 colour0(1.0, 0.0, 0.0);
    drawRay(ray0, colour0);

    // Normal for vector v1
    Ray ray1;
    ray1.origin = v1.position;
    ray1.direction = v1.normal * u;
    ray1.t = 3.0;
    glm::vec3 colour1(0.0, 1.0, 0.0);
    drawRay(ray1, colour1);

    // Normal for vector v2
    Ray ray2;
    ray2.origin = v2.position;
    ray2.direction = v2.normal * v;
    ray2.t = 3.0;
    glm::vec3 colour2(0.0, 0.0, 1.0);
    drawRay(ray2, colour2);

    // Result
    Ray result;
    result.origin = p;
    result.direction = hitInfo.normal;
    result.t = 3.0;
    glm::vec3 colour = w * colour0 + u * colour1 + v * colour2;
    drawRay(result, colour);
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h .
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo) const
{
    bool hit = false;
    // Intersect with all triangles of all meshes.
    for (const auto& mesh : m_pScene->meshes) {
        for (const auto& tri : mesh.triangles) {
            const auto v0 = mesh.vertices[tri[0]];
            const auto v1 = mesh.vertices[tri[1]];
            const auto v2 = mesh.vertices[tri[2]];
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                normalInterpolation(v0, v1, v2, ray, hitInfo);
               
                hitInfo.material = mesh.material;
                hit = true;
            }
        }
    }
    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres)
        hit |= intersectRayWithShape(sphere, ray, hitInfo);
    return hit;
}
