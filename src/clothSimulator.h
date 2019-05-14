#ifndef CGL_CLOTH_SIMULATOR_H
#define CGL_CLOTH_SIMULATOR_H

#include <nanogui/nanogui.h>

#include "camera.h"
#include "cloth.h"
#include "collision/collisionObject.h"

using namespace nanogui;

struct UserShader;
enum ShaderTypeHint { WIREFRAME = 0, NORMALS = 1, PHONG = 2 };

class ClothSimulator {
public:
    ClothSimulator(std::string project_root);
    ~ClothSimulator();

    void init();

    void loadCloth(Cloth *cloth);
    void loadClothParameters(ClothParameters *cp);
    void loadCollisionObjects(vector<CollisionObject *> *objects);
    virtual bool isAlive();
    void simulateRemotely();

    // File management

    std::string m_project_root;


    // Default simulation values

    int frames_per_sec = 60;//90;
    int simulation_steps = 2;//30;

    CGL::Vector3D gravity = CGL::Vector3D(0, -9.8, 0);

    Cloth *cloth;
    ClothParameters *cp;
    vector<CollisionObject *> *collision_objects;

    bool is_alive = true;

};

#endif // CGL_CLOTH_SIM_H