#include <cmath>
#include <glad/glad.h>

#include <CGL/vector3D.h>
#include <nanogui/nanogui.h>

#include "clothSimulator.h"
#include "leak_fix.h"

#include "camera.h"
#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"
#include "misc/camera_info.h"
#include "misc/file_utils.h"
// Needed to generate stb_image binaries. Should only define in exactly one source file importing stb_image.h.
#define STB_IMAGE_IMPLEMENTATION
#define GL_POINT_SMOOTH 0x0B10
#define GL_POINT_SMOOTH_HINT 0x0C51
#include "misc/stb_image.h"
#include <ctime>

using namespace nanogui;
using namespace std;


ClothSimulator::ClothSimulator(std::string project_root)
        : m_project_root(project_root) {
}

ClothSimulator::~ClothSimulator() {
    if (cloth) delete cloth;
    if (cp) delete cp;
    if (collision_objects) delete collision_objects;
}

void ClothSimulator::loadCloth(Cloth *cloth) { this->cloth = cloth; }

void ClothSimulator::loadClothParameters(ClothParameters *cp) { this->cp = cp; }

void ClothSimulator::loadCollisionObjects(vector<CollisionObject *> *objects) { this->collision_objects = objects; }

/**
 * Initializes the cloth simulation and spawns a new thread to separate
 * rendering from simulation.
 */
void ClothSimulator::init() {
}

bool ClothSimulator::isAlive() { return is_alive; }

void ClothSimulator::simulateRemotely() {
    clock_t begin = clock();

    for (int j = 0; j < 100; j++) {
        cout << "Iteration #: " << j << "\n";

        vector<Vector3D> external_accelerations = {gravity};

        for (int i = 0; i < simulation_steps; i++) {
            cloth->simulate(frames_per_sec, simulation_steps, cp, external_accelerations, collision_objects);
        }
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "This round of 100 iterations took: " << elapsed_secs << " seconds = " << (elapsed_secs / 60) << " minutes = " << (elapsed_secs / 60 / 60) << " hours.\n";
}
