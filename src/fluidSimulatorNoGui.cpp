#include <cmath>
#include <glad/glad.h>

#include <CGL/vector3D.h>
#include <nanogui/nanogui.h>

#include "fluidSimulatorNoGui.h"
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

using namespace nanogui;
using namespace std;

FluidSimulatorNoGui::~FluidSimulatorNoGui() {
    if (cloth) delete cloth;
    if (cp) delete cp;
    if (collision_objects) delete collision_objects;
}

void FluidSimulatorNoGui::loadCloth(Cloth *cloth) { this->cloth = cloth; }

void FluidSimulatorNoGui::loadClothParameters(ClothParameters *cp) { this->cp = cp; }

void FluidSimulatorNoGui::loadCollisionObjects(vector<CollisionObject *> *objects) { this->collision_objects = objects; }

/**
 * Initializes the cloth simulation and spawns a new thread to separate
 * rendering from simulation.
 */
void FluidSimulatorNoGui::init() {
}

bool FluidSimulatorNoGui::isAlive() { return is_alive; }

void FluidSimulatorNoGui::simulateRemotely() {
    is_paused = false;

    for (int i = 0; i < 100; i++) {
        cout << "Iteration #: " << i << "\n";

        vector<Vector3D> external_accelerations = {gravity};

        for (int i = 0; i < simulation_steps; i++) {
            cloth->simulate(frames_per_sec, simulation_steps, cp, external_accelerations, collision_objects);
        }
    }
    is_paused = true;
}
