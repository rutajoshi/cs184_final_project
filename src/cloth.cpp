#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"


using namespace std;

#define BOUNCE_DAMPING_FACTOR 0.9

static int count_steps;

bool check_vector(Vector3D vec) {
    bool result = true;
    result &= (!isnan(vec.x));
    result &= (!isinf(vec.x));

    result &=  (!isnan(vec.y));
    result &=  (!isinf(vec.y));

    result &=  (!isnan(vec.z));
    result &=  (!isinf(vec.z)); 
    return result;  
}

Cloth::Cloth(double width, double height, double depth, int num_width_points,
             int num_height_points, int num_depth_points, float thickness) {

    this->width = width;
    this->height = height;
    this->depth = depth;
    this->num_width_points = num_width_points;
    this->num_height_points = num_height_points;
    this->num_depth_points = num_depth_points;
    this->thickness = thickness;
    this->epsilon = 0.1; //0.01;

    buildGrid();
    buildClothMesh();
}

Cloth::~Cloth() {
    point_masses.clear();
    springs.clear();

    if (clothMesh) {
        delete clothMesh;
    }
}

void Cloth::buildGrid() {
    // TODO (Part 1): Build a grid of masses.
    // Make all the point masses
    double start_left = 0;
    double start_right = 0;
    double start_bottom = 0;
    for (int p = 0; p < num_depth_points; p++) {
        for (int c = 0; c < num_width_points; c++) {
            for (int r = 0; r < num_height_points; r++) {
                double x = start_left + r * (height / (num_height_points - 1));
                double y = start_right + c * (width / (num_width_points - 1));
                double z = start_bottom + p * (depth / (num_depth_points - 1));

                Vector3D pos;
                if (orientation == HORIZONTAL) {
                    //pos = Vector3D(x, 1, y);
                    pos = Vector3D(x, z, y);
                }
                else {
                    //double offset = ((rand() % 2) - 1.0) / 1000;
                    //pos = Vector3D(x, y, offset);
                    pos = Vector3D(x, y, z);
                }
                //vector<int> rc{ r, c };
                bool pinnedPoint = false;
                PointMass pm = PointMass(pos, pinnedPoint);
                // NOTE: density is from pinned2.json
                pm.mass = width * height * depth * 1000 / num_width_points / num_height_points / num_depth_points;
                point_masses.push_back(pm);
            }
        }
    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
    double mass = width * height * depth * cp->density / num_width_points / num_height_points / num_depth_points;
    double delta_t = 1.0f / frames_per_sec / simulation_steps; // 0.016; //
    count_steps += 1;

    Vector3D total_ext_force = Vector3D(0,0,0);
    for (int i = 0; i < external_accelerations.size(); i++) {
        total_ext_force += external_accelerations[i] * mass;
    }

    for (PointMass &pm : point_masses) {
        pm.forces += total_ext_force;
    }

    for (PointMass &pm : point_masses) {
        // For all particles i, apply forces to get a velocity update.
        // Then set the predicted position using only the velocity.
        pm.velocity += (pm.forces / mass) * delta_t;
        if (pm.collision_forces.norm() != 0) {
            // Remove collision forces from the last iteration so they only apply for one time step
            pm.collision_forces *= BOUNCE_DAMPING_FACTOR;
        }
        pm.predict_position = pm.position + delta_t * pm.velocity;
        // Set the last position to the current predicted position for use in hashing to find neighbors
        pm.last_position = pm.predict_position;

        // validity check the velocity and predicted position
        assert(check_vector(pm.velocity));
        assert(check_vector(pm.predict_position));
    }

    // Build a spatial map so you can easily find all the neighbors of a particle
    build_spatial_map();

    // For some number of iterations, do logic to update the predicted position
    for (int j = 0; j < 5; j++) {
        // For each particle, calculate lambda_i
        for (PointMass &pm : point_masses) {
            assert(check_vector(pm.predict_position));
            lambda_i(pm);
            assert(check_vector(pm.predict_position));

            // validity check the lambda
            assert(!isnan(pm.lambda));
            assert(!isinf(pm.lambda));

        }

        // For each particle, deal with self-collisions (repel it from other point masses)
        for (PointMass &pm : point_masses) {
            assert(check_vector(pm.predict_position));
            self_collide(pm, simulation_steps);
            assert(check_vector(pm.predict_position));
        }

        // For each particle, calculate delta position and do collision detection/response
        for (PointMass &pm : point_masses) {
            assert(check_vector(pm.predict_position));
            pm.delta_position = calculate_delta_p(pm);// CALCULATE delta_p here

            // Collide with all collision objects
            for (CollisionObject* c : *collision_objects){
                c->collide(pm);
            }

            pm.forces = pm.collision_forces;

            assert(check_vector(pm.predict_position));

            // test
            assert(check_vector(pm.delta_position));

        }

        // For each particle, update the predicted position using delta position
        for (PointMass &pm : point_masses) {
            // Update predicted positions
            pm.predict_position += pm.delta_position;

            // validity test the predicted positions
            assert(check_vector(pm.predict_position));
        }
    }

    for (PointMass &pm : point_masses) {
        // Update velocity
        assert(check_vector(pm.velocity));
        pm.velocity = (1.0 / delta_t) * (pm.predict_position - pm.position);
        assert(check_vector(pm.velocity));


        // Update vorticity
//        Vector3D force_vort = force_vorticity_i(pm);
//        pm.velocity += (force_vort / mass) * delta_t;
//
//        assert(check_vector(pm.velocity));
//
//        // Update viscosity (happens in place)
//        viscosity_constraint(pm);
//        assert(check_vector(pm.velocity));

        // Update position
        pm.position = pm.predict_position;

        // test
        assert(check_vector(pm.position));
    }

}


void Cloth::build_spatial_map() {
    for (const auto &entry : map) {
        delete(entry.second);
    }
    map.clear();

    // TODO (Part 4): Build a spatial map out of all of the point masses.
    for (PointMass &pm : point_masses) {
        float hash_pos = hash_position(pm.last_position);
        if (map.find(hash_pos) == map.end()) {
            // Insert a new pair
            vector<PointMass *> *vect = new vector<PointMass *>();
            vect->push_back(&pm);
            map.insert(std::pair<float, vector<PointMass *> *>(hash_pos, vect));
        } else {
            // Add this position to the list that goes with hash_pos
            map.find(hash_pos)->second->push_back(&pm);
        }
    }
}


// ##########################################################
double Cloth::calc_h() {
    double s = 1.0 / pow(num_width_points * num_height_points * num_depth_points, 1.0/3.0);
    return s;
//    return 0.25;
}


Vector3D Cloth::calculate_delta_p(PointMass &pm_i) {
    float hash_pos = hash_position(pm_i.last_position);

    // test
    if (isnan(hash_pos)) {
        std::cout << "\nHash pos is nan for position = " << pm_i.position << "\n";
    }

    auto getter = map.find(hash_pos);
    vector<PointMass *> *neighbors = getter->second;
    double h = calc_h();

    Vector3D delta_p = Vector3D(0,0,0);

    // Tensile instability constants
    double k = 0.2; //0.1;
    int n = 4;
    Vector3D delta_q = 0.2 * Vector3D(0.25, 0.25, 0.25); //Vector3D(0.03, 0.03, 0.03);
    double denom = kernel_poly6(delta_q, h);

    assert(check_vector(pm_i.predict_position));

    for (PointMass *neighbor : *neighbors) {
        if (neighbor == &pm_i) {
            continue;
        }
        Vector3D neighborToPm = (pm_i.predict_position - neighbor->predict_position);
        assert(check_vector(neighbor -> predict_position));
        assert(check_vector(neighborToPm));
        Vector3D term = spiky_kernel_grad(neighborToPm, h);


        // Tensile instability calculations
        double numer = kernel_poly6(neighborToPm, h);
        double s_corr = -k * pow(numer / denom, n);

        delta_p += (pm_i.lambda + neighbor->lambda + s_corr) * term;
    }
    delta_p = (1.0 / pm_i.rest_density) * delta_p;
    return delta_p;

}

double Cloth::kernel_poly6(Vector3D pos_dif, double h) {
    double r = pos_dif.norm();
    assert(!isnan(r));
    assert(!isinf(r));

    if (0 <= r && r <= h) {
        assert(r <= h);
        double mult = pow((pow(h,2) - pow(r,2)), 3);
        return 315. / 64. / M_PI / pow(h,9) * mult;
    }
    if (r <= h) {
        std::cout << "\n kernel_poly6 r= " << r<<" h = " <<h<<"\n"; 
    }
    assert(r > h);
    return 0;
}

double Cloth::calculate_density_neighbors(PointMass &pm) {
    float hash_pos = hash_position(pm.last_position);

    // validity test the hash position
    if (isnan(hash_pos)) {
        std::cout << "\nHash position is nan for position = " << pm.position << "\n";
    }

    auto getter = map.find(hash_pos);
    vector<PointMass *> *neighbors = getter->second;
    double h = calc_h();
    double sum = 0;
    for (PointMass *neighbor : *neighbors) {
        if (neighbor == &pm) {
            continue;
        }
        Vector3D neighborToPm = (pm.predict_position - neighbor->predict_position);
        sum += pm.mass * kernel_poly6(neighborToPm, h);
    }
    return sum;
}

Vector3D Cloth::spiky_kernel_grad(Vector3D pos_dif, double h) {
    double r = pos_dif.norm();
    assert(!isnan(r));
    assert(!isinf(r));

    if (0 < r && r <= h) {
        Vector3D r_hat = pos_dif / r;
        Vector3D result = (-45.0 / (M_PI * pow(h, 6))) * pow(h-r, 2) * r_hat;
        return result;
    }

    return Vector3D(0);

}

double Cloth::viscosity_kernel(Vector3D pos_dif, double h) {
    double r = pos_dif.norm();
    assert(!isnan(r));
    assert(!isinf(r));

    if (0 <= r && r <= h) {
        double mult = 45.0 / M_PI / pow(h, 6) * (h - r);
        return  mult;
    }
    assert( r > h);
    return 0.0;
}

double Cloth::delta_constraint_pk(PointMass &pm_i, PointMass &pm_k) {
    float hash_pos = hash_position(pm_i.last_position);
    auto getter = map.find(hash_pos);
    vector<PointMass *> *neighbors = getter->second;
    double h = calc_h();
    Vector3D sum; 

    if (&pm_k == &pm_i) {
        // k = i
        sum = Vector3D();
        for (PointMass *neighbor : *neighbors) {
            if (neighbor == &pm_i) {
                continue;
            } else {
                Vector3D neighborToPm = (pm_i.predict_position - neighbor->predict_position);
                assert(check_vector(neighborToPm));
                sum += spiky_kernel_grad(neighborToPm, h);

            }
        }
        sum /= (pm_i.rest_density);

    } else {
        // k = j
        Vector3D pi_pk = (pm_i.predict_position - pm_k.predict_position);

        assert(check_vector(pi_pk));

        sum = -1.0 * (1.0 / pm_i.rest_density) * spiky_kernel_grad(pi_pk, h);
        
    }
    double sum_norm = sum.norm();
    assert(!(isinf(sum_norm)));
    // assert(!(isnan(sum_norm)));
    if(isnan(sum_norm)) {
        std::cout << "\n pm_i " << pm_i.predict_position <<"\n pm_k" << pm_k.predict_position <<"\n" << "sum_norm " << sum <<"\n";
    }
    return sum_norm;

}

void Cloth::lambda_i(PointMass &pm) {
    double rho_i = calculate_density_neighbors(pm);
    double rho_o = pm.rest_density;
    double C_i = (rho_i / rho_o) - 1.0;

    assert(!isnan(rho_i) && !(isinf(rho_i)));

    double denom = 0.0;
    float hash_pos = hash_position(pm.last_position);
    auto getter = map.find(hash_pos);
    vector<PointMass *> *neighbors = getter->second;
    for (PointMass *neighbor : *neighbors) {
        double grad_Ci_pk_norm = delta_constraint_pk(pm, *neighbor);
        
        assert(!(isinf(grad_Ci_pk_norm)));
        assert(!isnan(grad_Ci_pk_norm));
        denom += pow(grad_Ci_pk_norm, 2);
    }
    denom += this->epsilon;

    assert(!isnan(denom) && !(isinf(denom)));

    double lambda = -1.0 * C_i / denom;

    pm.lambda = lambda;
}

Vector3D Cloth::vorticity_wi(PointMass &pm_i) {
    float hash_pos = hash_position(pm_i.last_position);
    auto getter = map.find(hash_pos);
    vector<PointMass *> *neighbors = getter->second;
    double h = calc_h();

    Vector3D w_i = Vector3D();
    for (PointMass *neighbor : *neighbors) {
        if (neighbor == &pm_i) {
            continue;
        } else {
            Vector3D neighborToPm = (pm_i.predict_position - neighbor->predict_position);
            assert(check_vector(neighborToPm));
            w_i += cross((pm_i.velocity - neighbor->velocity), (spiky_kernel_grad(neighborToPm, h)));
        }
    }

    return w_i;
}

Vector3D Cloth::location_vector(PointMass &pm_i) {
    // note: uses Bubbles Alive, Hong et al [2008]

    float hash_pos = hash_position(pm_i.last_position);
    auto getter = map.find(hash_pos);
    vector<PointMass *> *neighbors = getter->second;
    double h = calc_h();

    Vector3D p_pos_sum = Vector3D();
    for (PointMass *neighbor : *neighbors) {
        p_pos_sum += neighbor->predict_position;
    }

    Vector3D n = p_pos_sum / (neighbors->size() + epsilon) - pm_i.predict_position;

    n.normalize();

    return n;
}

Vector3D Cloth::force_vorticity_i(PointMass &pm_i) {
    Vector3D w_i = vorticity_wi(pm_i);
    Vector3D location_i = location_vector(pm_i);

    Vector3D force_i = cross(location_i, w_i);

    // TODO epsilon: user specified relaxation parameter
    double vort_epsilon = 0.0006;
    force_i *= vort_epsilon;

    return force_i;
}

void Cloth::viscosity_constraint(PointMass &pm_i) {
    double c = 0.01;

    float hash_pos = hash_position(pm_i.last_position);
    auto getter = map.find(hash_pos);
    vector<PointMass *> *neighbors = getter->second;
    double h = calc_h();

    Vector3D viscosity_sum = 0;
    for (PointMass *neighbor : *neighbors) {
        if (neighbor == &pm_i) {
            continue;
        } else {
            Vector3D neighborToPm = (pm_i.predict_position - neighbor->predict_position);
            double visc_kernel = viscosity_kernel(neighborToPm, h);
//          std::cout << "Distance between neighbors = " << neighborToPm.norm() << std::endl;
//        std::cout << "Viscosity kernel = " << visc_kernel << std::endl;
            viscosity_sum += (pm_i.velocity - neighbor->velocity) * visc_kernel;
        }
    }
//  std::cout << "Velocity before update = " << pm_i.velocity << std::endl;
    pm_i.velocity = pm_i.velocity + c*viscosity_sum;
//    std::cout << "Velocity after update = " << pm_i.velocity << std::endl;
}

// #######################################


void Cloth::self_collide(PointMass &pm, double simulation_steps) {
    // TODO (Part 4): Handle self-collision for a given point mass.
    float hash_pos = hash_position(pm.last_position);
    auto getter = map.find(hash_pos);
    vector<PointMass *> *neighbors = getter->second;
    Vector3D totalCorrection = Vector3D(0,0,0);
    int num_corrections = 0;

    for (PointMass *neighbor : *neighbors) {
        if (neighbor == &pm) {
            continue;
        }
        Vector3D neighborToPm = (pm.predict_position - neighbor->predict_position);
        assert(check_vector(neighborToPm));
        float dist = neighborToPm.norm();

        assert(check_vector(neighborToPm));

        /// Problem if both of the particles are in the same location
        if (dist == 0) {
            double small_e = 0.00001;
            neighborToPm = Vector3D();
            neighborToPm.x = totalCorrection.x + small_e;
            neighborToPm.y = totalCorrection.y + small_e;
            neighborToPm.z = totalCorrection.z + small_e;
        } 

        neighborToPm.normalize();
        assert(check_vector(neighborToPm));

        if (dist < 2 * thickness) {
            assert(check_vector(neighbor->predict_position));
            assert(!isnan(thickness));
            Vector3D corrected = neighbor->predict_position + (2 * thickness)*neighborToPm;
            assert(check_vector(corrected));
            Vector3D correction = corrected - pm.predict_position;
            assert(check_vector(correction));
            totalCorrection += correction;
            num_corrections += 1;
        }
    }

    if (num_corrections > 0) {
        assert(check_vector(totalCorrection));
        totalCorrection /= (num_corrections * simulation_steps);
        assert(check_vector(totalCorrection));
        assert(check_vector(pm.predict_position));
        pm.predict_position += totalCorrection;
        assert(check_vector(pm.predict_position));
    }
}

float Cloth::hash_position(Vector3D pos) {
    // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    double s = 1.0 / pow(num_width_points * num_height_points * num_depth_points, 1.0/3.0);
    Vector3D hash_indices = Vector3D(floor(pos.x / s), floor(pos.y / s), floor(pos.z / s));
    float hash_value = hash_indices.x * pow(3, 3) + hash_indices.y * pow(3, 2) + hash_indices.z * pow(3, 1);
    return hash_value;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
    PointMass *pm = &point_masses[0];
    for (int i = 0; i < point_masses.size(); i++) {
        pm->position = pm->start_position;
        pm->predict_position = pm->start_position;
        pm->last_position = pm->start_position;
        pm->velocity = Vector3D(0,0,0);
        pm++;
    }
}

void Cloth::buildClothMesh() {
    if (point_masses.size() == 0) return;

    ClothMesh *clothMesh = new ClothMesh();
    vector<Triangle *> triangles;

    // Create vector of triangles
    for (int y = 0; y < num_height_points - 1; y++) {
        for (int x = 0; x < num_width_points - 1; x++) {
            PointMass *pm = &point_masses[y * num_width_points + x];
            // Get neighboring point masses:
            /*                      *
             * pm_A -------- pm_B   *
             *             /        *
             *  |         /   |     *
             *  |        /    |     *
             *  |       /     |     *
             *  |      /      |     *
             *  |     /       |     *
             *  |    /        |     *
             *      /               *
             * pm_C -------- pm_D   *
             *                      *
             */

            float u_min = x;
            u_min /= num_width_points - 1;
            float u_max = x + 1;
            u_max /= num_width_points - 1;
            float v_min = y;
            v_min /= num_height_points - 1;
            float v_max = y + 1;
            v_max /= num_height_points - 1;

            PointMass *pm_A = pm                       ;
            PointMass *pm_B = pm                    + 1;
            PointMass *pm_C = pm + num_width_points    ;
            PointMass *pm_D = pm + num_width_points + 1;

            Vector3D uv_A = Vector3D(u_min, v_min, 0);
            Vector3D uv_B = Vector3D(u_max, v_min, 0);
            Vector3D uv_C = Vector3D(u_min, v_max, 0);
            Vector3D uv_D = Vector3D(u_max, v_max, 0);


            // Both triangles defined by vertices in counter-clockwise orientation
            triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
                                             uv_A, uv_C, uv_B));
            triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
                                             uv_B, uv_C, uv_D));
        }
    }

    // For each triangle in row-order, create 3 edges and 3 internal halfedges
    for (int i = 0; i < triangles.size(); i++) {
        Triangle *t = triangles[i];

        // Allocate new halfedges on heap
        Halfedge *h1 = new Halfedge();
        Halfedge *h2 = new Halfedge();
        Halfedge *h3 = new Halfedge();

        // Allocate new edges on heap
        Edge *e1 = new Edge();
        Edge *e2 = new Edge();
        Edge *e3 = new Edge();

        // Assign a halfedge pointer to the triangle
        t->halfedge = h1;

        // Assign halfedge pointers to point masses
        t->pm1->halfedge = h1;
        t->pm2->halfedge = h2;
        t->pm3->halfedge = h3;

        // Update all halfedge pointers
        h1->edge = e1;
        h1->next = h2;
        h1->pm = t->pm1;
        h1->triangle = t;

        h2->edge = e2;
        h2->next = h3;
        h2->pm = t->pm2;
        h2->triangle = t;

        h3->edge = e3;
        h3->next = h1;
        h3->pm = t->pm3;
        h3->triangle = t;
    }

    // Go back through the cloth mesh and link triangles together using halfedge
    // twin pointers

    // Convenient variables for math
    int num_height_tris = (num_height_points - 1) * 2;
    int num_width_tris = (num_width_points - 1) * 2;

    bool topLeft = true;
    for (int i = 0; i < triangles.size(); i++) {
        Triangle *t = triangles[i];

        if (topLeft) {
            // Get left triangle, if it exists
            if (i % num_width_tris != 0) { // Not a left-most triangle
                Triangle *temp = triangles[i - 1];
                t->pm1->halfedge->twin = temp->pm3->halfedge;
            } else {
                t->pm1->halfedge->twin = nullptr;
            }

            // Get triangle above, if it exists
            if (i >= num_width_tris) { // Not a top-most triangle
                Triangle *temp = triangles[i - num_width_tris + 1];
                t->pm3->halfedge->twin = temp->pm2->halfedge;
            } else {
                t->pm3->halfedge->twin = nullptr;
            }

            // Get triangle to bottom right; guaranteed to exist
            Triangle *temp = triangles[i + 1];
            t->pm2->halfedge->twin = temp->pm1->halfedge;
        } else {
            // Get right triangle, if it exists
            if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
                Triangle *temp = triangles[i + 1];
                t->pm3->halfedge->twin = temp->pm1->halfedge;
            } else {
                t->pm3->halfedge->twin = nullptr;
            }

            // Get triangle below, if it exists
            if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
                Triangle *temp = triangles[i + num_width_tris - 1];
                t->pm2->halfedge->twin = temp->pm3->halfedge;
            } else {
                t->pm2->halfedge->twin = nullptr;
            }

            // Get triangle to top left; guaranteed to exist
            Triangle *temp = triangles[i - 1];
            t->pm1->halfedge->twin = temp->pm2->halfedge;
        }

        topLeft = !topLeft;
    }

    clothMesh->triangles = triangles;
    this->clothMesh = clothMesh;
}
