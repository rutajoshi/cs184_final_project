#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"


using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;
  this->epsilon = 0.01;

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
  for (int c = 0; c < num_width_points; c++) {
      for (int r = 0; r < num_height_points; r++) {
          double x = start_left + r * (height / (num_height_points - 1));
          double y = start_right + c * (width / (num_width_points - 1));

          Vector3D pos;
          if (orientation == HORIZONTAL) {
              pos = Vector3D(x, 1, y);
          }
          else {
              double z = ((rand() % 2) - 1.0) / 1000;
              pos = Vector3D(x, y, z);
          }
          vector<int> rc{ r, c };
          bool pinnedPoint = false;

          point_masses.push_back(PointMass(pos, pinnedPoint));
      }
  }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  Vector3D total_ext_force = Vector3D(0,0,0);
  for (int i = 0; i < external_accelerations.size(); i++) {
      total_ext_force += external_accelerations[i] * mass;
  }

  for (PointMass &pm : point_masses) {
      pm.forces = total_ext_force;
  }

  for (PointMass &pm : point_masses) {
      pm.velocity += (pm.forces / mass) * delta_t;
      pm.predict_position = pm.position + delta_t * pm.velocity;

      // test
      assert(!isnan(pm.velocity.x) && !isinf(pm.velocity.x));
      assert(!isnan(pm.predict_position.x) && !isinf(pm.predict_position.x));

  }


  build_spatial_map();
  for (int j = 0; j < 4; j++) {
    for (PointMass &pm : point_masses) {
      lambda_i(pm);

      // test
      assert(!isnan(pm.lambda) && !isinf(pm.lambda));
      
    }


    for (PointMass &pm : point_masses) {
      pm.delta_position = calculate_delta_p(pm);// CALCULATE delta_p here

      // test
      assert(!isnan(pm.delta_position.y) && !isinf(pm.delta_position.y));

      // Collision detection and response
      for (CollisionObject* c : *collision_objects){
        c->collide(pm);
      }
    }

    for (PointMass &pm : point_masses) {
      // Update predicted positions
      pm.predict_position += pm.delta_position;

      // test
      assert(!isnan(pm.predict_position.x) && !isinf(pm.predict_position.x));
    }
  }

  for (PointMass &pm : point_masses) {
    // Update velocity
    pm.velocity = (1.0 / delta_t) * (pm.predict_position - pm.position);
    assert(!isnan(pm.velocity.x) && !isinf(pm.velocity.x));
    

    // Update vorticity
    Vector3D force_vort = force_vorticity_i(pm);
    pm.velocity += (force_vort / mass) * delta_t;

    assert(!isnan(pm.velocity.x) && !isinf(pm.velocity.x));

    // Update viscosity (happens in place)
    viscosity_constraint(pm);
    assert(!isnan(pm.velocity.x) && !isinf(pm.velocity.x));

    // Update position
    pm.position = pm.predict_position;

    // test
    assert(!isnan(pm.position.y) && !isinf(pm.position.y));
  }

}


void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.
  for (PointMass &pm : point_masses) {
      float hash_pos = hash_position(pm.position);
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
Vector3D Cloth::calculate_delta_p(PointMass &pm_i) {
  float hash_pos = hash_position(pm_i.position);

  // test
  if (isnan(hash_pos)) {
      std::cout << "\nHash pos is nan for position = " << pm_i.position << "\n";
  }

  auto getter = map.find(hash_pos);
  vector<PointMass *> *neighbors = getter->second;
  double h = 3 * width / num_width_points / 2;

  Vector3D delta_p = Vector3D();

  for (PointMass *neighbor : *neighbors) {
      if (neighbor == &pm_i) {
          continue;
      }
      Vector3D neighborToPm = (pm_i.position - neighbor->position);
      Vector3D term = spiky_kernel_grad(neighborToPm, h);

      // test
      if (isnan(term.x)) {
          std::cout << "\nSpiky kernel x is nan for position = " << pm_i.position << "\n";
      }

      delta_p = delta_p + (pm_i.lambda + neighbor->lambda) * term;

      // test
      if (isnan(pm_i.lambda)) {
          std::cout << "\nLambda is nan for position = " << pm_i.position << "\n";
      }
      // test
      if (isnan(neighbor->lambda)) {
          std::cout << "\nNeighbor lambda is nan for position = " << pm_i.position << "\n";
      }

      // test
      if (isnan(delta_p.y)) {
          std::cout << "\nDelta position is nan for position = " << pm_i.position << "\n";
      }
  }
  delta_p = (1.0 / pm_i.rest_density) * delta_p;
  return delta_p;

}

double Cloth::kernel_poly6(Vector3D pos_dif, double h) {
  double r = pos_dif.norm();
  if (0 <= r && r <= h) {
    double mult = pow((pow(h,2) - pow(r,2)), 3);
    return 315. / 64. / M_PI / pow(h,9) * mult;
  }
  return 0;
}

double Cloth::calculate_density_neighbors(PointMass &pm) {
  float hash_pos = hash_position(pm.position);

  // test
  if (isnan(hash_pos)) {
    std::cout << "\nHash position is nan for position = " << pm.position << "\n";
  }

  auto getter = map.find(hash_pos);
  vector<PointMass *> *neighbors = getter->second;
  double h = 3 * width / num_width_points / 2;
  double sum = 0;
  for (PointMass *neighbor : *neighbors) {
      if (neighbor == &pm) {
          continue;
      }
      Vector3D neighborToPm = (pm.position - neighbor->position);
      sum += kernel_poly6(neighborToPm, h);
  }
  return sum;
}

Vector3D Cloth::spiky_kernel_grad(Vector3D pos_dif, double h) {
  double r = pos_dif.norm();
  if (0 <= r && r <= h) {
    double mult = -45. /M_PI / pow(h,6) * pow((h - r), 2);
    return mult * pos_dif;
  }
  return Vector3D(0);

}

Vector3D Cloth::viscosity_kernel(Vector3D pos_dif, double h) {
  double r = pos_dif.norm();
  double operand = -pow(r,3)/(2*pow(h,3)) + pow(r,2) / pow(h,2) + h / (2*r) - 1;

  if (0 <= r && r <= h) {
    double mult = -15. / 2 / M_PI / pow(h,3);
    return mult * operand;
  }
  return Vector3D(0);

}

Vector3D Cloth::delta_constraint_pk(PointMass &pm_i, PointMass &pm_k) {
  float hash_pos = hash_position(pm_i.position);
  auto getter = map.find(hash_pos);
  vector<PointMass *> *neighbors = getter->second;
  double h = 3 * width / num_width_points / 2;

  if (&pm_k == &pm_i) {
    // k = i
    Vector3D sum = Vector3D();
    for (PointMass *neighbor : *neighbors) {
        if (neighbor == &pm_i) {
            continue;
        } else{
          Vector3D neighborToPm = (pm_i.position - neighbor->position);
          sum += spiky_kernel_grad(neighborToPm, h);

        }
    }
    sum /= (pm_i.rest_density);
    return sum;

  } else {
    // k = j
    Vector3D pi_pk = (pm_i.position - pm_k.position);
    return -1 * (1.0 / pm_i.rest_density) * spiky_kernel_grad(pi_pk, h);
  }

}

void Cloth::lambda_i(PointMass &pm) {
  double rho_i = calculate_density_neighbors(pm);
  double rho_o = pm.rest_density;
  double C_i = (rho_i / rho_o) - 1.0;

  double denom = 0.0;
  float hash_pos = hash_position(pm.position);
  auto getter = map.find(hash_pos);
  vector<PointMass *> *neighbors = getter->second;
  for (PointMass *neighbor : *neighbors) {
    Vector3D grad_Ci_pk = delta_constraint_pk(pm, *neighbor);
    denom += pow(grad_Ci_pk.norm(), 2);
  }

  denom += this->epsilon;

  double lambda = -1.0 * C_i / denom;

  pm.lambda = lambda;
}

Vector3D Cloth::vorticity_wi(PointMass &pm_i) {
  float hash_pos = hash_position(pm_i.position);
  auto getter = map.find(hash_pos);
  vector<PointMass *> *neighbors = getter->second;
  double h = 3 * width / num_width_points / 2;

  Vector3D w_i = Vector3D();
  for (PointMass *neighbor : *neighbors) {
      if (neighbor == &pm_i) {
          continue;
      } else {
        Vector3D neighborToPm = (pm_i.position - neighbor->position);
        w_i += cross((pm_i.velocity - neighbor->velocity), (spiky_kernel_grad(neighborToPm, h)));
      }
  }

  return w_i;
}

Vector3D Cloth::location_vector(PointMass &pm_i) {
  // note: uses Bubbles Alive, Hong et al [2008]

  float hash_pos = hash_position(pm_i.position);
  auto getter = map.find(hash_pos);
  vector<PointMass *> *neighbors = getter->second;
  double h = 3 * width / num_width_points / 2;

  Vector3D p_pos_sum = Vector3D();
  for (PointMass *neighbor : *neighbors) {
    p_pos_sum += neighbor->position;
  }

  Vector3D n = p_pos_sum / (neighbors->size() + epsilon) - pm_i.position;

  n.normalize();

  return n;
}

Vector3D Cloth::force_vorticity_i(PointMass &pm_i) {
  Vector3D w_i = vorticity_wi(pm_i);
  Vector3D location_i = location_vector(pm_i);

  Vector3D force_i = cross(location_i, w_i);

  // TODO epsilon: user specified relaxation parameter
  // force_i *= epsilon;

  return force_i;
}

void Cloth::viscosity_constraint(PointMass &pm_i) {
  double c = 0.01;

  float hash_pos = hash_position(pm_i.position);
  auto getter = map.find(hash_pos);
  vector<PointMass *> *neighbors = getter->second;
  double h = 3 * width / num_width_points / 2;

  Vector3D viscosity_sum = 0;
  for (PointMass *neighbor : *neighbors) {
      if (neighbor == &pm_i) {
          continue;
      } else {
        Vector3D neighborToPm = (pm_i.position - neighbor->position);
        viscosity_sum += cross((pm_i.velocity - neighbor->velocity), viscosity_kernel(neighborToPm, h));
      }
  }

  pm_i.velocity = pm_i.velocity + c*viscosity_sum;
}

// #######################################


void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.
  if (pm.pinned) {
      return;
  }

  float hash_pos = hash_position(pm.position);
  auto getter = map.find(hash_pos);
  vector<PointMass *> *neighbors = getter->second;
  Vector3D totalCorrection = Vector3D(0,0,0);
  int num_corrections = 0;

  for (PointMass *neighbor : *neighbors) {
      if (neighbor == &pm) {
          continue;
      }
      Vector3D neighborToPm = (pm.position - neighbor->position);
      float dist = neighborToPm.norm();
      neighborToPm.normalize();
      if (dist < 2 * thickness) {
          Vector3D corrected = neighbor->position + (2 * thickness)*neighborToPm;
          Vector3D correction = corrected - pm.position;
          totalCorrection += correction;
          num_corrections += 1;
      }
  }

  if (num_corrections > 0) {
      totalCorrection /= (num_corrections * simulation_steps);
      pm.position += totalCorrection;
  }
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
  double w = 3 * width / num_width_points;
  double h = 3 * height / num_height_points;
  double t = max(w, h);
  Vector3D truncated = Vector3D((pos.x - fmod(pos.x, w)) / w, (pos.y - fmod(pos.y, h)) / h, (pos.z - fmod(pos.z, t)) / t);
  float hash_value = truncated.x * pow(3, 3) + truncated.y * pow(3, 2) + truncated.z * pow(3, 1);
  return hash_value;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
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
