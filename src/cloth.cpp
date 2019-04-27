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
  // TODO (Part 1): Build a grid of masses and springs.

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
            if (std::find(pinned.begin(), pinned.end(), rc) != pinned.end()) {
                pinnedPoint = true;
            }
            point_masses.push_back(PointMass(pos, pinnedPoint));
        }
    }

    // Make all the springs
    // for (int r = 0; r < num_height_points; r++) {
    //     for (int c = 0; c < num_width_points; c++) {
    //
    //         // For each point, (r,c), generate the 6 possible springs
    //         if (c > 0) {
    //             // Make left structural constraint spring
    //             // PointMass leftStruc = point_masses[r * num_width_points + (c - 1)];
    //             springs.push_back(Spring(&point_masses[r * num_width_points + c],
    //                                      &point_masses[r * num_width_points + (c - 1)],
    //                                      STRUCTURAL));
    //         }
    //
    //         if (r > 0) {
    //             // Make above structural constraint spring
    //             springs.push_back(Spring(&point_masses[r * num_width_points + c],
    //                                      &point_masses[(r - 1) * num_width_points + c],
    //                                      STRUCTURAL));
    //         }
    //
    //         if (r > 0 && c > 0) {
    //             // Make the diagonal upper left shearing constraint spring
    //             springs.push_back(Spring(&point_masses[r * num_width_points + c],
    //                                      &point_masses[(r - 1) * num_width_points + (c - 1)],
    //                                      SHEARING));
    //         }
    //
    //         if (r > 0 && c < num_width_points - 1) {
    //             // Make the diagonal right shearing constraint spring
    //             springs.push_back(Spring(&point_masses[r * num_width_points + c],
    //                                      &point_masses[(r - 1) * num_width_points + (c + 1)],
    //                                      SHEARING));
    //         }
    //
    //         if (c > 1) {
    //             // Make the 2left bending constraint spring
    //             springs.push_back(Spring(&point_masses[r * num_width_points + c],
    //                                      &point_masses[r * num_width_points + (c - 2)],
    //                                      BENDING));
    //         }
    //
    //         if (r > 1) {
    //             // Make the 2above bending constraint spring
    //             springs.push_back(Spring(&point_masses[r * num_width_points + c],
    //                                      &point_masses[(r - 2) * num_width_points + c],
    //                                      BENDING));
    //         }
    //     }
    // }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.
  Vector3D total_ext_force = Vector3D(0,0,0);
  for (int i = 0; i < external_accelerations.size(); i++) {
      total_ext_force += external_accelerations[i] * mass;
  }

  for (PointMass &pm : point_masses) {
      pm.forces = total_ext_force;
  }

  // for (Spring &sp : springs) {
  //     double Fs = 0;
  //     if ((sp.spring_type == STRUCTURAL && cp->enable_structural_constraints) ||
  //         (sp.spring_type == SHEARING && cp->enable_shearing_constraints)) {
  //         Fs = cp->ks * ((sp.pm_a->position - sp.pm_b->position).norm() - sp.rest_length);
  //     } else if (sp.spring_type == BENDING && cp->enable_bending_constraints) {
  //         Fs = 0.2 * cp->ks * ((sp.pm_a->position - sp.pm_b->position).norm() - sp.rest_length);
  //     } else {
  //         continue;
  //     }
  //     Vector3D direction = (sp.pm_b->position - sp.pm_a->position);
  //     direction.normalize();
  //     Vector3D spring_force = Fs * direction;
  //     sp.pm_a->forces += spring_force;
  //     sp.pm_b->forces -= spring_force;
  // }

  // TODO (Part 2): Use Verlet integration to compute new point mass positions
  for (PointMass &pm : point_masses) {
      if (!pm.pinned) {
          Vector3D xt = Vector3D(pm.position);
          Vector3D old_last = Vector3D(pm.last_position);
          pm.last_position = xt;
          Vector3D acc = pm.forces / mass;
          pm.position = xt + (1 - cp->damping/100.0)*(xt - old_last) + acc * pow(delta_t, 2);
      }
  }

  build_spatial_map();

  // TODO (Part 4): Handle self-collisions.
  for (PointMass &pm : point_masses) {
      self_collide(pm, simulation_steps);
  }

  // TODO (Part 3): Handle collisions with other primitives.
  for (PointMass &pm : point_masses) {
      for (CollisionObject* co : *collision_objects) {
          co->collide(pm);
      }
  }

  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  // for (Spring &sp : springs) {
  //     double dist = (sp.pm_a->position - sp.pm_b->position).norm();
  //     if (dist > 1.1 * sp.rest_length) {
  //         // If both are pinned, move neither
  //         if (sp.pm_a->pinned && sp.pm_b->pinned) {
  //             continue;
  //         }
  //         // If one of the two is pinned, only move the not-pinned one
  //         else if (sp.pm_a->pinned) {
  //             Vector3D directionAtoB = sp.pm_b->position - sp.pm_a->position;
  //             directionAtoB.normalize();
  //             sp.pm_b->position = sp.pm_a->position + 1.1 * sp.rest_length * directionAtoB;
  //         } else if (sp.pm_b->pinned) {
  //             Vector3D directionBtoA = sp.pm_a->position - sp.pm_b->position;
  //             directionBtoA.normalize();
  //             sp.pm_a->position = sp.pm_b->position + 1.1 * sp.rest_length * directionBtoA;
  //         }
  //         // If neither is pinned, move both
  //         else {
  //             Vector3D directionAtoB = sp.pm_b->position - sp.pm_a->position;
  //             directionAtoB.normalize();
  //             Vector3D halfway = (sp.pm_a->position + sp.pm_b->position) / 2.0;
  //             sp.pm_a->position = halfway - 0.55 * sp.rest_length * directionAtoB;
  //             sp.pm_b->position = halfway + 0.55 * sp.rest_length * directionAtoB;
  //         }
  //     }
  // }

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
