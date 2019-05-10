#include "iostream"
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "../leak_fix.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.01

#define BOUNCE_DAMPING_FACTOR 0.1

int sign(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

void Plane::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with planes.

  // find out if points are on the same side
  int position_side = sign(dot(normal, pm.position - point));
  int predict_position_side = sign(dot(normal, pm.predict_position - point));

  // if they are not the same, apply correction
  if (position_side != predict_position_side) {
      // direction of current point mass
      Vector3D motion_direction = pm.predict_position - pm.position;

      motion_direction.normalize();

      // find t value (intersection of direction with plane)
      double t = (dot(point - pm.position, normal)) / (dot(motion_direction, normal));

      // get tangent point
      Vector3D tangent_point = pm.position + t * motion_direction;

      // offset from surface of plane
      Vector3D correction_point;
      Vector3D correction_vector;
      if (position_side >= 0) {
          // position is on same side of normal
          correction_point = tangent_point + normal * SURFACE_OFFSET;
      } else {
          // position in opposite direction of normal
          correction_point = tangent_point - normal * SURFACE_OFFSET;
      }
      correction_vector = correction_point - pm.position;

      // apply correction to position
      pm.predict_position = pm.position + (1.0 - friction) * correction_vector;
  }

  int corrected_predict_position_side = sign(dot(normal, pm.predict_position - point));

  // check that correction moves the point to the correct side
  assert(corrected_predict_position_side == position_side);

  // check for bottom plane
//  assert(pm.predict_position.y > -0.15);

//  }
}

void Plane::render(GLShader &shader) {
  nanogui::Color color(0.7f, 0.7f, 0.7f, 1.0f);

  Vector3f sPoint(point.x, point.y, point.z);
  Vector3f sNormal(normal.x, normal.y, normal.z);
  Vector3f sParallel(normal.y - normal.z, normal.z - normal.x,
                     normal.x - normal.y);
  sParallel.normalize();
  Vector3f sCross = sNormal.cross(sParallel);
  sCross.normalize();

  MatrixXf positions(3, 4);
  MatrixXf normals(3, 4);
  MatrixXf vert(1, 4);
  MatrixXf h(1, 4);

  Vector3f top = sPoint + 1 * (sCross);
  Vector3f bot = sPoint + 1 * (- sParallel);
  Vector3f lef = sPoint + 1 * (sParallel);
  Vector3f rig = sPoint + 1 * (-sCross );

  // top.normalize();
  // bot.normalize();
  // lef.normalize();
  // rig.normalize();



  positions.col(0) << top;
  positions.col(1) << bot;
  positions.col(2) << lef;
  positions.col(3) << rig;

  normals.col(0) << sNormal;
  normals.col(1) << sNormal;
  normals.col(2) << sNormal;
  normals.col(3) << sNormal;

  for (int i = 0; i < 4 ; i++) {
    h.col(i) << -1;
    vert.col(i) << 0.0;
  }
  shader.uploadAttrib("is_vertex", vert);
  shader.uploadAttrib("height", h);
  shader.uploadAttrib("xpos", h);
  shader.uploadAttrib("zpos", h);

  if (shader.uniform("u_color", false) != -1) {
    shader.setUniform("u_color", color);
  }
  shader.uploadAttrib("in_position", positions);
  if (shader.attrib("in_normal", false) != -1) {
    shader.uploadAttrib("in_normal", normals);
  }

  shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);
#ifdef LEAK_PATCH_ON
  shader.freeAttrib("in_position");
  shader.freeAttrib("is_vertex");
  if (shader.attrib("in_normal", false) != -1) {
    shader.freeAttrib("in_normal");
  }
#endif
}
