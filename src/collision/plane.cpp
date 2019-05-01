#include "iostream"
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "../leak_fix.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

#define BOUNCE_DAMPING_FACTOR 0.1

void Plane::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with planes.
  Vector3D inter_position = pm.position;

  double t_pos = dot(point - pm.predict_position, normal) / dot(-normal, normal);
  double t_lastpos = dot(point - inter_position, normal) / dot(-normal, normal);

  // then the points are on opposite sides of the plane
  if (t_pos * t_lastpos <= 0) {
//      std::cout<<"\n collision occured \n";
      Vector3D correctionPoint;
      if (t_pos <= 0) {
        Vector3D tangentPoint = pm.predict_position + t_pos * (-normal);
        correctionPoint = tangentPoint + SURFACE_OFFSET * normal;
      } else {
        Vector3D tangentPoint = pm.predict_position + t_pos * (-normal);
        correctionPoint = tangentPoint - SURFACE_OFFSET * normal;
      }
      Vector3D correction = correctionPoint - inter_position;
      pm.predict_position = inter_position + (1 - friction) * correction;

//      assert(pm.predict_position.y <= 1.4);
//      assert(pm.predict_position.y >= -0.2);

      double new_t_pos = dot(point - pm.predict_position, normal) / dot(-normal, normal);
      double new_t_lastpos = dot(point - inter_position, normal) / dot(-normal, normal);
      assert(new_t_pos * new_t_lastpos >= 0);

      // If collision --> apply normal force in the opposite direction
//      pm.num_collisions += 1;
//      double cos_theta = dot(-pm.velocity, normal) / pm.velocity.norm();
//      double theta = acos(cos_theta);
//      double delta_t = 1.0f / 90; // 0.016; //
//      Vector3D normal_force = (-2 * pm.mass * pm.velocity * cos_theta / delta_t);
//      pm.collision_forces += normal_force;
  }
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
