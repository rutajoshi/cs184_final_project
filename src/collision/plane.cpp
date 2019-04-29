#include "iostream"
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "../leak_fix.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with planes.
  Vector3D inter_position = pm.predict_position - pm.delta_position;

  double t_pos = dot(point - pm.predict_position, normal) / dot(-normal, normal);
  double t_lastpos = dot(point - inter_position, normal) / dot(-normal, normal);

  // then the points are on opposite sides of the plane
  if (t_pos * t_lastpos <= 0) {
      Vector3D correctionPoint;
      if (t_pos < 0) {
        Vector3D tangentPoint = pm.predict_position + t_pos * (-normal);
        correctionPoint = tangentPoint + SURFACE_OFFSET * normal;
      } else {
        Vector3D tangentPoint = pm.predict_position + t_pos * (normal);
        correctionPoint = tangentPoint - SURFACE_OFFSET * normal;
      }
      Vector3D correction = correctionPoint - inter_position;
      pm.predict_position = inter_position + (1 - friction) * correction;
      
  }
}

void Plane::render(GLShader &shader) {
  return;
  nanogui::Color color(0.7f, 0.7f, 0.7f, 1.0f);

  Vector3f sPoint(point.x, point.y, point.z);
  Vector3f sNormal(normal.x, normal.y, normal.z);
  Vector3f sParallel(normal.y - normal.z, normal.z - normal.x,
                     normal.x - normal.y);
  sParallel.normalize();
  Vector3f sCross = sNormal.cross(sParallel);

  MatrixXf positions(3, 4);
  MatrixXf normals(3, 4);

  positions.col(0) << sPoint + 2 * (sCross + sParallel);
  positions.col(1) << sPoint + 2 * (sCross - sParallel);
  positions.col(2) << sPoint + 2 * (-sCross + sParallel);
  positions.col(3) << sPoint + 2 * (-sCross - sParallel);

  normals.col(0) << sNormal;
  normals.col(1) << sNormal;
  normals.col(2) << sNormal;
  normals.col(3) << sNormal;

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
  if (shader.attrib("in_normal", false) != -1) {
    shader.freeAttrib("in_normal");
  }
#endif
}
