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

  double t_pos = dot(point - pm.predict_position, normal) / dot(-normal, normal);
  double t_lastpos = dot(point - pm.position, normal) / dot(-normal, normal);

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
      Vector3D correction = correctionPoint - pm.position;
      pm.predict_position = pm.position + (1 - friction) * correction;
      
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

  MatrixXf positions(3, 4);
  MatrixXf normals(3, 4);
  MatrixXf vert(1, 4);

  positions.col(0) << sPoint + 2 * (sCross + sParallel);
  positions.col(1) << sPoint + 2 * (sCross - sParallel);
  positions.col(2) << sPoint + 2 * (-sCross + sParallel);
  positions.col(3) << sPoint + 2 * (-sCross - sParallel);

  normals.col(0) << sNormal;
  normals.col(1) << sNormal;
  normals.col(2) << sNormal;
  normals.col(3) << sNormal;

  for (int i = 0; i < 4 ; i++) {
    vert.col(i) << 0.0;
  }
  shader.uploadAttrib("is_vertex", vert);

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

  return;









  // nanogui::Color color(1.0f, 1.0f, 0.2f, 1.0f);


  // Vector3f sPoint(point.x, point.y, point.z);
  // Vector3f sNormal(normal.x, normal.y, normal.z);
  // Vector3f sParallel(normal.y - normal.z, normal.z - normal.x,
  //                    normal.x - normal.y);
  // sParallel.normalize();
  // Vector3f sCross = sNormal.cross(sParallel);

  // MatrixXf positions(4, 8);

  positions.col(0) << sPoint , 1.0 ;
  positions.col(1) << sPoint + 6 * (sCross + sParallel), 1.0;
  
  positions.col(2) << sPoint, 1.0  ;
  positions.col(3) << sPoint + 6 * (sCross - sParallel), 1.0 ;
  
  positions.col(4) << sPoint, 1.0 ;
  positions.col(5) << sPoint + 6 * (-sCross + sParallel), 1.0 ;
  
  positions.col(6) << sPoint , 1.0  ;
  positions.col(7) << sPoint + 6 * (-sCross - sParallel), 1.0 ;




  
  

  shader.uploadAttrib("in_position", positions, false);
  // Commented out: the wireframe shader does not have this attribute
  //shader.uploadAttrib("in_normal", normals);

  shader.drawArray(GL_LINES, 0, 8);
  #ifdef LEAK_PATCH_ON
    shader.freeAttrib("in_position");
  #endif

  return;}

  // Draw springs as lines

//   int si = 0;

//   for (int i = 0; i < cloth->springs.size(); i++) {
//     Spring s = cloth->springs[i];

//     if ((s.spring_type == STRUCTURAL && !cp->enable_structural_constraints) ||
//         (s.spring_type == SHEARING && !cp->enable_shearing_constraints) ||
//         (s.spring_type == BENDING && !cp->enable_bending_constraints)) {
//       continue;
//     }

//     Vector3D pa = s.pm_a->position;
//     Vector3D pb = s.pm_b->position;

//     Vector3D na = s.pm_a->normal();
//     Vector3D nb = s.pm_b->normal();

//     positions.col(si) << pa.x, pa.y, pa.z, 1.0;
//     positions.col(si + 1) << pb.x, pb.y, pb.z, 1.0;

//     normals.col(si) << na.x, na.y, na.z, 0.0;
//     normals.col(si + 1) << nb.x, nb.y, nb.z, 0.0;

//     si += 2;
//   }

//   //shader.setUniform("u_color", nanogui::Color(1.0f, 1.0f, 1.0f, 1.0f), false);
//   shader.uploadAttrib("in_position", positions, false);
//   // Commented out: the wireframe shader does not have this attribute
//   //shader.uploadAttrib("in_normal", normals);

//   shader.drawArray(GL_LINES, 0, num_springs * 2);

// #ifdef LEAK_PATCH_ON
//   shader.freeAttrib("in_position");
// #endif



  

  

//   normals.col(0) << sNormal;
//   normals.col(1) << sNormal;
//   normals.col(2) << sNormal;
//   normals.col(3) << sNormal;

//   if (shader.uniform("u_color", false) != -1) {
//     shader.setUniform("u_color", color);
//   }
//   shader.uploadAttrib("in_position", positions);
//   if (shader.attrib("in_normal", false) != -1) {
//     shader.uploadAttrib("in_normal", normals);
//   }

//   shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);
// #ifdef LEAK_PATCH_ON
//   shader.freeAttrib("in_position");
//   if (shader.attrib("in_normal", false) != -1) {
//     shader.freeAttrib("in_normal");
//   }
// #endif
// }
