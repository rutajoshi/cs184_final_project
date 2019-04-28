#ifndef POINTMASS_H
#define POINTMASS_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

using namespace CGL;

// Forward declarations
class Halfedge;

struct PointMass {
  PointMass(Vector3D position, bool pinned)
      : pinned(pinned), start_position(position), position(position), last_position(position),
        predict_position(position), velocity(Vector3D()), delta_position(Vector3D()) {}

  Vector3D normal();

  // static values
  bool pinned;
  Vector3D start_position;
  double rest_density = 1.0;

  // dynamic values
  Vector3D position;
  Vector3D last_position;
  Vector3D predict_position;
  Vector3D velocity;
  Vector3D forces;
  Vector3D delta_position;
  double lambda = 0.0;

  // mesh reference
  Halfedge *halfedge;
};

#endif /* POINTMASS_H */
