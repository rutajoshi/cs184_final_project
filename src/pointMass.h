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
      : pinned(pinned), start_position(position), position(position),
        last_position(position), velocity(Vector3D()) {}

  Vector3D normal();

  // static values
  bool pinned;
  Vector3D start_position;

  // dynamic values
  Vector3D position;
  Vector3D predict_position;
  Vector3D velocity;
  Vector3D forces;

  // mesh reference
  Halfedge *halfedge;
};

#endif /* POINTMASS_H */
