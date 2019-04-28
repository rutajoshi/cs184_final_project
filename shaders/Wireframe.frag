#version 330

uniform vec4 u_color;

in vec4 v_position;
in vec4 v_normal;

out vec4 out_color;

void main() {
  vec2 pt = gl_PointCoord - vec2(0.5);
  if(pt.x*pt.x+pt.y*pt.y > 0.25)
      discard;
  out_color = u_color;
}
