#version 330

uniform vec4 u_color;

in vec4 v_position;
in vec4 v_normal;
in float out_vertex;
in float out_height;
in float out_z;
in float out_x;

out vec4 out_color;

void main() {
  out_color = (vec4(1, 1, 1, 0) + v_normal) / 2;
  out_color.a = 1;
}
