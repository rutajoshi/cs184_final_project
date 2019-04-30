#version 330

uniform vec4 u_color;

in vec4 v_position;
in vec4 v_normal;
in float out_vertex;
in float out_height;

out vec4 out_color;

void main() {
  vec2 pt = gl_PointCoord - vec2(0.5);
  if(out_vertex == 1.0 && pt.x*pt.x+pt.y*pt.y > 0.25)
      discard;
  if(out_vertex == 1.0) 
  	out_color = vec4(0, 0, out_height, 1);
  if(out_vertex != 1.0) 
  	out_color = vec4(.1,0,0,.001);
}
