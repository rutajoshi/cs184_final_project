#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  float ka = 0.1;
  float kd = 1.0;
  float ks = 0.5;
  vec3 Ia = vec3(1.0, 1.0, 1.0);
  float p = 32.0;

  vec3 l = u_light_pos - v_position.xyz;
  float radius = length(l);

  float angle = dot(v_normal.xyz, normalize(l));

  vec3 h = (normalize(u_cam_pos - v_position.xyz) + normalize(u_light_pos - v_position.xyz)) / 2.0;

  float angleh = dot(v_normal.xyz, normalize(h));

  out_color.xyz = ka * Ia +
                  kd * (u_light_intensity / pow(radius, 2.0)) * max(0.0, angle) +
                  ks * (u_light_intensity / pow(radius, 2.0)) * pow(max(0.0, angleh), p);
  out_color.a = 1.0;
}
