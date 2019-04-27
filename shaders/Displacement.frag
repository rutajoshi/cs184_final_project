#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_3;
uniform vec2 u_texture_3_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
  return texture(u_texture_3, uv).r;
}

void main() {
  // YOUR CODE HERE
  vec3 b = cross(v_normal.xyz, v_tangent.xyz);
  mat3 tbn = mat3(v_tangent.xyz, b, v_normal.xyz);
  float dU = (h(vec2(v_uv.x + 1.0 / u_texture_3_size.x, v_uv.y)) - h(v_uv)) * u_height_scaling * u_normal_scaling;
  float dV = (h(vec2(v_uv.x, v_uv.y + 1.0 / u_texture_3_size.y)) - h(v_uv)) * u_height_scaling * u_normal_scaling;
  vec3 no = vec3(-dU, -dV, 1.0);
  vec3 nd = tbn * no;

  float ka = 0.1;
  float kd = 1.0;
  float ks = 0.5;
  vec3 Ia = vec3(1.0, 1.0, 1.0);
  float p = 32.0;

  float radius = length(u_light_pos - v_position.xyz);
  float angle = dot(normalize(u_light_pos - v_position.xyz), nd);
  vec3 h = ((u_cam_pos - v_position.xyz) + (u_light_pos - v_position.xyz)) / 2.0;
  float angleh = dot(nd, normalize(h));

  out_color.xyz = ka * Ia +
                  kd * (u_light_intensity / pow(radius, 2.0)) * max(0.0, angle) +
                  ks * (u_light_intensity / pow(radius, 2.0)) * pow(max(0.0, angleh), p);
  out_color.a = 1.0;
}
