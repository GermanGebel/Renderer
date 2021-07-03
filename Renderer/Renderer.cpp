#define _USE_MATH_DEFINES
#define OMP_NUM_THREADS 1000

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include "geometry.h"

using namespace std;

const int width = 512;
const int height = 384;

vector<Vec3f> framebuffer(width* height);

struct Light {
    Light(const Vec3f& p, const float i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Material {
    Material(const float r, const Vec4f& a, const Vec3f& color, const float spec) : k_refraction(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : k_refraction(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}
    Vec4f albedo;
    //float k_diffusion;
    float specular_exponent;
    float k_refraction;
    Vec3f diffuse_color;
    //Vec3f color;    
};

struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f& c, const float r, const Material& m) : center(c), radius(r), material(m) {}

    bool isIntersect(const Vec3f& origin, const Vec3f& dir, float& t) const {
        Vec3f co = center - origin;
        co = origin - center;
        float a = dir * dir;
        float halb_b = (co * dir);
        float c = (co * co) - (radius * radius);
        float discriminant = halb_b * halb_b - a * c;
        if (discriminant < 0) return false;
        float t1 = (-halb_b - sqrt(discriminant)) / a;
        float t2 = (-halb_b + sqrt(discriminant)) / a;
        t = t1;
        if (t < 0) t = t2;
        if (t < 0) return false;
        return true;
    }
};

struct Rectangle {
    Vec3f center;
    Vec3f normal;
    float width;
    float height;
    Material material;

    Rectangle(Vec3f c, Vec3f n, float w, float h, Material m) :center(c), normal(n), width(w), height(h), material(m) {}

    bool isIntersect(Vec3f& origin, Vec3f& dir, float& t)
    {
        //if (fabs(dir.x) > 1e-3) t = - (origin.x + 10) / dir.x;
        //if (fabs(dir.y) > 1e-3) t = - (origin.y + 4) / dir.y;
        //if (fabs(dir.z) > 1e-3) t = -(origin.z + 30) / dir.z;
        ///*if (fabs(dir.y) > 1e-3) {
        //    float t = -(origin.y + 4) / dir.y;
        //    Vec3f xt = origin + dir * t;
        //    if (t > 0 && fabs(xt.x) < 10 && xt.z < 0 && xt.z>-30 && t < min_distance) {
        //        min_distance = t;
        //        hit = xt;
        //        N = Vec3f(0.f, 1.f, 0.f);
        //        material = white_wall;
        //    }
        //}*/
    }
};

vector<Sphere> spheres;
vector<Rectangle> rectangles;
vector<Light>  lights;

Vec3f reflect(const Vec3f& I, const Vec3f& N) {
    return I - N * 2.f * (I * N);
}

Vec3f refract(const Vec3f& I, const Vec3f& N, const float eta_t, const float eta_i = 1.f) {
    float cosi = -max(-1.f, min(1.f, I * N));
    if (cosi < 0) return refract(I, -N, eta_i, eta_t);
    float eta = eta_i / eta_t;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(1, 0, 0) : I * eta + N * (eta * cosi - sqrtf(k));
}

bool scene_intersect(const Vec3f& origin, const Vec3f& dir, Vec3f& hit, Vec3f& N, Material& material) {
    float min_distance = numeric_limits<float>::max();
    for (size_t i = 0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].isIntersect(origin, dir, dist_i) && dist_i < min_distance) {
            min_distance = dist_i;
            hit = origin + dir * dist_i;
            N = (spheres[i].center - hit).normalize();
            material = spheres[i].material;
        }
    }

    // walls
    Material red_wall = Material(0, Vec4f(1, 0, 0.0, 0.0), Vec3f(0.849, 0.260, 0.270), 0);
    Material green_wall = Material(0, Vec4f(1, 0, 0.0, 0.0), Vec3f(0.260, 0.849, 0.270), 0);
    Material white_wall = Material(0, Vec4f(1, 0, 0.0, 0.0), Vec3f(0.507, 0.507, 0.513), 0);
    Material mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);

    /*for (size_t i = 0; i < rectangles.size(); i++) {
        float dist_i;
        if (rectangles[i].isIntersect(origin, dir, dist_i) && dist_i < min_distance) {
            min_distance = dist_i;
            hit = origin + dir * dist_i;
            N = (spheres[i].center - hit).normalize();
            material = spheres[i].material;
        }
    }*/

    //bottom
    if (fabs(dir.y) > 1e-3) {
        float d = -(origin.y + 4) / dir.y;
        Vec3f pt = origin + dir * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.z < 0 && pt.z>-30 && d < min_distance) {
            min_distance = d;
            hit = pt;
            N = Vec3f(0.f, 1.f, 0.f);
            material = mirror;    
        }
    }
    // back
    if (fabs(dir.z) > 1e-3) {
        float d = -(origin.z + 30) / dir.z; 
        Vec3f pt = origin + dir * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.y > -4 && pt.y < 10 && d < min_distance) {
            min_distance = d;
            hit = pt;
            N = Vec3f(0.f, 0.f, 1.f);
            material = white_wall;    
        }
    }
    
    // left
    if (fabs(dir.x) > 1e-3) {
        float d = -(origin.x + 10) / dir.x; 
        Vec3f pt = origin + dir * d;
        if (d > 0 && pt.y > -4 && pt.y < 10 && pt.z < 0 && pt.z>-30 && d < min_distance) {
            min_distance = d;
            hit = pt;
            N = Vec3f(1.f, 0.f, 0.f);
            material = red_wall;   
        }
    }

    // top
    if (fabs(dir.y) > 1e-3) {
        float d = -(origin.y - 10) / dir.y;
        Vec3f xt = origin + dir * d;
        if (d > 0 && fabs(xt.x) < 10 && xt.z < 0 && xt.z>-30 && d < min_distance) {
            min_distance = d;
            hit = xt;
            N = Vec3f(0.f, -1.f, 0.f);
            material = white_wall;
        }
    }
    // right
    if (fabs(dir.x) > 1e-3) {
        float d = -(origin.x - 10) / dir.x; 
        Vec3f xt = origin + dir * d;
        if (d > 0 && xt.y > -4 && xt.y < 10 && xt.z < 0 && xt.z>-30 && d < min_distance) {
            min_distance = d;
            hit = xt;
            N = Vec3f(-1.f, 0.f, 0.f);
            material = green_wall;   
        }
    }

    return min_distance < 1000;
}

Vec3f trace_ray(const Vec3f& orig, const Vec3f& dir, size_t depth = 0) {
    Vec3f point, N;
    Material material;

    if (depth > 4 || !scene_intersect(orig, dir, point, N, material)) {
        return Vec3f(0, 0, 0); // background color
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f refract_dir = refract(dir, N, material.k_refraction).normalize();
    Vec3f reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // offset the original point to avoid occlusion by the object itself
    Vec3f refract_orig = refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    Vec3f reflect_color = trace_ray(reflect_orig, reflect_dir, depth + 1);
    Vec3f refract_color = trace_ray(refract_orig, refract_dir, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        // shadows
        Vec3f shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // checking if the point lies in the shadow of the lights[i] offset the original point to avoid occlusion by the object itself
        Vec3f shadow_pt, shadow_N;
        Material tmp_material;
        if (scene_intersect(shadow_orig, light_dir, shadow_pt, shadow_N, tmp_material) && (shadow_pt - shadow_orig).norm() < light_distance) continue;

        diffuse_light_intensity += lights[i].intensity * abs(light_dir * N) / (light_distance * light_distance);
        specular_light_intensity += powf(abs(-reflect(-light_dir, N) * dir), material.specular_exponent) * lights[i].intensity / (light_distance * light_distance);
    }
    return (material.diffuse_color * diffuse_light_intensity * material.albedo[0]
        + Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1]
        + reflect_color * material.albedo[2]
        + refract_color * material.albedo[3]) * (1 / M_PI);
}


void initScene() {
    Material ivory(1.0, Vec4f(0.6, 0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3), 50.);
    Material glass(2.5, Vec4f(0.8, 0.5, 0.1, 1.0), Vec3f(1, 1, 1), 12.);
    Material red_rubber = { 1.0, {1, 0, 0.0, 0.0}, {0.3, 0.1, 0.1},   10. };
    Material mirror = { 1.0, {0.0, 10.0, 0.8, 0.0}, {1.0, 1.0, 1.0}, 1425. };
    Material white_wall = Material(0, Vec4f(1, 0, 0.0, 0.0), Vec3f(0.507, 0.507, 0.513), 0);

    spheres.push_back(Sphere(Vec3f(7, 7, -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-5, 0, -10), 2, mirror));
    spheres.push_back(Sphere(Vec3f(-3, 1, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f(6, -2, -13), 2, glass));


    lights.push_back(Light(Vec3f(0, 9.9, -15), 1000));
    lights.push_back(Light(Vec3f(0, 2, 30), 2000));
    //lights.push_back(Light(Vec3f(30, 20, 30), 1.7));
};

void saveImageAsPPM() {
    ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm", ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i) {
        Vec3f& c = framebuffer[i];
        float maximum = max(c[0], max(c[1], c[2]));
        if (maximum > 1) c = c * (1. / maximum);
        for (size_t j = 0; j < 3; j++) {
            ofs << (char)(255 * max(0.f, min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
};
void saveAsTXT() {
    ofstream ofs; // save the framebuffer to file
    ofs.open("./out.txt", ios::binary);
    for (size_t c = 0; c < 3; c++) {
        switch (c)
        {
        case 0:
            ofs << "RED\n";
            break;
        case 1:
            ofs << "GREEN\n";
            break;
        case 2:
            ofs << "BLUE\n";
            break;
        default:
            break;
        }
        for (size_t i = 0; i < height * width; ++i) {
            ofs << (255 * max(0.f, min(1.f, framebuffer[i][c]))) << ' ';
            if ((i + 1) % width == 0) {
                ofs << '\n';
            }
        }
        ofs << '\n';
    }

    ofs.close();
}


void render() {
    const float fov = M_PI / 6;

    int pixels_amount = width * height;
    int count = 0;
    Vec3f camera_position = Vec3f(0, 2, 20);



#pragma omp parallel for
    for (size_t j = 0; j < height; j++) { // actual rendering loop
        for (size_t i = 0; i < width; i++) {
            float dir_x = (i + 0.5) - width / 2.;
            float dir_y = - (j + 0.5) + height / 2.;    // this flips the image at the same time
            float dir_z = - height / (2. * tan(fov / 2.));
            framebuffer[i + j * width] = trace_ray(camera_position, (Vec3f(dir_x, dir_y, dir_z)).normalize());
            //cout << count++ << "/" << pixels_amount << endl;
        }
    }
}

int main() {
    cout << "Initializing scene..." << endl;
    initScene();

    cout << "Rendering..." << endl;
    auto begin = std::chrono::system_clock::now();
    render();
    auto end = chrono::system_clock::now();
    double render_time = (double)chrono::duration_cast<chrono::seconds>(end - begin).count();
    cout << "Rendering time: " << render_time << "s" << endl;

    cout << "Saving image..." << endl;
    saveImageAsPPM();
    saveAsTXT();

    return 0;
}