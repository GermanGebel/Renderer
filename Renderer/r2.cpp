#define _USE_MATH_DEFINES

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include "geometry.h"

using namespace std;

const int width = 900;
const int height = 600;

vector<Vec3f> framebuffer(width* height);

struct Vector3 {
    float x, y, z;
    Vector3(float x, float y, float z) : x(x), y(y), z(z) {};
    Vector3() : x(0), y(0), z(0) {};

    float dot(Vector3 v) { return x * v.x + y * v.y + z * v.z; }
    
    Vector3 operator-(Vector3 v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    Vector3 operator+(Vector3 v) {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    Vector3 operator*(float a) {
        x *= a;
        y *= a;
        z *= a;
    }

};

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

    bool isIntersect(const Vec3f& orig, const Vec3f& dir, float& t0) const {
        Vec3f L = center - orig;
        float tca = L * dir;
        float d2 = L * L - tca * tca;
        if (d2 > radius * radius) return false;
        float thc = sqrtf(radius * radius - d2);
        t0 = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
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

    bool isIntersect(const Vec3f& orig, const Vec3f& dir, float& t, Vec3f& hit)
    {
        if (abs(normal.x) == 1.f){
            t = center.x < 0 ? (orig.x + center.x) / dir.x : (orig.x + center.x) / dir.x;
            hit = orig + dir * t;

            if (t > 0 && hit.y < center.y + width / 2 && hit.y > center.y - width / 2 && hit.z < center.z + height / 2 && hit.z > center.z - height / 2) {
                return true;
            }
        }
        if (abs(normal.y) == 1.f) {

            t = center.y < 0 ? (orig.y + center.y) / dir.y  : (orig.y + center.y) / dir.y;
            hit = orig + dir * t;

            if (t > 0 && hit.x < center.x + width / 2 && hit.x > center.x - width / 2 && hit.z < center.z + height / 2 && hit.z > center.z - height / 2) {
                return true;
            }
        }
        if (abs(normal.z) == 1.f) {
            t = center.z < 0 ? (orig.z + center.z) / dir.z : (orig.z + center.z) / dir.z;
            hit = orig + dir * t;

            if (t > 0 && hit.x < center.x + width / 2 && hit.x > center.x - width / 2 && hit.y < center.y + height / 2 && hit.y > center.y - height / 2) {
                return true;
            }
        }
        return false;                
    }
};
//
//struct Box {
//    Vec3f boxMin, boxMax;
//    Box(Vec3f bmin, Vec3f bmax): boxMin(bmin), boxMax(bmax) {}
//
//    bool isIntersect(Vec3f orig, Vec3f dir) {
//        auto x0 = min((boxMin.x - orig.x) / dir.x, (boxMax.x - orig.x) / dir.x);
//        auto x1 = max((boxMin.x - orig.x) / dir.x, (boxMax.x - orig.x) / dir.x);
//
//        auto y0 = min((boxMin.y - orig.y) / dir.y, (boxMax.y - orig.y) / dir.y);
//        auto y1 = max((boxMin.y - orig.y) / dir.y, (boxMax.y - orig.y) / dir.y);
//
//        auto z0 = min((boxMin.z - orig.z) / dir.z, (boxMax.z - orig.z) / dir.z);
//        auto z1 = max((boxMin.z - orig.z) / dir.z, (boxMax.z - orig.z) / dir.z);
//
//        float t_min = max(min(x0, x1), min(y0, y1), min(z0, z1));
//        float t_max = min(max(x0, x1), max(y0, y1), max(z0, z1));
//
//        if (t_max <= t_min) return false;
//        return true;
//    }
//
//    
//};

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

bool scene_intersect(const Vec3f& orig, const Vec3f& dir, Vec3f& hit, Vec3f& N, Material& material) {
    float spheres_dist = numeric_limits<float>::max();

    for (size_t i = 0; i < spheres.size(); i++) {
        float t;
        if (spheres[i].isIntersect(orig, dir, t) && t < spheres_dist) {
            spheres_dist = t;
            hit = orig + dir * t;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    float rect_dist = numeric_limits<float>::max();
    
    for (int i = 0; i < rectangles.size(); i++) {
        float t;
        if (rectangles[i].isIntersect(orig, dir, t, hit) && t < spheres_dist) {
            rect_dist = t;
            N = rectangles[i].normal;
            material = rectangles[i].material;
        }
    }
    
    return min(spheres_dist, rect_dist) < 1000;
}


Vec3f trace_ray(const Vec3f& orig, const Vec3f& dir, size_t depth = 0) {
    Vec3f point, N;
    Material material;

    if (depth > 4 || !scene_intersect(orig, dir, point, N, material)) {
        return Vec3f(0, 0, 0); // background color
    }
    
    Vec3f reflect_color = Vec3f(0, 0, 0);
    Vec3f refract_color = Vec3f(0, 0, 0);

    if (material.specular_exponent != 0) {
        Vec3f reflect_dir = reflect(dir, N).normalize();
        Vec3f reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // offset the original point to avoid occlusion by the object itself
        reflect_color = trace_ray(reflect_orig, reflect_dir, depth + 1);
    }

    if (material.k_refraction != 0) {
        Vec3f refract_dir = refract(dir, N, material.k_refraction).normalize();
        Vec3f refract_orig = refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
        refract_color = trace_ray(refract_orig, refract_dir, depth + 1);
    }

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        // shadows
        Vec3f shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmp_material;
        if (scene_intersect(shadow_orig, light_dir, shadow_pt, shadow_N, tmp_material) && (shadow_pt - shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity += lights[i].intensity * max(0.f, light_dir * N);
        specular_light_intensity += powf(max(0.f, -reflect(-light_dir, N) * dir), material.specular_exponent) * lights[i].intensity;
    }
    return (material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1] + reflect_color * material.albedo[2] + refract_color * material.albedo[3]);
}

void saveImage() {
    ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm",ios::binary);
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

void initScene() {
    Material ivory(1.0, Vec4f(0.6, 0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3), 50.);
    Material glass(1.5, Vec4f(0.0, 0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8), 125.);
    Material red_rubber(1.0, Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1), 10.);
    Material mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);

    Material red_wall(0, Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.9f, 0.0f, 0.0f), 0);
    Material green_wall(0, Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.f, 0.9f, 0.f), 0);
    Material white_wall(0, Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.9f, 0.8f, 0.7f), 0);

    /*spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));*/
    //spheres.push_back(Sphere(Vec3f(0, 0, -10), 2, mirror));
    spheres.push_back(Sphere(Vec3f(0, 0, -10), 1, red_rubber));
    /*spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, glass));*/

    // bottom
    rectangles.push_back(Rectangle(Vec3f(0, -7, -10), Vec3f(0, 1, 0), 20, 20, white_wall));
    // top 
    rectangles.push_back(Rectangle(Vec3f(0, 7, -10), Vec3f(0, -1, 0), 20, 20, white_wall));
    // back
    rectangles.push_back(Rectangle(Vec3f(0, 0, -20), Vec3f(0, 0, 1), 20, 14, white_wall));
    // left
    rectangles.push_back(Rectangle(Vec3f(10, 0, -10), Vec3f(1, 0, 0), 20, 14, red_wall));
    // right
    rectangles.push_back(Rectangle(Vec3f(10, 0, -10), Vec3f(-1, 0, 0), 20, 14, green_wall));


    lights.push_back(Light(Vec3f(0, 5, -15), 2));

    //lights.push_back(Light(Vec3f(0, 5, -10), 2.1));
    /*lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f(30, 20, 30), 1.7));*/
};

void render() {    
    const float fov = M_PI / 3.;
    
    int pixels_amount = width * height;
    int count = 0;


#pragma omp parallel for
    for (size_t j = 0; j < height; j++) { // actual rendering loop
        for (size_t i = 0; i < width; i++) {

            float dir_x = (i + 0.5) - width / 2.;
            float dir_y = -(j + 0.5) + height / 2.;    // this flips the image at the same time
            float dir_z = -height / (2. * tan(fov / 2.));

            framebuffer[i + j * width] = trace_ray(Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z).normalize());
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
    saveImage();

    return 0;
}