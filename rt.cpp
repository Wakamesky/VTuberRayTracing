#include <iostream>
#include <fstream>
#include <vector>
#include <optional>
#include <cmath>

struct V {
    double x;
    double y;
    double z;
    V(double v = 0)
        : V(v, v, v) {}
    V(double x, double y, double z)
        : x(x), y(y), z(z) {}
    double operator[](int i) const {
        return (&x)[i];
    }
};

V operator+(V a, V b) {
    return V(a.x + b.x, a.y + b.y, a.z + b.z);
}

V operator-(V a, V b) {
    return V(a.x - b.x, a.y - b.y, a.z - b.z);
}

V operator*(V a, V b) {
    return V(a.x * b.x, a.y * b.y, a.z * b.z);
}

V operator/(V a, V b) {
    return V(a.x / b.x, a.y / b.y, a.z / b.z);
}

V operator-(V v) {
    return V(-v.x, -v.y, -v.z);
}

double dot(V a, V b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

V cross(V a, V b) {
    return V(a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x);
}

V normalize(V v) {
    return v / sqrt(dot(v, v));
}

struct Ray {
    V o;
    V d;
};

struct Sphere;
struct Hit {
    double t;
    const Sphere* sphere;
};

struct Sphere {
    V p;
    double r;
    std::optional<Hit> intersect(const Ray& ray, double tmin, double tmax) const {
        const V op = p - ray.o;
        const double b = dot(op, ray.d);
        const double det = b * b - dot(op, op) + r * r;
        if (det < 0) { return {}; }
        
        const double t1 = b - sqrt(det);
        if (tmin < t1 && t1 < tmax) { return Hit{t1, this}; }
        const double t2 = b + sqrt(det);
        if (tmin < t2 && t2 < tmax) { return Hit{t2, this}; }

        return {};
    }
};

struct Scene {
    std::vector<Sphere> spheres {
        { V(), 1 }
    };
    std::optional<Hit> intersect(const Ray& ray, double tmin, double tmax) const {
        std::optional<Hit> minh;
        for(const auto& sphere : spheres) {
            const auto h = sphere.intersect(ray, tmin, tmax);
            if(!h) { continue; }
            minh = h;
            tmax = minh->t;
        }
        return minh;        
    }
};

int main() {
    const int w = 1200;
    const int h = 800;
    Scene scene;
    std::ofstream ofs("result.ppm");
    ofs << "P3\n" << w << " " << h << "\n255\n";
    for(int i = 0; i < w * h; i++) {
        const int x = i % w;
        const int y = i / w;
        Ray ray;
        ray.o = V(2.*(double)x/(double)w-1.0, 2.*(double)y/(double)h-1.0, 5.0);
        ray.d = V(0, 0, -1);

        const auto h = scene.intersect(ray, 0, 1e+10);
        if (h) {
            ofs << "255 0 255\n";
        } else {
            ofs << "0 0 0\n";
        }
    }
    return 0;
}