#include <iostream>
#include <fstream>
#include <vector>
// 9:00
// 値のほかに無効の情報を持てる
#include <optional>
// sqrtが使えないので...
#include <cmath>
// std::mt19937 のため
#include <random>
#include <algorithm>
#include <tuple>
#include <omp.h>

// 03:48
// vector構造体
// 位置, 方向, 色を表す
// OK
struct V {
    double x;
    double y;
    double z;
    // コンストラクタの定義
    // V() -> (0, 0, 0)
    // v = 0 はデフォルト引数
    // V(v) -> (v, v, v)
    V(double v = 0)
        : V(v, v, v) {}
    // V(x, y, z) -> (x, y, z)
    V(double x, double y, double z)
        : x(x), y(y), z(z) {}
    // Vへの作用素
    double operator[](int i) const {
        return (&x)[i];
    }
};

// 04:32
// ベクトル四則演算の定義
// OK
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

// ベクトル内積 (戻り値はスカラ量なのでdouble)
double dot(V a, V b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// ベクトル外積
V cross(V a, V b) {
    return V(a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x);
}

// 正規化　[0, 1] の範囲に
V normalize(V v) {
    return v / sqrtl(dot(v, v));
}


std::tuple<V, V> tangentSpace(const V& n) {
    const double s = std::copysign(1, n.z);
    const double a = -1 / (s + n.z);
    const double b = n.x*n.y*a;
    return {
        V(1 + s * n.x*n.x*a,s*b,-s * n.x), 
        V(b,s + n.y*n.y*a,-n.y)
    };
}


// 16:20
// 浮動小数点数 [0, 1] を [0, 255] に変換
// ガンマ補正 ?
// OK
int tonemap(double v) {
    return std::min(std::max(int(std::pow(v, 1/2.2)*255), 0), 255);
};

// 26:01
// [0, 1) の一様乱数生成
struct Random {
    std::mt19937 engine;
    std::uniform_real_distribution<double> dist;
    Random() {};
    Random(int seed) {
        engine.seed(seed);
        dist.reset();
    }
    double next() { return dist(engine); }
};

// 05:10
// OK
struct Ray {
    V o;    // 原点 (origin)
    V d;    // 方向 (direction)
};

struct Sphere;
// 08:20
// Hit : 交差点の情報を表す
// OK
struct Hit {
    double t;   // レイの原点から交差点までの距離
    V p;
    V n;
    const Sphere* sphere;   // 当たった球へのポインタ
};

// シーン内のオブジェクト（球）
// OK
struct Sphere {
    V p;        // p : 中心座標 (point)
    double r;   // r : 球半径   (radius)
    // 22:35
    V R;        // R : 反射率   (reflectance)
    // 26:17
    V Le;       // Le : 照度    (illuminance)
    // 08:50
    std::optional<Hit> intersect(const Ray& ray, double tmin, double tmax) const {
        // 11:03
        const V op = p - ray.o;
        const double b = dot(op, ray.d);
        const double det = b * b - dot(op, op) + r * r;
        if (det < 0) { return {}; }
        
        const double t1 = b - sqrt(det);
        if (tmin < t1 && t1 < tmax) { return Hit{t1, {}, {}, this}; }
        const double t2 = b + sqrt(det);
        if (tmin < t2 && t2 < tmax) { return Hit{t2, {}, {}, this}; }

        return {};
    }
};

// シーン内のオブジェクト
// OK
struct Scene {
    std::vector<Sphere> spheres{
        // 17:19
        // 球2つでオーバーラップするときをテスト
        // 1
        /*{ V(-.5, 0, 0), 1, V(1, 0, 0)},  // Rの追加22:58
        { V( .5, 0, 0), 1, V(0, 1, 0)}   // Rの追加22:58*/

        // 2
        { V(1e5+1.0, 40.8, 81.6), 1e5, V(.75, .25, .25) },    // left
        { V(-1e5+99, 40.8, 81.6), 1e5, V(.25, .25, .75) },    // right
        { V(50, 40.8, 1e5),       1e5, V(.75, .75, .75) },    // front
        // { V(50, 40.8, -1e5 + 170), 1e5 },
        { V(50, 1e5, 81.6),       1e5, V(.75) },              // bottom
        { V(50, -1e5+81.6, 81.6), 1e5, V(.75) },              // top
        { V(27, 16.5, 47),       16.5, V(.999) },             // left - ball
        { V(73, 16.5, 78),       16.5, V(.999) },             // right - ball
        { V(50, 681.6-.27, 81.6), 600, V(), V(12) },          // light

        // 3
        // { V(0, 1, 0), 1  },
        // { V(-4, 1, 0), 1},
        // { V(4, 1, 0) , 1}
    };

    // 09:16
    // 交差したらHit型の値を返し、交差しなければstd::nullopt（交差していない情報）を返す
    // レイの存在範囲：原点からの距離[tmin, tmax]
    std::optional<Hit> intersect(const Ray& ray, double tmin, double tmax) const {
        // 各球に対してintersect関数を呼ぶ
        // tが最小のものを選択して返す
        std::optional<Hit> minh;
        for(const auto& sphere : spheres) {
            const auto h = sphere.intersect(ray, tmin, tmax);
            if(!h) { continue; }
            minh = h;
            tmax = minh->t;
        }
        if (minh) {
            const auto* s = minh->sphere;
            minh->p = ray.o + ray.d * minh->t;
            minh->n = (minh->p - s->p) / s->r;
        }
        return minh;        
    }
};

int main() {
    // Image size
    const int w = 1200;
    const int h = 800;

    // Samples per pixel
    const int spp = 100;
    // 18:47
    // Camera parameters
    // 1
    /*const V eye(5, 5, 5);       // カメラ位置
    const V center(0, 0, 0);    // 注視点
    const V up(0, 1, 0);        // 上を表すvector
    const double fov = 30*M_PI/180;*/      // 垂直方向視野角 (radian)

    // 2
    const V eye(50, 52, 295.6);
    const V center = eye + V(0, -0.042612, -1);
    const V up(0, 1, 0);
    const double fov = 30 * M_PI / 180;

    /*// 3
    const V eye(13, 2, 3);
    const V center;
    const V up(0, 1, 0);
    const double fov = 20 * M_PI / 180;
    */
   
    const double aspect = (double)w/(double)h;  // 画面のアスペクト比(注釈に基づきdoubleキャスト)
    // Basis vectors for camera coordinates
    // カメラ座標系の基底計算
    const auto wE = normalize(eye - center);
    const auto uE = normalize(cross(up, wE));
    const auto vE = cross(wE, uE);

    Scene scene;
    std::vector<V> I(w*h);
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < w * h; i++) {
        thread_local Random rng(42 + omp_get_thread_num());
        for (int j = 0; j < spp; j++) {
            // 05:54
            // z 軸平行なレイ
            const int x = i % w;
            const int y = h - i / w;
            Ray ray;
            // 19:59
            // ray.o = V(2.*(double)x/(double)w-1.0, 2.*(double)y/(double)h-1.0, 5.0);
            // ray.d = V(0, 0, -1);
            ray.o = eye;
            ray.d = [&]() {
                const double tf = std::tan(fov*.5);
                const double rpx = 2.*(x+rng.next())/w-1.0;
                const double rpy = 2.*(y+rng.next())/h-1.0;
                const V w = normalize(V(aspect * tf * rpx, tf * rpy, -1));
                return uE * w.x + vE * w.y + wE * w.z;
            }();

        V L(0), th(1);        
        for(int depth = 0; depth < 10; depth++) {
            // Intersection
            const auto h = scene.intersect(ray, 1e-4, 1e+10);
            if(!h) {
                break;
            }
            // Add contribution
            L = L + th * h->sphere->Le;
            // Update next direction
            ray.o = h->p;
            ray.d = [&]() {
                // Sample direction in local coordinates
                const auto n = dot(h->n, -ray.d) > 0 ? h->n : -h->n;
                const auto& [u, v] = tangentSpace(n);
                const auto d = [&]() {
                    const double r = sqrt(rng.next());
                    const double t = 2*M_PI*rng.next();
                    const double x = r*cos(t);
                    const double y = r*sin(t);
                    return V(x, y,
                            std::sqrt(
                                std::max(.0,1-x*x-y*y)));
                }();
                // Convert to world coordinates
                return u*d.x+v*d.y+n*d.z;
            }();
            // Update throughput
            th = th * h->sphere->R;
            if (std::max({ th.x,th.y,th.z }) == 0) {
                break;
            }
        }
        I[i] = I[i] + L / spp;
        }
    }
    std::ofstream ofs("result.ppm");
    // P3フォーマット
    // 1行目: P3
    // 2行目: WIDTH HEIGHT
    // 3行目: MAX_VALUE
    // 4行目以降: 配列[R, G, B], 値域[0, 255]
    ofs << "P3\n" << w << " " << h << "\n255\n";
    for (const auto& i : I) {
        ofs << tonemap(i.x) << " "
            << tonemap(i.y) << " "
            << tonemap(i.x) << "\n";
    }
    return 0;
}