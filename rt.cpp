#include <iostream>
#include <fstream>
#include <vector>
// 9:00
// 値のほかに無効の情報を持てる
#include <optional>

// sqrtが使えないので...
#include <cmath>

// 03:48
// vector構造体
// 位置, 方向, 色を表す
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
    return v / sqrt(dot(v, v));
}

// 16:20
// 浮動小数点数 [0, 1] を [0, 255] に変換
// ガンマ補正 ?
int tonemap(double v) {
    return std::min(std::max(int(std::pow(v, 1/2.2)*255), 0), 255);
};


// 05:10
struct Ray {
    V o;    // 原点 (origin)
    V d;    // 方向 (direction)
};

struct Sphere;
// 08:20
// Hit : 交差点の情報を表す
struct Hit {
    double t;   // レイの原点から交差点までの距離
    V p;
    V n;
    const Sphere* sphere;   // 当たった球へのポインタ
};

// シーン内のオブジェクト（球）
struct Sphere {
    V p;        // p : 中心座標 (point)
    double r;   // r : 球半径   (radius)
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
struct Scene {
    std::vector<Sphere> spheres{
        // 17:19
        // 球2つでオーバーラップするときをテスト
        // 1
        // { V(-.5, 0, 0), 1 },
        // { V( .5, 0, 0), 1 }

        // 2
        { V(1e5 + 1, 40.8, 81.6),   1e5 },
        { V(-1e5 + 99, 40.8, 81.6), 1e5 },
        { V(50, 40.8, 1e5),         1e5 },
        // { V(50, 40.8, -1e5 + 170), 1e5 },
        { V(50, 1e5, 81.6),         1e5 },
        { V(50, -1e5 + 81.6, 81.6), 1e5 },
        { V(27, 16.5, 47),         16.5 },
        { V(73, 16.5, 78),         16.5 },
        { V(50, 93, 50),           20.5 }
        // { V(50, 681.6-.27, 81.6), 600 },
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
    // 18:47
    // Camera parameters
    /*// 1
    const V eye(5, 5, 5);       // カメラ位置
    const V center(0, 0, 0);    // 注視点
    const V up(0, 1, 0);        // 上を表すvector
    const double fov = 30*M_PI/180;      // 垂直方向視野角 (radian)
    */

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
    std::ofstream ofs("result.ppm");
    // P3フォーマット
    // 1行目: P3
    // 2行目: WIDTH HEIGHT
    // 3行目: MAX_VALUE
    // 4行目以降: 配列[R, G, B], 値域[0, 255]
    ofs << "P3\n" << w << " " << h << "\n255\n";
    for(int i = 0; i < w * h; i++) {
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
            const double rpx = 2.*(double)x/(double)w-1;
            const double rpy = 2.*(double)y/(double)h-1;
            const V w = normalize(V(aspect * tf * rpx, tf * rpy, -1));
            return uE * w.x + vE * w.y + wE * w.z;
        }();

        // 12:00
        const auto h = scene.intersect(ray, 0, 1e+10);
        if (h) {
            const auto n = h->n;
            ofs << tonemap(n.x) << " "
                << tonemap(n.y) << " "
                << tonemap(n.z) << "\n";
            // マゼンタ
            // ofs << "0 0 255\n";
        } else {
            // 黒
            ofs << "0 0 0\n";
        }
    }
    return 0;
}