#include <iostream>
#include <fstream>

int main() {
    const int w = 1200;
    const int h = 800;
    std::ofstream ofs("result.ppm");
    // P3フォーマット
    // 1行目: P3
    // 2行目: WIDTH HEIGHT
    // 3行目: MAX_VALUE
    // 4行目以降: 配列[R, G, B], 値域[0, 255]
    ofs << "P3\n" << w << " " << h << "\n255\n";
    for(int i = 0; i < w * h; i++){
        // マゼンタ
        ofs << "255 0 255\n";
    }
    return EXIT_SUCCESS;
}