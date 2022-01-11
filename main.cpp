#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

float frand() {
    return (float)std::rand() / RAND_MAX * 2 - 1;
}

const size_t N = 48;

struct Stars {
    std::array<float, N> px, py, pz;
    std::array<float, N> vx, vy, vz;
    std::array<float, N> mass;
}stars;

void init() {
    for (size_t i = 0; i < 48; i++) {
        stars.px[i]=frand();
        stars.py[i]=frand();
        stars.pz[i]=frand();
        stars.vx[i]=frand();
        stars.vy[i]=frand();
        stars.vz[i]=frand();
        stars.mass[i]=frand()+1;
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;

void step() {
    for (size_t i = 0; i < N; i++) {
        float px_tmp = stars.px[i];
        float py_tmp = stars.py[i];
        float pz_tmp = stars.pz[i];
        float vx_tmp = 0.0f;
        float vy_tmp = 0.0f;
        float vz_tmp = 0.0f;
        for (size_t j = 0; j < N; j++) {
            float dx = stars.px[j] - px_tmp;
            float dy = stars.py[j] - py_tmp;
            float dz = stars.pz[j] - pz_tmp;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            d2 *= std::sqrt(d2);
            vx_tmp += dx * stars.mass[j] * G * dt / d2;
            vy_tmp += dy * stars.mass[j] * G * dt / d2;
            vz_tmp += dz * stars.mass[j] * G * dt / d2;
        }
        stars.vx[i] += vx_tmp;
        stars.vy[i] += vy_tmp;
        stars.vz[i] += vz_tmp;
    }
    for (size_t i = 0; i < N; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for (size_t i = 0; i < N; i++) {
        float px_tmp = stars.px[i];
        float py_tmp = stars.py[i];
        float pz_tmp = stars.pz[i];
        float mass_tmp = stars.mass[i];
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2.0f;
        for (size_t j = 0; j < N; j++) {
            float dx = stars.px[j] - px_tmp;
            float dy = stars.py[j] - py_tmp;
            float dz = stars.pz[j] - pz_tmp;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars.mass[j] * mass_tmp * G / std::sqrt(d2) / 2.0f;
        }
    }
    return energy;
}

template <class Func>
long benchmark(Func const &func) {
    auto t0 = std::chrono::steady_clock::now();
    func();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    return dt.count();
}

int main() {
    init();
    printf("Initial energy: %f\n", calc());
    auto dt = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
