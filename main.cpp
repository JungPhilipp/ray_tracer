// smallpt, a Path Tracer by Kevin Beason, 2008
//        Remove "-fopenmp" for g++ version < 4.2
// Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
// Usage: time ./smallpt 5000 && xv image.ppm
#include <array>
#include <cmath>
#include <concepts>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <vector>

using precession_t = double;

struct Vec
{
    // position, also color (r,g,b)
    precession_t x = 0;
    precession_t y = 0;
    precession_t z = 0;
    constexpr auto operator+(const Vec& other) const -> Vec
    {
        return { x + other.x, y + other.y, z + other.z };
    }
    constexpr auto operator-(const Vec& other) const -> Vec
    {
        return { x - other.x, y - other.y, z - other.z };
    }
    constexpr auto operator*(precession_t factor) const -> Vec
    {
        return { x * factor, y * factor, z * factor };
    }
    constexpr auto mult(const Vec& other) const -> Vec
    {
        return { x * other.x, y * other.y, z * other.z };
    }
    constexpr auto& normalize()
    {
        // TODO implement fast inverse sqrt
        return *this = *this * (1 / std::sqrt(x * x + y * y + z * z));
    }
    constexpr auto dot(const Vec& other) const -> precession_t
    {
        return x * other.x + y * other.y + z * other.z;
    }
    // cross
    constexpr auto operator%(const Vec& other) const -> Vec
    {
        return { y * other.z - z * other.y,
                 z * other.x - x * other.z,
                 x * other.y - y * other.x };
    }
};
struct Ray
{
    Vec orientation;
    Vec direction;
};
enum class Reflection
{
    Diffuse,
    Specular,
    Refraction
};

// material types, used in radiance()
struct Sphere
{
    precession_t radius;
    Vec position;
    Vec emission;
    Vec color;
    Reflection reflection;
};

/**
 * @brief Returns distance, 0 if not hit
 */
constexpr auto
intersect(const Sphere& sphere, const Ray& ray) -> precession_t
{
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    const auto op = Vec{ sphere.position - ray.orientation };
    const auto epsilon = 1e-4;
    const auto b = op.dot(ray.direction);
    auto det = b * b - op.dot(op) + sphere.radius * sphere.radius;
    if (det < 0)
        return 0;
    else
        det = std::sqrt(det);

    const auto t = b - det;
    if (t > epsilon)
        return t;
    else
    {
        const auto t1 = b + det;
        if (t1 > epsilon)
            return t1;
        else
            return 0;
    }
}
constexpr auto spheres = std::array{
    // Left
    Sphere{ 1e5, { 1e5 + 1, 40.8, 81.6 }, {}, { .75, .25, .25 }, Reflection::Diffuse },
    // Right
    Sphere{ 1e5, { -1e5 + 99, 40.8, 81.6 }, {}, { .25, .25, .75 }, Reflection::Diffuse },
    // Back
    Sphere{ 1e5, { 50, 40.8, 1e5 }, {}, { .75, .75, .75 }, Reflection::Diffuse },
    // Front
    Sphere{ 1e5, { 50, 40.8, -1e5 + 170 }, {}, {}, Reflection::Diffuse },
    // Bottom
    Sphere{ 1e5, { 50, 1e5, 81.6 }, {}, { .75, .75, .75 }, Reflection::Diffuse },
    // Top
    Sphere{ 1e5, { 50, -1e5 + 81.6, 81.6 }, {}, { .75, .75, .75 }, Reflection::Diffuse },
    // Mirror
    Sphere{ 16.5, { 27, 16.5, 47 }, {}, Vec{ 1, 1, 1 } * .999, Reflection::Specular },
    // Glass
    Sphere{ 16.5, { 73, 16.5, 78 }, {}, Vec{ 1, 1, 1 } * .999, Reflection::Refraction },
    // Lite
    Sphere{ 600, { 50, 681.6 - .27, 81.6 }, { 12, 12, 12 }, {}, Reflection::Diffuse }
};
template<std::floating_point T>
constexpr auto
clamp(T x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}
template<std::floating_point T>
constexpr auto
toInt(T x)
{
    return static_cast<int>(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

constexpr auto
intersect_with_spheres(const Ray& ray) -> std::tuple<bool, precession_t, size_t>
{
    const auto inf = 1e20;
    auto t = inf;
    auto id = size_t(0);
    auto i = size_t(0);
    for (auto sphere = rbegin(spheres); sphere != rend(spheres); sphere++, i++)
    {
        if (auto distance = intersect(spheres[i], ray); distance != 0 && distance < t)
        {
            t = distance;
            id = i;
        }
    }
    return { t < inf, t, id };
}
constexpr auto
radiance(const Ray& r, int depth, unsigned short* Xi) -> Vec
{
    const auto [intersects, distance, intersected_object_idx] = intersect_with_spheres(r);
    if (!intersects)
        return {};                                     // if miss, return black
    const auto& obj = spheres[intersected_object_idx]; // the hit object

    const auto x = r.orientation + r.direction * distance;
    const auto n = (x - obj.position).normalize();
    const auto nl = n.dot(r.direction) < 0 ? n : n * -1;
    auto f = obj.color;
    const auto p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
    if (++depth > 5)
    {
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return obj.emission; // R.R.
    }
    switch (obj.reflection)
    {
        case Reflection::Diffuse:
        {
            const auto r1 = 2 * M_PI * erand48(Xi);
            const auto r2 = erand48(Xi);
            const auto r2s = std::sqrt(r2);
            const auto w = nl;
            const auto u =
              ((std::abs(w.x) > .1 ? Vec{ 0, 1 } : Vec{ 1 }) % w).normalize();
            const auto v = w % u;
            const auto d =
              (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalize();
            return obj.emission + f.mult(radiance(Ray{ x, d }, depth, Xi));
        }
        case Reflection::Specular:
        {
            return obj.emission +
                   f.mult(radiance(
                     Ray{ x, r.direction - n * 2 * n.dot(r.direction) }, depth, Xi));
        }
        case Reflection::Refraction:
        {
            Ray reflRay{
                x, r.direction - n * 2 * n.dot(r.direction)
            };                               // Ideal dielectric REFRACTION
            const auto into = n.dot(nl) > 0; // Ray from outside going in?
            const auto nc = 1;
            const auto nt = 1.5;
            const auto nnt = into ? nc / nt : nt / nc;
            const auto ddn = r.direction.dot(nl);
            const auto cos2t = 0;
            // Total internal reflection
            if (const auto cos2t = 1 - nnt * nnt * (1 - ddn * ddn); cos2t < 0)
                return obj.emission + f.mult(radiance(reflRay, depth, Xi));

            const auto tdir =
              (r.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + std::sqrt(cos2t))))
                .normalize();

            const auto a = nt - nc;
            const auto b = nt + nc;
            const auto R0 = a * a / (b * b);
            const auto c = 1 - (into ? -ddn : tdir.dot(n));
            const auto Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re,
                       P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
            return obj.emission +
                   f.mult(depth > 2
                            ? (erand48(Xi) < P ? // Russian roulette
                                 radiance(reflRay, depth, Xi) * RP
                                               : radiance(Ray{ x, tdir }, depth, Xi) * TP)
                            : radiance(reflRay, depth, Xi) * Re +
                                radiance(Ray{ x, tdir }, depth, Xi) * Tr);
        }
    }
}
auto
main(int argc, char* argv[]) -> int
{
    constexpr auto width = 1024;
    constexpr auto height = 768;
    // # samples
    const auto samples = argc == 2 ? std::stoi(argv[1]) / 4 : 1;
    const auto cam = Ray{ Vec{ 50, 52, 295.6 }, Vec{ 0, -0.042612, -1 }.normalize() };

    constexpr auto cx = Vec{ width * .5135 / height };
    const auto cy = (cx % cam.direction).normalize() * .5135;
    auto c = std::vector<Vec>(width * height);

#pragma omp parallel for schedule(dynamic, 1) // OpenMP
    for (int y = 0; y < height; y++)
    { // Loop over image rows
        fprintf(
          stderr, "\rRendering (%d spp) %5.2f%%", samples * 4, 100. * y / (height - 1));
        for (unsigned short x = 0, Xi[3] = { 0, 0, y * y * y }; x < width;
             x++) // Loop cols
            for (int sy = 0, i = (height - y - 1) * width + x; sy < 2;
                 sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++)
                { // 2x2 subpixel cols
                    auto r = Vec{};
                    for (int s = 0; s < samples; s++)
                    {
                        const auto r1 = 2 * erand48(Xi),
                                   dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        const auto r2 = 2 * erand48(Xi),
                                   dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        auto d = Vec{ cx * (((sx + .5 + dx) / 2 + x) / width - .5) +
                                      cy * (((sy + .5 + dy) / 2 + y) / height - .5) +
                                      cam.direction };
                        d.normalize();
                        r = r + radiance(Ray{ cam.orientation + d * 140, d }, 0, Xi) *
                                  (1. / samples);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec{ clamp(r.x), clamp(r.y), clamp(r.z) } * .25;
                }
    }
    FILE* f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
    for (const auto& color : c)
        fprintf(f, "%d %d %d ", toInt(color.x), toInt(color.y), toInt(color.z));
}
