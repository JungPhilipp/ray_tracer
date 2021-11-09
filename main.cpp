// smallpt, a Path Tracer by Kevin Beason, 2008
//        Remove "-fopenmp" for g++ version < 4.2
// Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
// Usage: time ./smallpt 5000 && xv image.ppm
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct Vec
{
    // position, also color (r,g,b)
    double x = 0;
    double y = 0;
    double z = 0;
    constexpr Vec operator+(const Vec& other) const
    {
        return Vec{ x + other.x, y + other.y, z + other.z };
    }
    constexpr Vec operator-(const Vec& other) const
    {
        return Vec{ x - other.x, y - other.y, z - other.z };
    }
    constexpr Vec operator*(double factor) const
    {
        return Vec{ x * factor, y * factor, z * factor };
    }
    constexpr Vec mult(const Vec& other) const
    {
        return Vec{ x * other.x, y * other.y, z * other.z };
    }
    constexpr Vec& normalize()
    {
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }
    constexpr double dot(const Vec& other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }
    // cross
    constexpr Vec operator%(Vec& other)
    {
        return Vec{ y * other.z - z * other.y,
                    z * other.x - x * other.z,
                    x * other.y - y * other.x };
    }
};
struct Ray
{
    Vec orientation;
    Vec direction;
};
enum class Refl_t
{
    DIFF,
    SPEC,
    REFR
};

// material types, used in radiance()
struct Sphere
{
    double radius;
    Vec position;
    Vec emission;
    Vec color;
    Refl_t reflection;

    constexpr double intersect(const Ray& ray) const
    { // returns distance, 0 if nohit
        Vec op =
          position - ray.orientation; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        const double epsilon = 1e-4;
        const double b = op.dot(ray.direction);
        double det = b * b - op.dot(op) + radius * radius;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);

        const double t = b - det;
        if (t > epsilon)
            return t;
        else
        {
            const double t1 = b + det;
            if (t1 > epsilon)
                return t1;
            else
                return 0;
        }
    }
};
constexpr Sphere spheres[] = {
    // Scene: radius, position, emission, color, material
    Sphere{ 1e5,
            Vec{ 1e5 + 1, 40.8, 81.6 },
            Vec{},
            Vec{ .75, .25, .25 },
            Refl_t::DIFF }, // Left
    Sphere{ 1e5,
            Vec{ -1e5 + 99, 40.8, 81.6 },
            Vec{},
            Vec{ .25, .25, .75 },
            Refl_t::DIFF }, // Rght
    Sphere{ 1e5,
            Vec{ 50, 40.8, 1e5 },
            Vec(),
            Vec{ .75, .75, .75 },
            Refl_t::DIFF },                                                 // Back
    Sphere{ 1e5, Vec{ 50, 40.8, -1e5 + 170 }, Vec{}, Vec(), Refl_t::DIFF }, // Frnt
    Sphere{ 1e5,
            Vec{ 50, 1e5, 81.6 },
            Vec(),
            Vec{ .75, .75, .75 },
            Refl_t::DIFF }, // Botm
    Sphere{ 1e5,
            Vec{ 50, -1e5 + 81.6, 81.6 },
            Vec{},
            Vec{ .75, .75, .75 },
            Refl_t::DIFF }, // Top
    Sphere{ 16.5,
            Vec{ 27, 16.5, 47 },
            Vec{},
            Vec{ 1, 1, 1 } * .999,
            Refl_t::SPEC }, // Mirr
    Sphere{ 16.5,
            Vec{ 73, 16.5, 78 },
            Vec{},
            Vec{ 1, 1, 1 } * .999,
            Refl_t::REFR }, // Glas
    Sphere{ 600,
            Vec{ 50, 681.6 - .27, 81.6 },
            Vec{ 12, 12, 12 },
            Vec(),
            Refl_t::DIFF } // Lite
};
constexpr double
clamp(double x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}
constexpr int
toInt(double x)
{
    return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}
constexpr bool
intersect(const Ray& r, double& t, int& id)
{
    double n = sizeof(spheres) / sizeof(Sphere);
    double d = 0;
    double inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}
constexpr Vec
radiance(const Ray& r, int depth, unsigned short* Xi)
{
    double t = 0; // distance to intersection
    int id = 0;   // id of intersected object
    if (!intersect(r, t, id))
        return Vec();                // if miss, return black
    const Sphere& obj = spheres[id]; // the hit object
    Vec x = r.orientation + r.direction * t, n = (x - obj.position).normalize(),
        nl = n.dot(r.direction) < 0 ? n : n * -1, f = obj.color;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
    if (++depth > 5)
    {
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return obj.emission; // R.R.
    }
    if (obj.reflection == Refl_t::DIFF)
    { // IdealRefl_t::DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec{ 0, 1 } : Vec{ 1 }) % w).normalize(),
            v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalize();
        return obj.emission + f.mult(radiance(Ray{ x, d }, depth, Xi));
    }
    else if (obj.reflection == Refl_t::SPEC) // IdealRefl_t::SPECULAR reflection
        return obj.emission +
               f.mult(
                 radiance(Ray{ x, r.direction - n * 2 * n.dot(r.direction) }, depth, Xi));
    Ray reflRay{
        x, r.direction - n * 2 * n.dot(r.direction)
    };                         // Ideal dielectricRefl_t::REFRACTION
    bool into = n.dot(nl) > 0; // Ray from outside going in?
    double nc = 1;
    double nt = 1.5;
    double nnt = into ? nc / nt : nt / nc;
    double ddn = r.direction.dot(nl);
    double cos2t = 0;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
        return obj.emission + f.mult(radiance(reflRay, depth, Xi));
    Vec tdir =
      (r.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b),
           c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re,
           RP = Re / P, TP = Tr / (1 - P);
    return obj.emission +
           f.mult(depth > 2 ? (erand48(Xi) < P ? // Russian roulette
                                 radiance(reflRay, depth, Xi) * RP
                                               : radiance(Ray{ x, tdir }, depth, Xi) * TP)
                            : radiance(reflRay, depth, Xi) * Re +
                                radiance(Ray{ x, tdir }, depth, Xi) * Tr);
}
int
main(int argc, char* argv[])
{
    int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : 1;     // # samples
    Ray cam{ Vec{ 50, 52, 295.6 }, Vec{ 0, -0.042612, -1 }.normalize() }; // cam pos, dir
    Vec cx = Vec{ w * .5135 / h }, cy = (cx % cam.direction).normalize() * .5135, r,
        *c = new Vec[w * h];
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
    for (int y = 0; y < h; y++)
    { // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = { 0, 0, y * y * y }; x < w; x++) // Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * erand48(Xi),
                               dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi),
                               dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.direction;
                        d.normalize();
                        r = r + radiance(Ray{ cam.orientation + d * 140, d }, 0, Xi) *
                                  (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec{ clamp(r.x), clamp(r.y), clamp(r.z) } * .25;
                }
    }
    FILE* f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
