#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>

#ifdef __EMSCRIPTEN__
# include <emscripten.h>
# include <emscripten/bind.h>
#endif

struct Particle
{
    Particle(float x, float y, float z): r(x,y,z), v(0.0, 0.0, 0.0), f(0.0, 0.0, 0.0), rho(0.0), p(0.0){;}
    Eigen::Vector3d r;
    Eigen::Vector3d v;
    Eigen::Vector3d f;
    float rho, p;
    double getX(){return r[0];}
    double getY(){return r[1];}
    //double getZ(){return float(r[2]);}
};

std::vector<Particle> particles;

std::vector<float> serialized_coordinate;
void init_serialized_coordinate(void){
    serialized_coordinate.clear();
    for(const auto &ptcl: particles) {
        serialized_coordinate.push_back(float(ptcl.r[0]));
        serialized_coordinate.push_back(float(ptcl.r[1]));
        serialized_coordinate.push_back(float(ptcl.r[2]));
    }
    std::printf("%s: size: %d", __func__, serialized_coordinate.size());
}

uintptr_t get_serialized_coordinate_buffer_ptr() {
    return reinterpret_cast<uintptr_t>(serialized_coordinate.data());
}
size_t get_serialized_coordinate_buffer_size() {
    return serialized_coordinate.size();
}
size_t get_num_particles() {
    return particles.size();
}
std::vector<float> rho_buffer;
void init_rho_vector(void) {
    rho_buffer.clear();
    for(const auto &ptcl: particles) {
        rho_buffer.push_back(float(ptcl.rho));
    }
}
size_t get_rho_buffer_ptr() {
    return reinterpret_cast<uintptr_t>(rho_buffer.data());
}
size_t get_rho_buffer_size() {
    return rho_buffer.size();
}


const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
const static float VIEW_WIDTH = 1.5 * float(WINDOW_WIDTH);
const static float VIEW_HEIGHT= 1.5 * float(WINDOW_HEIGHT);

const static Eigen::Vector3d G(0.f, -10.f, 0.f);
const static float REST_DENS = 300.f;
const static float GAS_CONST = 2000.f;  // const for eq of state
const static float H = 16.f;            // kernel radius
const static float HSQ = H*H;       // radius^2 for optimization
const static float MASS = 2.5f;     //assume all particles have the same mass
const static float VISC = 200.f;    // viscosity constant
const static float DT = 0.0007f;    // integration timestep

// smoothing kernels defined in Müller and their gradients
// adapted to 2D per "SPH Based Shallow Water Simulation" by Solenthaler et al.
// const static float POLY6 = 4.f / (M_PI * pow(H, 8.f));
// const static float SPIKY_GRAD = -10.f / (M_PI * pow(H, 5.f));
// const static float VISC_LAP = 40.f / (M_PI * pow(H, 5.f));

const static float POLY6 = 315.f / (64.f * M_PI * pow(H, 9.f));
const static float SPIKY_GRAD = -45.f / (M_PI * pow(H, 6.f));
const static float VISC_LAP = 45.f / (M_PI * pow(H, 6.f));

// simulation parameters
const static float EPS = H; // boundary epsilon
const static float BOUND_DAMPING = -0.5f;

size_t num_particles_limit = 100;
void set_num_particles(const size_t n){
    num_particles_limit = n;
}

void InitSPH(void)
{
    particles.clear();
    //const size_t num_particles_limit = 3000;
    //float z = 0.0;
    for(float y = EPS; y < VIEW_HEIGHT - EPS * 2.0; y += H) {
        for(float z = VIEW_WIDTH/4; z <= VIEW_WIDTH/2; z+=H) {

            for(float x = VIEW_WIDTH/ 4; x <= VIEW_WIDTH / 2; x += H) {
            //for(float x = VIEW_WIDTH/ 4; x <= VIEW_WIDTH *3 / 4; x += H) {
                if (particles.size() < num_particles_limit) {
                    //float jitter = static_cast<float>(arc4random())/static_cast<float>(RAND_MAX);
                    float jitter = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
                    float jitter2 = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
                    particles.push_back(Particle(x+jitter, y, z+jitter2));
                } else {
                    break;
                }
            }

        }
    }
    std::printf("From InitSPH: %d particles initialized", particles.size());
    std::cout << particles.size() << std::endl;
    return;
}
void ComputeDensityPressure(void)
{
    //for(auto &pi : particles) {
    #pragma omp parallel for
    for(size_t i = 0; i < particles.size(); i++) {
        Particle &pi = particles[i];
        pi.rho = 0.0;
        for(auto &pj : particles) {
            Eigen::Vector3d rij = pj.r - pi.r;
            float r2 = rij.squaredNorm();

            if (r2 < HSQ) {
                pi.rho += MASS * POLY6 * pow(HSQ-r2, 3.f);
            }
        }
        pi.p = GAS_CONST * (pi.rho - REST_DENS);
    }
}

void ComputeForces(void)
{
    //for(auto &pi : particles) {
    #pragma omp parallel for
    for(size_t i = 0; i < particles.size(); i++) {
        Particle &pi = particles[i];
        Eigen::Vector3d fpress(0.f, 0.f, 0.f);
        Eigen::Vector3d fvisc(0.f, 0.f, 0.f);
        for(auto &pj : particles) {
            if (&pi == &pj){
                continue;
            }
            Eigen::Vector3d rij = pj.r - pi.r;
            float r = rij.norm();
            if (r < H) {
                fpress += -rij.normalized() * MASS * (pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD * pow(H - r, 3.f);
                fvisc += VISC * MASS * (pj.v - pi.v) / pj.rho * VISC_LAP * (H-r);
            }
        }
        Eigen::Vector3d fgrav = G*MASS / pi.rho;
        pi.f = fpress + fvisc + fgrav;
    }
}

void Integrate(void)
{
    //for(auto &p : particles) {
    #pragma omp parallel for
    for(size_t i = 0; i < particles.size(); i++) {
        Particle &p = particles[i];
        p.v += DT * p.f / p.rho;
        p.r += DT * p.v;
        if (p.r(0) - EPS < 0.f){
            p.v(0) *= BOUND_DAMPING;
            p.r(0) = EPS;
        }
        if (p.r(0) + EPS > VIEW_WIDTH) {
            p.v(0) *= BOUND_DAMPING;
            p.r(0) = VIEW_WIDTH - EPS;
        }
        if (p.r(1) - EPS < 0.f){
            p.v(1) *= BOUND_DAMPING;
            p.r(1) = EPS;
        }
        if (p.r(1) + EPS > VIEW_HEIGHT) {
            p.v(1) *= BOUND_DAMPING;
            p.r(1) = VIEW_HEIGHT - EPS;
        }

        if (p.r(2) - EPS < 0.f){
            p.v(2) *= BOUND_DAMPING;
            p.r(2) = EPS;
        }
        if (p.r(2) + EPS > VIEW_WIDTH / 4) {
            p.v(2) *= BOUND_DAMPING;
            p.r(2) = VIEW_WIDTH/4 - EPS;
        }
    }
}

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(my_module) {
    emscripten::function("InitSPH", &InitSPH);
    emscripten::function("Integrate", &Integrate);
    emscripten::function("ComputeForces", &ComputeForces);
    emscripten::function("ComputeDensityPressure", &ComputeDensityPressure);
    emscripten::function("init_serialized_coordinate", &init_serialized_coordinate);
    emscripten::function("get_serialized_coordinate_buffer_ptr", &get_serialized_coordinate_buffer_ptr);
    emscripten::function("get_serialized_coordinate_buffer_size", &get_serialized_coordinate_buffer_size);
    emscripten::function("get_num_particles", &get_num_particles);
    emscripten::function("set_num_particles", &set_num_particles);
    emscripten::class_<Particle>("Particle")
        .constructor<float, float, float>()
        .function("getX", &Particle::getX)
        .function("getY", &Particle::getY);
}
#else
int main(void) {
    set_num_particles(100);
    InitSPH();
    init_serialized_coordinate();
    std::cout << "ptr: " << get_serialized_coordinate_buffer_ptr() << std::endl;
    std::cout << "size: " << get_serialized_coordinate_buffer_size() << std::endl;
}
#endif
