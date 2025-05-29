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

struct SPHParams {
    float x_limit, y_limit, z_limit;
    float rho0 = 1000.0; // [kg/m^3]
    float h = 0.026; // [m]    kernel redius.
    float mass = 8.0e-3;    // [kg] particle mass
    float stiffness = 10000;    // [Pa] 
    float dt = 0.001;   // [s]  simulation time step.
};

struct SPHSimulator {
    SPHSimulator(size_t num_particles):
    num_particles(num_particles)
    {}
        //x_limit(VIEW_WIDTH), y_limit(VIEW_HEIGHT), z_limit(VIEW_WIDTH),
        //h(H), rho0(REST_DENS), mass(MASS), stiffness(GAS_CONST), dt(DT)
    //float x_limit, y_limit, z_limit;
    //float rho0 = 1000.0; // [kg/m^3]
    //float h = 0.026; // [m]    kernel redius.
    //float mass = 8.0e-3;    // [kg] particle mass
    //float stiffness = 10000;    // [Pa] 
    //float dt = 0.001;   // [s]  simulation time step.


    void init_simulation(bool serialize_coordinate = false){
        this->particles.clear();
        //const size_t num_particles_limit = 3000;
        //float z = 0.0;
        srand(0);
        for(float y = EPS; y < VIEW_HEIGHT - EPS * 2.0; y += H) {
            for(float z = VIEW_WIDTH/4; z <= VIEW_WIDTH/2; z+=H) {

                for(float x = VIEW_WIDTH/ 4; x <= VIEW_WIDTH / 2; x += H) {
                //for(float x = VIEW_WIDTH/ 4; x <= VIEW_WIDTH *3 / 4; x += H) {
                    if (this->particles.size() < this->num_particles) {
                        //float jitter = static_cast<float>(arc4random())/static_cast<float>(RAND_MAX);
                        float jitter = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
                        float jitter2 = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
                        this->particles.push_back(Particle(x+jitter, y, z+jitter2));
                    } else {
                        break;
                    }
                }

            }
        }
        std::printf("From SPHSimulator::init_simulation: %d particles initialized", particles.size());
        std::cout << particles.size() << std::endl;
        if (serialize_coordinate == true) {
            this->update_serialized_coordinate();
        }
        return;
    }
    void step(const size_t n_step = 1, const bool serialize_coordinate = false) {
        for(size_t i = 0; i < n_step; i++) {
            this->compute_density_pressure();
            this->compute_forces();
            this->integrate();
        }
        if (serialize_coordinate == true) {
            this->update_serialized_coordinate();
        }
    }
    void compute_density_pressure(void) {
        #pragma omp parallel for
        for(size_t i = 0; i < this->particles.size(); i++) {
            Particle &pi = this->particles[i];
            pi.rho = 0.0;
            for(auto &pj : this->particles) {
                Eigen::Vector3d rij = pj.r - pi.r;
                float r2 = rij.squaredNorm();

                if (r2 < HSQ) {
                    pi.rho += MASS * POLY6 * pow(HSQ-r2, 3.f);
                }
            }
            pi.p = GAS_CONST * (pi.rho - REST_DENS);
        }
    }
    void compute_forces(void) {
        #pragma omp parallel for
        for(size_t i = 0; i < this->particles.size(); i++) {
            Particle &pi = this->particles[i];
            Eigen::Vector3d fpress(0.f, 0.f, 0.f);
            Eigen::Vector3d fvisc(0.f, 0.f, 0.f);
            for(auto &pj : this->particles) {
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
    void integrate(void) {
        #pragma omp parallel for
        for(size_t i = 0; i < this->particles.size(); i++) {
            Particle &p = this->particles[i];
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
    size_t get_num_particles() {
        return this->particles.size();
    }
    void update_serialized_coordinate(void) {
        this->serialized_coordinates.clear();
        for(const auto &ptcl: this->particles) {
            this->serialized_coordinates.push_back(float(ptcl.r[0]));
            this->serialized_coordinates.push_back(float(ptcl.r[1]));
            this->serialized_coordinates.push_back(float(ptcl.r[2]));
        }
    }
    uintptr_t get_serialized_coordinate_buffer_ptr() {
        return reinterpret_cast<uintptr_t>(this->serialized_coordinates.data());
    }
    size_t get_serialized_coordinate_buffer_size() {
        return this->serialized_coordinates.size();
    }
    std::vector<Particle> particles;
    size_t num_particles;
    std::vector<float> serialized_coordinates;
};

void InitSPH(void)
{
    srand(0);
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


bool compare_results(const SPHSimulator &sim, const std::vector<Particle> &particles) {
    bool ret = true;
    if (sim.particles.size() != particles.size()) {
        std::cout << "size not matched" << std::endl;
        return false;
    }
    for(size_t i = 0; i < sim.particles.size(); i++) {
        if (sim.particles[i].r != particles[i].r) {
            std::printf("index: %d\n", i);
            std::cout << sim.particles[i].r << std::endl;
            std::printf("----------\n");
            std::cout << particles[i].r << std::endl;
            ret = false;
            break;
        }
    }
    std::printf("All particles matched!\n");
    return ret;
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

    emscripten::class_<SPHSimulator>("SPHSimulator")
        .constructor<unsigned long>()
        .function("init_simulation", &SPHSimulator::init_simulation)
        .function("step", &SPHSimulator::step)
        .function("get_num_particles", &SPHSimulator::get_num_particles)
        .function("get_serialized_coordinate_buffer_ptr", &SPHSimulator::get_serialized_coordinate_buffer_ptr)
        .function("get_serialized_coordinate_buffer_size", &SPHSimulator::get_serialized_coordinate_buffer_size);
}
#else

int main(void) {
    set_num_particles(100);
    InitSPH();
    init_serialized_coordinate();
    std::cout << "ptr: " << get_serialized_coordinate_buffer_ptr() << std::endl;
    std::cout << "size: " << get_serialized_coordinate_buffer_size() << std::endl;

    SPHSimulator sim(num_particles_limit);
    sim.init_simulation();
    assert(compare_results(sim, particles));

    for(int j = 0; j < 10; j++) {
        for(int i = 0; i < 10; i++) {
            sim.compute_density_pressure();
            sim.compute_forces();
            sim.integrate();

            ComputeDensityPressure();
            ComputeForces();
            Integrate();
        }
        assert(compare_results(sim, particles));

    }
    return 0;
}
#endif
