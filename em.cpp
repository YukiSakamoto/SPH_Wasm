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


const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
const static float VIEW_WIDTH = 1.5 * float(WINDOW_WIDTH);
const static float VIEW_HEIGHT= 1.5 * float(WINDOW_HEIGHT);

const static Eigen::Vector3d G(0.f, -10.f, 0.f);
const static float REST_DENS = 300.f;
//const static float GAS_CONST = 2000.f;  // const for eq of state
const static float GAS_CONST = 2000.f;  // const for eq of state
const static float H = 16.f;            // kernel radius
const static float HSQ = H*H;       // radius^2 for optimization
const static float MASS = 2.5f;     //assume all particles have the same mass
//const static float VISC = 200.f;    // viscosity constant
const static float VISC = 300.f;    // viscosity constant
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
const static float BOUND_DAMPING = -0.8f;


struct SPHSimulator {
    SPHSimulator(size_t num_particles):
    num_particles(num_particles), dt(DT), time(0.0), 
    x_limit(VIEW_WIDTH), y_limit(VIEW_HEIGHT), z_limit(VIEW_WIDTH),
    rho0(REST_DENS), mass(MASS), stiffness(GAS_CONST), h(H), dx(H*0.5), epsilon(EPS), visc(VISC)
    {}

    void init_simulation(bool serialize_coordinate = false){
        this->time = 0;
        this->particles.clear();
        //const size_t num_particles_limit = 3000;
        //float z = 0.0;
        srand(0);
        //for(float y = this->epsilon; y < this->y_limit - this->epsilon * 2.0; y += dx) {
        for(float y = this->epsilon; y < this->y_limit - this->epsilon * 2.0; y += dx) {
            for(float z = this->epsilon; z <= this->z_limit/2; z+=dx) {

                for(float x = this->epsilon; x <= this->x_limit / 2; x += dx) {
                //for(float x = VIEW_WIDTH/ 4; x <= VIEW_WIDTH *3 / 4; x += H) {
                    if (this->particles.size() < this->num_particles) {
                        //float jitter = static_cast<float>(arc4random())/static_cast<float>(RAND_MAX);
                        float jitter = this->dx * static_cast<float>(rand())/static_cast<float>(RAND_MAX);
                        float jitter2 = this->dx * static_cast<float>(rand())/static_cast<float>(RAND_MAX);
                        this->particles.push_back(Particle(x+jitter, y, z+jitter2));
                    } else {
                        break;
                    }
                }

            }
        }
        this->define_kernel_constants();
        std::printf("From SPHSimulator::init_simulation: %d particles initialized", particles.size());
        std::cout << particles.size() << std::endl;
        if (serialize_coordinate == true) {
            this->initialize_serialized_coordinate();
        }
        std::printf("Parameters\n");
        std::printf("Limit: x: %f, y: %f, z: %f\n", x_limit, y_limit, z_limit);
        std::printf("dx: %f\n", this->dx);
        std::printf("Rho0: %f\n", rho0);
        std::printf("Mass: %f\n", mass);
        std::printf("Stiffness: %f\n", stiffness);
        std::printf("kernel radius: %f\n", h);
        std::printf("epsilon: %f\n", this->epsilon);
        std::printf("visc: %f\n", this->visc);

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
        float hsq = std::pow(this->h, 2);
        #pragma omp parallel for
        for(size_t i = 0; i < this->particles.size(); i++) {
            Particle &pi = this->particles[i];
            pi.rho = 0.0;
            for(auto &pj : this->particles) {
                Eigen::Vector3d rij = pj.r - pi.r;
                float r2 = rij.squaredNorm();

                if (r2 < hsq) {
                    pi.rho += this->mass * this->poly6 * pow(hsq-r2, 3.f);
                }
            }
            //pi.p = std::max(0.f, this->stiffness * (pi.rho - this->rho0));
            pi.p = this->stiffness * (pi.rho - this->rho0);
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
                if (r < this->h) {
                    //fpress += -rij.normalized() * this->mass * (pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD * pow(this->h - r, 3.f);
                    //fpress += -rij.normalized() * this->mass * (pi.p + pj.p) / (2.f * pj.rho) * this->spiky_grad * pow(this->h - r, 2.f);
                    fpress += rij.normalized() * this->mass * (pi.p + pj.p) / (2.f * pj.rho) * this->spiky_grad * pow(this->h - r, 2.f);
                    fvisc += this->visc * this->mass * (pj.v - pi.v) / pj.rho * this->visc_lap * (this->h -r);
                }
            }
            //Eigen::Vector3d fgrav = G*this->mass / pi.rho;
            Eigen::Vector3d fgrav = G * pi.rho;
            pi.f = fpress + fvisc + fgrav;
        }
    }
    void integrate(void) {
        #pragma omp parallel for
        for(size_t i = 0; i < this->particles.size(); i++) {
            Particle &p = this->particles[i];
            p.v += this->dt * p.f / p.rho;
            p.r += this->dt * p.v;
            this->time += dt;
            if (p.r(0) - this->epsilon < 0.f){
                p.v(0) *= BOUND_DAMPING;
                p.r(0) = this->epsilon;
            }
            if (p.r(0) + this->epsilon > this->x_limit) {
                p.v(0) *= BOUND_DAMPING;
                p.r(0) = this->x_limit - this->epsilon;
            }
            if (p.r(1) - this->epsilon < 0.f){
                p.v(1) *= BOUND_DAMPING;
                p.r(1) = this->epsilon;
            }
            if (p.r(1) + this->epsilon > this->y_limit) {
                p.v(1) *= BOUND_DAMPING;
                p.r(1) = this->y_limit - this->epsilon;
            }

            if (p.r(2) - this->epsilon < 0.f){
                p.v(2) *= BOUND_DAMPING;
                p.r(2) = this->epsilon;
            }
            if (p.r(2) + this->epsilon > this->z_limit) {
                p.v(2) *= BOUND_DAMPING;
                p.r(2) = this->z_limit - this->epsilon;
            }
        }
    }
    size_t get_num_particles() {
        return this->particles.size();
    }
    void initialize_serialized_coordinate(void) {
        this->serialized_coordinates.clear();
        for(const auto &ptcl: this->particles) {
            this->serialized_coordinates.push_back(float(ptcl.r[0]));
            this->serialized_coordinates.push_back(float(ptcl.r[1]));
            this->serialized_coordinates.push_back(float(ptcl.r[2]));
        }
    }
    void update_serialized_coordinate(void) {
        this->serialized_coordinates.clear();
        for(size_t i = 0; i < this->particles.size(); i++) {
            this->serialized_coordinates[3*i + 0] = this->particles[i].r[0];
            this->serialized_coordinates[3*i + 1] = this->particles[i].r[1];
            this->serialized_coordinates[3*i + 2] = this->particles[i].r[2];
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

    float dt;
    float time;
    float x_limit, y_limit, z_limit;
    float rho0;
    float mass;
    float stiffness;
    float h;
    float dx;
    float epsilon;
    float visc;

    float poly6, spiky_grad, visc_lap;
    void define_kernel_constants(void) {
        this->poly6 = 315.f / (64.f * M_PI * pow(this->h, 9.f));
        this->spiky_grad = -45.f / (M_PI * pow(this->h, 6.f));
        this->visc_lap = 45.f / (M_PI * pow(this->h, 6.f));
    }
};

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(my_module) {

    emscripten::class_<SPHSimulator>("SPHSimulator")
        .constructor<unsigned long>()
        .function("init_simulation", &SPHSimulator::init_simulation)
        .function("step", &SPHSimulator::step)
        .function("get_num_particles", &SPHSimulator::get_num_particles)
        .function("get_serialized_coordinate_buffer_ptr", &SPHSimulator::get_serialized_coordinate_buffer_ptr)
        .function("get_serialized_coordinate_buffer_size", &SPHSimulator::get_serialized_coordinate_buffer_size)
        .property("x_limit", &SPHSimulator::x_limit)
        .property("y_limit", &SPHSimulator::y_limit)
        .property("z_limit", &SPHSimulator::z_limit)
        .property("rho0", &SPHSimulator::rho0)
        .property("mass", &SPHSimulator::mass)
        .property("stiffness", &SPHSimulator::stiffness)
        .property("h", &SPHSimulator::h)
        .property("dx", &SPHSimulator::dx)
        .property("dt", &SPHSimulator::dt)
        .property("epsilon", &SPHSimulator::epsilon)
        .property("visc", &SPHSimulator::visc)
        ;
}
#else

int main(void) {

    SPHSimulator sim(2000);
    sim.init_simulation();

    for(int j = 0; j < 10; j++) {
        for(int i = 0; i < 10; i++) {
            sim.compute_density_pressure();
            sim.compute_forces();
            sim.integrate();

        }

    }
    return 0;
}
#endif
