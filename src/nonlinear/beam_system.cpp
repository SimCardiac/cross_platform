#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>

// Constants
const double m = 0.1;          // Mass of each ball (kg)
const double k_spring = 1000.0; // Spring constant (N/m)
const double L0 = 0.2;         // Rest length (m)
const double g = 9.81;         // Gravity (m/s^2)
const double T_total = 50.0;    // Total simulation time (s)
const double dt = 0.001;       // Time step (s)

constexpr int N = 500;


struct Segment {
    int x, y;
};

struct Vector2 {
    double x, y;
    Vector2 operator+(const Vector2& other) const { return {x + other.x, y + other.y}; }
    Vector2 operator-(const Vector2& other) const { return {x - other.x, y - other.y}; }
    Vector2 operator*(double scalar) const { return {x * scalar, y * scalar}; }
    Vector2 operator/(double scalar) const { return {x / scalar, y / scalar}; }
    double norm() const { return std::sqrt(x*x + y*y); }
};

struct Ball {
    Vector2 pos;
    Vector2 pos_current;
    Vector2 vel;
    Vector2 force;
};

void save_vtp(const std::string& filename, const std::vector<Ball>& balls, const std::vector<Segment>& segments) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Ensure dot is used as decimal separator and use scientific notation
    file.imbue(std::locale::classic());
    file << std::fixed << std::setprecision(10);

    file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <PolyData>\n";
    file << "    <Piece NumberOfPoints=\"" << balls.size() << "\" NumberOfLines=\"" << segments.size() << "\">\n";
    
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& ball : balls) {
        file << ball.pos.x << " " << ball.pos.y << " " << 0.0 << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    file << "      <Lines>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& seg : segments) {
        file << seg.x << " " << seg.y << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 0; i < segments.size(); ++i) {
        file << (i + 1) * 2 << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Lines>\n";

    file << "    </Piece>\n";
    file << "  </PolyData>\n";
    file << "</VTKFile>\n";
    file.close();
}

void save_pvd(const std::string& filename, const std::vector<std::pair<double, std::string>>& steps) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <Collection>\n";
    for (const auto& step : steps) {
        file << "    <DataSet timestep=\"" << step.first << "\" group=\"\" part=\"0\" file=\"" << step.second << "\"/>\n";
    }
    file << "  </Collection>\n";
    file << "</VTKFile>\n";
    file.close();
}

double calculate_total_energy(const std::vector<Ball>& balls, const std::vector<Segment>& segments, const std::vector<double>& initial_length) {
    double kinetic_energy = 0.0;
    double elastic_energy = 0.0;
    double gravitational_energy = 0.0;

    // Kinetic & Gravitational
    for (const auto& ball : balls) {
        double v2 = ball.vel.x * ball.vel.x + ball.vel.y * ball.vel.y;
        kinetic_energy += 0.5 * m * v2;
        gravitational_energy -= m * g * ball.pos.y; // y-axis points down, so potential decreases with y
    }

    // Elastic
    for (size_t i = 0; i < segments.size(); ++i) {
        Vector2 d = balls[segments[i].y].pos - balls[segments[i].x].pos;
        double len = d.norm();
        double stretch = len - initial_length[i];
        elastic_energy += 0.5 * k_spring * stretch * stretch;
    }

    return kinetic_energy + elastic_energy + gravitational_energy;
}

int main() {
    // int N = 5; // Removed shadowing variable
    double v0 = 1.0;
    std::vector<int> fixed_points(N);
    std::vector<Ball> balls(N);
    std::vector<Segment> segments(N - 1);
    std::vector<double> initial_length(N - 1);
    std::vector<double> current_length(N - 1);
    
    // Explicit initialization for Ball 0
    for (int i = 1; i < N; ++i) {
        balls[i].pos = {0.0, i * L0};
        balls[i].vel = {0.0, 0.0};
        balls[i].force = {0.0, 0.0};
        fixed_points[i] = 0; // Not fixed
    }
    fixed_points[0] = 1; // Fixed point at the top
    balls[N-1].vel = {1.0, 0.0};
    for (int i = 1; i < N; ++i) {
        segments[i-1] = {i-1, i};
        initial_length[i-1] = L0;
        current_length[i-1] = L0;
        printf("Segment %d: Ball %d to Ball %d\n", i-1, segments[i-1].x, segments[i-1].y);
    }

    int steps = static_cast<int>(T_total / dt);
    std::vector<std::pair<double, std::string>> pvd_steps;
    
    // Save initial state
    std::string filename = "beam_system_00000.vtp";
    save_vtp("build/" + filename, balls, segments);
    pvd_steps.push_back({0.0, filename});

    for (int step = 1; step <= steps; ++step) {
        for (int i = 0; i < N; ++i) {
            balls[i].force = {0.0, m * g}; // Gravity (y is down)
        }

        // Spring Forces
        // iterate over segments
        for (int i = 0; i < N - 1; ++i) {
            int idx1 = segments[i].x;;
            int idx2 = segments[i].y;
            Vector2 dir = balls[idx2].pos - balls[idx1].pos;
            double len = dir.norm();
            current_length[i] = len;
            Vector2 dir_normalized = dir / len;
            double fspring = k_spring * (len - initial_length[i]); // Corrected sign: positive for tension
            Vector2 force_vec = dir_normalized * fspring;
            balls[idx1].force = balls[idx1].force + force_vec;
            balls[idx2].force = balls[idx2].force - force_vec;
        }

        // 2. Update Position and Velocity (Symplectic Euler)
        for (int i = 0; i < N; ++i) {
            if (fixed_points[i]) continue; // Skip fixed points

            // v_new = v_old + a_old * dt
            Vector2 accel = balls[i].force / m;
            balls[i].vel = balls[i].vel + accel * dt;

            // x_new = x_old + v_new * dt
            balls[i].pos = balls[i].pos + balls[i].vel * dt;
        }
        
        if (step % 50 == 0) { // Save every 50 steps (0.05s)
            std::stringstream ss;
            ss << "beam_system_" << std::setw(5) << std::setfill('0') << step << ".vtp";
            std::string filename = ss.str();
            save_vtp("build/" + filename, balls, segments);
            pvd_steps.push_back({step * dt, filename});
            
            if (step % 50 == 0) {
                double energy = calculate_total_energy(balls, segments, initial_length);
                printf("Step %d: Energy = %.6e, Saved %s\n", step, energy, ("build/" + filename).c_str());
            }
        }
    }
    
    save_pvd("build/beam_system.pvd", pvd_steps);
    std::cout << "Saved PVD file to build/beam_system.pvd" << std::endl;

    return 0;
}
