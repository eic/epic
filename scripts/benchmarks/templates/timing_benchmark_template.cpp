#include <DD4hep/Detector.h>
#include <DD4hep/FieldTypes.h>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace dd4hep;
using namespace std;

int main() {
  {
    Detector& detector = Detector::getInstance();
    detector.fromXML("{xml_file}");

    auto field = detector.field();
    if (!field.isValid()) {
      {
        cerr << "ERROR: No field found in detector description" << endl;
        return 1;
      }
    }

    // Generate random test points
    random_device rd;
    mt19937 gen(42); // Fixed seed for reproducibility
    uniform_real_distribution<> r_dist({r_min}, {r_max});
    uniform_real_distribution<> phi_dist(0, 2 * M_PI);
    uniform_real_distribution<> z_dist({z_min}, {z_max});

    vector<tuple<double, double, double>> test_points;
    test_points.reserve({n_points});

    for (int i = 0; i < {n_points}; ++i) {
      {
        double r   = r_dist(gen);
        double phi = phi_dist(gen);
        double z   = z_dist(gen);
        double x   = r * cos(phi);
        double y   = r * sin(phi);
        test_points.emplace_back(x, y, z);
      }
    }

    // Warm up
    double pos[3], field_val[3];
    for (int i = 0; i < 1000; ++i) {
      {
        auto [x, y, z] = test_points[i % test_points.size()];
        pos[0]         = x;
        pos[1]         = y;
        pos[2]         = z;
        field.magneticField(pos, field_val);
      }
    }

    // Timing test
    auto start = chrono::high_resolution_clock::now();

    double sum_bx = 0, sum_by = 0, sum_bz = 0;
    for (const auto& point : test_points) {
      {
        auto [x, y, z] = point;
        pos[0]         = x;
        pos[1]         = y;
        pos[2]         = z;
        field.magneticField(pos, field_val);
        sum_bx += field_val[0];
        sum_by += field_val[1];
        sum_bz += field_val[2];
      }
    }

    auto end      = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

    // Output results
    cout << "{{" << endl;
    cout << "  \\" n_points\\": " << {n_points} << "," << endl;
    cout << "  \\" total_time_us\\": " << duration.count() << "," << endl;
    cout << "  \\" time_per_evaluation_ns\\": " << (duration.count() * 1000.0 / {n_points}) << ","
         << endl;
    cout << "  \\" evaluations_per_second\\": " << ({n_points} * 1e6 / duration.count()) << ","
         << endl;
    cout << "  \\" sum_field\\": [" << sum_bx << ", " << sum_by << ", " << sum_bz << "]," << endl;
    cout << "  \\" field_magnitude_avg\\": "
         << sqrt(sum_bx * sum_bx + sum_by * sum_by + sum_bz * sum_bz) / {n_points} << endl;
    cout << "}}" << endl;

    return 0;
  }
}
