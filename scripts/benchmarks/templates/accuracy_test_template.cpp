#include <DD4hep/Detector.h>
#include <DD4hep/FieldTypes.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace dd4hep;

int main() {
  {
    Detector& detector = Detector::getInstance();
    detector.fromXML("{xml_file}");

    auto field = detector.field();
    if (!field.isValid()) {
      {
        std::cerr << "ERROR: No field found" << std::endl;
        return 1;
      }
    }

    // Test field properties at key points
    double pos[3], field_val[3];

    // Test at origin
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;
    field.magneticField(pos, field_val);
    double field_at_origin = sqrt(field_val[0] * field_val[0] + field_val[1] * field_val[1] +
                                  field_val[2] * field_val[2]);

    // Test cylindrical symmetry (for BrBz fields)
    double asymmetry = 0.0;
    for (int phi_deg = 0; phi_deg < 360; phi_deg += 45) {
      {
        double phi = phi_deg * M_PI / 180.0;
        double r   = 50.0; // 50 cm
        pos[0]     = r * cos(phi);
        pos[1]     = r * sin(phi);
        pos[2]     = 0;
        field.magneticField(pos, field_val);
        double field_mag = sqrt(field_val[0] * field_val[0] + field_val[1] * field_val[1] +
                                field_val[2] * field_val[2]);
        asymmetry += abs(field_mag - field_at_origin) / field_at_origin;
      }
    }
    asymmetry /= 8.0; // Average over 8 points

    // Test field gradient
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;
    field.magneticField(pos, field_val);
    double field_center = field_val[2]; // Bz at center

    pos[2] = 10.0; // 10 cm offset
    field.magneticField(pos, field_val);
    double field_offset = field_val[2];
    double gradient     = (field_offset - field_center) / 10.0; // T/cm

    std::cout << "{{" << std::endl;
    std::cout << "  \\" field_at_origin_T\\": " << field_at_origin << "," << std::endl;
    std::cout << "  \\" cylindrical_asymmetry\\": " << asymmetry << "," << std::endl;
    std::cout << "  \\" field_gradient_T_per_cm\\": " << gradient << std::endl;
    std::cout << "}}" << std::endl;

    return 0;
  }
}
