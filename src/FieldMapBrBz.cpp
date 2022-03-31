#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/FieldTypes.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
namespace fs = std::filesystem;

#include "FileLoaderHelper.h"

using namespace dd4hep;


// implementation of the field map
class FieldMapBrBz : public dd4hep::CartesianField::Object
{
public:
    FieldMapBrBz(const std::string &field_type = "magnetic");
    void Configure(double rmin, double rmax, double rstep, double zmin, double zmax, double zstep);
    void LoadMap(const std::string &map_file, double scale);
    void GetIndices(double r, double z, int &ir, int &iz, double &dr, double &dz);
    void SetTransform(const Transform3D &tr) { trans = tr; trans_inv = tr.Inverse(); }

    virtual void fieldComponents(const double *pos, double *field);


private:
    Transform3D trans, trans_inv;
    double rmin, rmax, rstep, zmin, zmax, zstep;
    std::vector<std::vector<std::vector<double>>> Bvals;
};

// constructor
FieldMapBrBz::FieldMapBrBz(const std::string &field_type)
{
    std::string ftype = field_type;
    for (auto &c : ftype) { c = tolower(c); }

    // set type
    if (ftype == "magnetic") {
        type = CartesianField::MAGNETIC;
    } else if (ftype == "electric") {
        type = CartesianField::ELECTRIC;
    } else {
        type = CartesianField::UNKNOWN;
        std::cout << "FieldMapBrBz Warning: Unknown field type " << field_type << "!" << std::endl;
    }
}

void FieldMapBrBz::Configure(double r1, double r2, double rs, double z1, double z2, double zs)
{
    rmin = r1;
    rmax = r2;
    rstep = rs;
    zmin = z1;
    zmax = z2;
    zstep = zs;

    int nr = int((r2 - r1)/rs) + 2;
    int nz = int((z2 - z1)/zs) + 2;

    Bvals.resize(nr);
    for (auto &B2 : Bvals) {
        B2.resize(nz);
        for (auto &B : B2) {
            B.resize(2, 0.);
        }
    }
}

void FieldMapBrBz::GetIndices(double r, double z, int &ir, int &iz, double &dr, double &dz)
{
    // boundary check
    if (r > rmax || r < rmin || z > zmax || z < zmin) {
        ir = -1;
        iz = -1;
        return;
    }
    // get indices
    double idr, idz;
    dr = std::modf((r - rmin)/rstep, &idr);
    dz = std::modf((z - zmin)/zstep, &idz);

    ir = static_cast<int>(idr);
    iz = static_cast<int>(idz);
}

// load data
void FieldMapBrBz::LoadMap(const std::string &map_file, double scale)
{
    std::string line;
    std::ifstream input(map_file);
    if (!input) {
        std::cout << "FieldMapBrBz Error: file \"" << map_file << "\" cannot be read." << std::endl;
    }

    double r, z, br, bz;
    int ir, iz;
    double dr, dz;
    while (std::getline(input, line).good()) {
        std::istringstream iss(line);
        iss >> r >> z >> br >> bz;
        GetIndices(r, z, ir, iz, dr, dz);
        if (ir < 0 || iz < 0) {
            std::cout << "FieldMapBrBz Warning: coordinates out of range ("
                      << r << ", " << z << "), skipped it." << std::endl;
        } else {
            Bvals[ir][iz] = {br*scale, bz*scale};
            // ROOT::Math::XYZPoint p(r, 0, z);
            // std::cout << p << " -> " << trans*p << std::endl;
            // std::cout << ir << ", " << iz << ", " << br << ", " << bz << std::endl;
        }
    }
}

// get field components
void FieldMapBrBz::fieldComponents(const double *pos, double *field)
{
    // coordinate conversion
    auto p = trans_inv*ROOT::Math::XYZPoint(pos[0], pos[1], pos[2]);

    // coordinates conversion
    const double r = sqrt(p.x()*p.x() + p.y()*p.y());
    const double z = p.z();
    const double phi = atan2(p.y(), p.x());

    int ir, iz;
    double dr, dz;
    GetIndices(r, z, ir, iz, dr, dz);

    // out of the range
    if (ir < 0 || iz < 0) { return; }

    // p1    p3
    //    p
    // p0    p2
    auto &p0 = Bvals[ir][iz];
    auto &p1 = Bvals[ir][iz + 1];
    auto &p2 = Bvals[ir + 1][iz];
    auto &p3 = Bvals[ir + 1][iz + 1];

    // linear interpolation
    double Br = p0[0] * (1-dr) * (1-dz)
              + p1[0] * (1-dr) *    dz
              + p2[0] *    dr  * (1-dz)
              + p3[0] *    dr  *    dz;

    double Bz = p0[1] * (1-dr) * (1-dz)
              + p1[1] * (1-dr) *    dz
              + p2[1] *    dr  * (1-dz)
              + p3[1] *    dr  *    dz;

    // convert Br Bz to Bx By Bz
    auto B = trans*ROOT::Math::XYZPoint(Br*sin(phi), Br*cos(phi), Bz);
    field[0] += B.x()*tesla;
    field[1] += B.y()*tesla;
    field[2] += B.z()*tesla;
    return;
}


// assign the field map to CartesianField
static Ref_t create_field_map_brbz(Detector & /*lcdd*/, xml::Handle_t handle)
{
    xml_comp_t x_par(handle);

    if (!x_par.hasAttr(_Unicode(field_map))) {
        throw std::runtime_error("FieldMapBrBz Error: must have an xml attribute \"field_map\" for the field map.");
    }

    CartesianField field;
    std::string field_type = x_par.attr<std::string>(_Unicode(field_type));

    // dimensions
    xml_comp_t x_dim = x_par.dimensions();

    // min, max, step
    xml_comp_t r_dim = x_dim.child(_Unicode(transverse));
    xml_comp_t z_dim = x_dim.child(_Unicode(longitudinal));

    std::string field_map_file = x_par.attr<std::string>(_Unicode(field_map));
    std::string field_map_url = x_par.attr<std::string>(_Unicode(url));

    EnsureFileFromURLExists(field_map_url,field_map_file);

    double field_map_scale = x_par.attr<double>(_Unicode(scale));

    if( !fs::exists(fs::path(field_map_file))  ) {
        printout(ERROR, "FieldMapBrBz", "file " + field_map_file + " does not exist");
        printout(ERROR, "FieldMapBrBz", "use a FileLoader plugin before the field element");
        std::quick_exit(1);
    }

    auto map = new FieldMapBrBz(field_type);
    map->Configure(r_dim.rmin(), r_dim.rmax(), r_dim.step(), z_dim.zmin(), z_dim.zmax(), z_dim.step());

    // translation, rotation
    static double deg2r = ROOT::Math::Pi()/180.;
    RotationZYX rot(0., 0., 0.);
    if (x_dim.hasChild(_Unicode(rotation))) {
        xml_comp_t rot_dim = x_dim.child(_Unicode(rotation));
        rot = RotationZYX(rot_dim.z()*deg2r, rot_dim.y()*deg2r, rot_dim.x()*deg2r);
    }

    Translation3D trans(0., 0., 0.);
    if (x_dim.hasChild(_Unicode(translation))) {
        xml_comp_t trans_dim = x_dim.child(_Unicode(translation));
        trans = Translation3D(trans_dim.x(), trans_dim.y(), trans_dim.z());
    }
    map->SetTransform(trans*rot);

    map->LoadMap(field_map_file, field_map_scale);
    field.assign(map, x_par.nameStr(), "FieldMapBrBz");

    return field;
}

DECLARE_XMLELEMENT(FieldMapBrBz, create_field_map_brbz)

