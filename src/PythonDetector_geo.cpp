// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025 ePIC Collaboration
//
// C++ trampoline plugin that delegates geometry construction to a Python function.
//
// Usage in compact XML:
//   <detector id="X_ID" name="X" type="epic_PythonDetector"
//             module="mydetector_geo" function="create_detector"
//             readout="XHits">
//     ...any sub-elements your Python function reads...
//   </detector>
//
// The Python function must have the signature:
//   def create_detector(description, xml_element, sens) -> dd4hep.DetElement
//
// Python module search path:
//   Python uses its normal module search path (for example, as configured by
//   sys.path / $PYTHONPATH). Optionally add a "pythonpath" attribute to the
//   <detector> element to prepend an additional search directory.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "TPython.h"

#include <cstdint>
#include <sstream>
#include <string>

using namespace dd4hep;

// Global DetElement that receives the result from Python via SwapWithObjAtAddr.
// Using a module-level static means re-entrant calls are not supported, but that
// matches the single-threaded geometry-build model of DD4hep.
static dd4hep::DetElement _epic_python_det_result;

// Expand ${VAR} environment variable references in a string path.
static std::string expand_env(const std::string& input) {
  std::string out;
  for (std::size_t i = 0; i < input.size();) {
    if (input[i] == '$' && i + 1 < input.size() && input[i + 1] == '{') {
      std::size_t end = input.find('}', i + 2);
      if (end != std::string::npos) {
        std::string var = input.substr(i + 2, end - i - 2);
        const char* val = std::getenv(var.c_str());
        if (val)
          out += val;
        else
          out += input.substr(i, end - i + 1);
        i = end + 1;
        continue;
      }
    }
    out += input[i++];
  }
  return out;
}

static Ref_t create_python_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t x_det  = e;
  std::string mod  = x_det.attr<std::string>(dd4hep::xml::Strng_t("module"));
  std::string func = x_det.attr<std::string>(dd4hep::xml::Strng_t("function"));

  // Pass addresses of C++ objects so Python can reconstruct cppyy proxies.
  auto desc_addr   = reinterpret_cast<std::intptr_t>(&description);
  auto xml_addr    = reinterpret_cast<std::intptr_t>(&e);
  auto sens_addr   = reinterpret_cast<std::intptr_t>(&sens);
  auto result_addr = reinterpret_cast<std::intptr_t>(&_epic_python_det_result);

  std::string pypath;
  if (x_det.hasAttr(dd4hep::xml::Strng_t("pythonpath"))) {
    pypath = expand_env(x_det.attr<std::string>(dd4hep::xml::Strng_t("pythonpath")));
  }

  std::ostringstream script;
  script << "import importlib.util, sys, cppyy\n";
  if (!pypath.empty()) {
    script << "if '" << pypath << "' not in sys.path:\n"
           << "    sys.path.insert(0, '" << pypath << "')\n";
  }
  script << "_m    = importlib.import_module('" << mod << "')\n"
         << "_desc = cppyy.bind_object(" << desc_addr << ", 'dd4hep::Detector')\n"
         << "_xml  = cppyy.bind_object(" << xml_addr << ", 'dd4hep::xml::Handle_t')\n"
         << "_sens = cppyy.bind_object(" << sens_addr << ", 'dd4hep::SensitiveDetector')\n"
         << "try:\n"
         << "    _res = _m." << func << "(_desc, _xml, _sens)\n"
         << "except Exception:\n"
         << "    import traceback; traceback.print_exc()\n"
         << "    raise\n"
         // Swap the result into the global C++ DetElement so C++ can retrieve it.
         << "cppyy.gbl.ROOT.Internal.SwapWithObjAtAddr['dd4hep::DetElement']("
         << "_res, " << result_addr << ")\n";

  printout(DEBUG, "PythonDetector", "Calling %s.%s()", mod.c_str(), func.c_str());
  if (!TPython::Exec(script.str().c_str())) {
    except("PythonDetector", "Python call to %s.%s() failed", mod.c_str(), func.c_str());
  }

  return _epic_python_det_result;
}

DECLARE_DETELEMENT(epic_PythonDetector, create_python_detector)
