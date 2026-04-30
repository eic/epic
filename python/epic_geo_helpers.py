# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2025 ePIC Collaboration
#
# Shared helpers for Python-based DD4hep geometry plugins.
#
# Usage in a plugin module:
#   from epic_geo_helpers import *
#
# Provides:
#   - DD4hep header / cppyy initialisation (idempotent)
#   - C++ helper namespace `epic_python` with hasAttr / getAttrDbl / getAttrBool / collValid
#   - Python wrapper functions: _handle, _collValid, _hasAttr, _getAttrBool
#   - Type aliases for the most-used DD4hep geometry types

import cppyy
from ROOT import gInterpreter, TMath  # noqa: F401  (re-exported via __all__)

# ---------------------------------------------------------------------------
# Load DD4hep headers into cppyy (idempotent – safe to call multiple times)
# ---------------------------------------------------------------------------
gInterpreter.ProcessLine('#include "DD4hep/DetFactoryHelper.h"')

# ---------------------------------------------------------------------------
# Thin C++ helpers for DD4hep XML attribute access.
#
# Background: dd4hep::xml::Element::hasAttr / attr<T> require a const XmlChar*
# argument (XmlChar = char16_t with the Xerces-C backend).  cppyy cannot
# auto-convert a Python str to char16_t*, so direct calls raise TypeError.
# The helpers below accept const char* and convert via dd4hep::xml::Strng_t.
#
# xml_coll_t::operator bool() is also not mapped to Python __bool__ by cppyy,
# so a separate collValid() helper is provided.
# ---------------------------------------------------------------------------
cppyy.cppdef(r"""
#include "DD4hep/DetFactoryHelper.h"
namespace epic_python {
  inline bool hasAttr(dd4hep::xml::Element el, const char* name) {
    return el.hasAttr(dd4hep::xml::Strng_t(name));
  }
  inline std::string getAttrStr(dd4hep::xml::Element el, const char* name) {
    return el.attr<std::string>(dd4hep::xml::Strng_t(name));
  }
  inline double getAttrDbl(dd4hep::xml::Element el, const char* name) {
    return el.attr<double>(dd4hep::xml::Strng_t(name));
  }
  inline bool getAttrBool(dd4hep::xml::Element el, const char* name) {
    return el.attr<bool>(dd4hep::xml::Strng_t(name));
  }
  // xml_coll_t::operator bool() is NOT mapped to Python __bool__; use this.
  inline bool collValid(const dd4hep::xml::Collection_t& c) {
    return static_cast<bool>(c);
  }
}
""")


# ---------------------------------------------------------------------------
# Python wrapper functions
# ---------------------------------------------------------------------------

def _handle(el):
    """Cast any xml Element type to Handle_t (required by xml_coll_t constructor).

    cppyy does not auto-upcast xml_det_t / xml_comp_t to Handle_t even though
    C++ would do so implicitly.
    """
    return cppyy.gbl.dd4hep.xml.Handle_t(el.ptr())


def _collValid(c) -> bool:
    """Return True if the xml_coll_t iterator is still valid.

    Use instead of ``while c:`` — cppyy maps Python __bool__ to the default
    object truthiness, not to C++ operator bool(), so ``while c:`` loops
    forever and eventually segfaults on a NULL pointer.
    """
    return cppyy.gbl.epic_python.collValid(c)


def _hasAttr(el, name: str) -> bool:
    """Return True if the XML element has the named attribute."""
    return cppyy.gbl.epic_python.hasAttr(el, name)


def _getAttrBool(el, name: str) -> bool:
    """Return the named XML attribute as a bool."""
    return cppyy.gbl.epic_python.getAttrBool(el, name)


# ---------------------------------------------------------------------------
# DD4hep type aliases
# ---------------------------------------------------------------------------
_dd4hep      = cppyy.gbl.dd4hep
Assembly     = _dd4hep.Assembly
DetElement   = _dd4hep.DetElement
PlacedVolume = _dd4hep.PlacedVolume
Position     = _dd4hep.Position
Transform3D  = _dd4hep.Transform3D
RotationY    = cppyy.gbl.ROOT.Math.RotationY
Tube         = _dd4hep.Tube
Volume       = _dd4hep.Volume

_U        = _dd4hep.xml.Strng_t       # equivalent of the _U("tag") macro
_toString = _dd4hep.xml._toString     # equivalent of the _toString(n, "fmt") macro

xml_det_t  = cppyy.gbl.xml_det_t
xml_coll_t = cppyy.gbl.xml_coll_t
xml_comp_t = cppyy.gbl.xml_comp_t

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
__all__ = [
    # ROOT helper
    "TMath",
    # Python wrapper functions
    "_handle", "_collValid", "_hasAttr", "_getAttrBool",
    # DD4hep geometry types
    "Assembly", "DetElement", "PlacedVolume", "Position",
    "Transform3D", "RotationY", "Tube", "Volume",
    # XML helpers
    "_U", "_toString", "xml_det_t", "xml_coll_t", "xml_comp_t",
]
