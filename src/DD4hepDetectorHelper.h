// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sakib Rahman

#pragma once

#include <DD4hep/DetElement.h>
#include <DDRec/DetectorData.h>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include "DD4hep/DetFactoryHelper.h"


namespace DD4hepDetectorHelper {

template <typename T>
T& ensureExtension(dd4hep::DetElement& elt) {
  T* ext = elt.extension<T>(false);
  if (ext == nullptr) {
    ext = new T();
  }
  elt.addExtension<T>(ext);
  return *ext;
}

inline void xmlToProtoSurfaceMaterial(const xml_comp_t& x_material,
                                      dd4hep::rec::VariantParameters& params,
                                      const std::string& baseTag) {
  using namespace std::string_literals;
  using namespace dd4hep;

  // Add the layer material flag
  params.set(baseTag, true);
  // prepare everything here
  std::string mSurface = x_material.attr<std::string>("surface");
  std::string mBinning = x_material.attr<std::string>("binning");
  boost::char_separator<char> sep(",");
  boost::tokenizer binTokens(mBinning, sep);
  const auto n = std::distance(binTokens.begin(), binTokens.end());
  if (n == 2) {
    // Fill the bins
    auto bin = binTokens.begin();
    std::string bin0 = *(bin);
    std::string bin1 = *(++bin);
    size_t nBins0 = x_material.attr<int>("bins0");
    size_t nBins1 = x_material.attr<int>("bins1");
    // Add the material tags
    std::string btmSurface = baseTag + "_"s + mSurface;
    params.set<bool>(btmSurface, true);
    params.set<int>(btmSurface + "_"s + bin0, nBins0);
    params.set<int>(btmSurface + "_"s + bin1, nBins1);
  }
}

}  // namespace DD4hepDetectorHelper

