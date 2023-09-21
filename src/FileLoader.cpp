// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Factories.h>
#include <DD4hep/Primitives.h>
#include <DD4hep/Printout.h>
#include <XML/Utilities.h>

#include <fmt/core.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include "FileLoaderHelper.h"

using namespace dd4hep;

void usage(int argc, char** argv)
{
  std::cerr << "Usage: -plugin <name> -arg [-arg]                                                  \n"
               "     cache:<string>           cache location (may be read-only)                    \n"
               "     file:<string>            file location                                        \n"
               "     url:<string>             url location                                         \n"
               "     cmd:<string>             download command with {0} for url, {1} for output    \n"
               "\tArguments given: "
            << arguments(argc, argv) << std::endl;
  std::exit(EINVAL);
}

// Plugin to download files
long load_file(Detector& /* desc */, int argc, char** argv)
{
  // argument parsing
  std::string cache, file, url;
  std::string cmd("curl --retry 5 -f {0} -o {1}");
  for (int i = 0; i < argc && argv[i]; ++i) {
    if (0 == std::strncmp("cache:", argv[i], 6))
      cache = (argv[i] + 6);
    else if (0 == std::strncmp("file:", argv[i], 5))
      file = (argv[i] + 5);
    else if (0 == std::strncmp("url:", argv[i], 4))
      url = (argv[i] + 4);
    else if (0 == std::strncmp("cmd:", argv[i], 4))
      cmd = (argv[i] + 4);
    else {
      std::cerr << "Unexpected argument \"" << argv[i] << "\"" << std::endl;
      usage(argc, argv);
    }
  }
  printout(DEBUG, "FileLoader", "arg cache: " + cache);
  printout(DEBUG, "FileLoader", "arg file: " + file);
  printout(DEBUG, "FileLoader", "arg url: " + url);
  printout(DEBUG, "FileLoader", "arg cmd: " + cmd);

  // if file or url is empty, do nothing
  if (file.empty()) {
    printout(WARNING, "FileLoader", "no file specified");
  }
  if (url.empty()) {
    printout(WARNING, "FileLoader", "no url specified");
  }

  EnsureFileFromURLExists(url, file, cache, cmd);

  return 1;
}

DECLARE_APPLY(epic_FileLoader, load_file)
