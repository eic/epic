#pragma once

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Factories.h>
#include <DD4hep/Primitives.h>
#include <DD4hep/Printout.h>

#include <fmt/core.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <regex>
#include <string>

namespace fs = std::filesystem;

using namespace dd4hep;

// Function to download files
inline void EnsureFileFromURLExists(std::string url, std::string file, std::string cache_str = "",
                                    std::string cmd = "curl --silent --retry 5 --location --fail {0} --output {1}")
{
  // parse cache for environment variables
  auto pos = std::string::npos;
  while ((pos = cache_str.find('$')) != std::string::npos) {
    auto after = cache_str.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                             "abcdefghijklmnopqrstuvwxyz"
                                             "0123456789"
                                             "_",
                                             pos + 1);
    if (after == std::string::npos)
      after = cache_str.size(); // cache ends on env var
    const std::string env_name(cache_str.substr(pos + 1, after - pos - 1));
    auto              env_ptr = std::getenv(env_name.c_str());
    const std::string env_value(env_ptr != nullptr ? env_ptr : "");
    cache_str.erase(pos, after - pos);
    cache_str.insert(pos, env_value);
    printout(INFO, "FileLoader", "$" + env_name + " -> " + env_value);
  }

  // tokenize cache on regex
  std::regex cache_sep(":");
  std::sregex_token_iterator cache_iter(cache_str.begin(), cache_str.end(), cache_sep, -1);
  std::sregex_token_iterator cache_end;
  std::vector<std::string> cache_vec(cache_iter, cache_end);

  // create file path
  fs::path file_path(file);

  // create hash from url, hex of unsigned long long
  std::string hash = fmt::format("{:016x}", dd4hep::detail::hash64(url)); // TODO: Use c++20 std::fmt

  // create file parent path, if not exists
  fs::path parent_path = file_path.parent_path();
  if (!fs::exists(parent_path)) {
    if (fs::create_directories(parent_path) == false) {
      printout(ERROR, "FileLoader", "parent path " + parent_path.string() + " cannot be created");
      printout(ERROR, "FileLoader", "hint: try running 'mkdir -p " + parent_path.string() + "'");
      std::_Exit(EXIT_FAILURE);
    }
  }

  // if file exists and is symlink to correct hash
  fs::path hash_path(parent_path / hash);
  if (fs::exists(file_path) && fs::equivalent(file_path, hash_path)) {
    printout(INFO, "FileLoader", "link " + file + " -> hash " + hash + " already exists");
    return;
  }

  // if hash does not exist, we try to retrieve file from cache
  if (!fs::exists(hash_path)) {
    // recursive loop into cache directories
    bool success = false;
    for (auto cache : cache_vec) {
      fs::path cache_path(cache);
      printout(INFO, "FileLoader", "cache " + cache_path.string());
      if (fs::exists(cache_path)) {
        for (auto const& dir_entry : fs::recursive_directory_iterator(cache_path)) {
          if (!dir_entry.is_directory())
            continue;
          fs::path cache_dir_path = cache_path / dir_entry;
          printout(INFO, "FileLoader", "checking " + cache_dir_path.string());
          fs::path cache_hash_path = cache_dir_path / hash;
          if (fs::exists(cache_hash_path)) {
            // symlink hash to cache/.../hash
            printout(INFO, "FileLoader", "file " + file + " with hash " + hash + " found in " + cache_hash_path.string());
            try {
              fs::create_symlink(cache_hash_path, hash_path);
              success = true;
            } catch (const fs::filesystem_error&) {
              printout(ERROR, "FileLoader",
                       "unable to link from " + hash_path.string() + " to " + cache_hash_path.string());
              printout(ERROR, "FileLoader", "hint: this may be resolved by removing directory " + parent_path.string());
              printout(ERROR, "FileLoader", "hint: or in that directory removing the file or link " + cache_hash_path.string());
              std::_Exit(EXIT_FAILURE);
            }
            break;
          }
        }
      }
      if (success) break;
    }
  }

  // if hash does not exist, we try to retrieve file from url
  if (!fs::exists(hash_path)) {
    cmd = fmt::format(cmd, url, hash_path.c_str()); // TODO: Use c++20 std::fmt
    printout(INFO, "FileLoader", "downloading " + file + " as hash " + hash + " with " + cmd);
    // run cmd
    auto ret = std::system(cmd.c_str());
    if (!fs::exists(hash_path)) {
      printout(ERROR, "FileLoader", "unable to run the download command " + cmd);
      printout(ERROR, "FileLoader", "the return value was ", ret);
      printout(ERROR, "FileLoader", "hint: check the command and try running manually");
      printout(ERROR, "FileLoader", "hint: allow insecure connections on some systems with the flag -k");
      std::_Exit(EXIT_FAILURE);
    }
  }

  // check if file already exists
  if (fs::exists(file_path)) {
    // file already exists
    if (fs::is_symlink(file_path)) {
      // file is symlink
      if (fs::equivalent(hash_path, fs::read_symlink(file_path))) {
        // link points to correct path
        return;
      } else {
        // link points to incorrect path
        if (fs::remove(file_path) == false) {
          printout(ERROR, "FileLoader", "unable to remove symlink " + file_path.string());
          printout(ERROR, "FileLoader", "we tried to create a symlink " + file_path.string() + " to the actual resource, " +
                                        "but a symlink already exists there and points to an incorrect location");
          printout(ERROR, "FileLoader", "hint: this may be resolved by removing directory " + parent_path.string());
          printout(ERROR, "FileLoader", "hint: or in that directory removing the file or link " + cache_hash_path.string());
          std::_Exit(EXIT_FAILURE);
        }
      }
    } else {
      // file exists but not symlink
      printout(ERROR, "FileLoader", "file " + file_path.string() + " already exists but is not a symlink");
      printout(ERROR, "FileLoader", "we tried to create a symlink " + file_path.string() + " to the actual resource, " +
                                    "but a file already exists there and we will not remove it automatically");
      printout(ERROR, "FileLoader", "hint: backup the file, remove it manually, and retry");
      std::_Exit(EXIT_FAILURE);
    }
  }
  // file_path now does not exist

  // symlink file_path to hash_path
  try {
    // use new path from hash so file link is local
    fs::create_symlink(fs::path(hash), file_path);
  } catch (const fs::filesystem_error&) {
    printout(ERROR, "FileLoader", "unable to link from " + file_path.string() + " to " + hash_path.string());
    printout(ERROR, "FileLoader", "check permissions and retry");
    std::_Exit(EXIT_FAILURE);
  }

  // final check of the file size
  if (fs::file_size(file_path) == 0) {
    printout(ERROR, "FileLoader", "zero file size of symlink from " + file_path.string() + " to " + hash_path.string());
    printout(ERROR, "FileLoader", "hint: check whether the file " + hash_path.string() + " has any content");
    printout(ERROR, "FileLoader", "hint: check whether the URL " + url + " has any content");
    std::_Exit(EXIT_FAILURE);
  }
}
