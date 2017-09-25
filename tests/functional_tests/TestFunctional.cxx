/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "TestFunctional.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <cstdio>
#include <thread>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <unistd.h>

#include <stdlib.h>
#include <src/support/utilities/SerializeDeformationState.h>
#include <src/launch/deformetrica.h>

namespace bpt = boost::property_tree;
namespace bfs = boost::filesystem;

namespace testing
{
namespace internal
{
enum GTestColor {
  COLOR_DEFAULT,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW
};

extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
}
}
#define PRINTF(...)  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__); } while(0)

// C++ stream interface
class TestCout : public std::stringstream
{
 public:
  ~TestCout()
  {
    PRINTF("%s",str().c_str());
  }
};

#define TEST_COUT  TestCout()


namespace def {
namespace test {

void TestFunctional::SetUp() {
  Test::SetUp();
}


#ifdef USE_CUDA
bool use_cuda = true;
#else
bool use_cuda = false;
#endif

#ifdef USE_DOUBLE_PRECISION
bool use_double_precision = true;
#else
bool use_double_precision = false;
#endif

TEST_F(TestFunctional, LoadAndRun) {

  auto file_exist = [](const std::string& f) {
    std::ifstream infile(f);
    return infile.good();
  };

  auto deformetrica_run = [](const std::string args, const std::string path, const std::string stdout_file, const std::string stderr_file) {
    pid_t pid;
    auto child = [&]() {
      std::vector<std::string> strs;
      boost::split(strs, args, boost::is_any_of("\t "));

      int argc = strs.size();
      char *argv[argc + 1];

      int i = 0;
      for (auto &str: strs)
        argv[i++] = (char *) str.c_str();

      /* CHANGING THE CURRENT WORKING DIRECTORY */
      bfs::current_path(path);

      /* Redirect std::cout streaming to empty/fake streaming */
      std::ofstream stdout(stdout_file);
      std::ofstream stderr(stderr_file);
      std::streambuf *coutbuf = std::cout.rdbuf();
      std::streambuf *cerrbuf = std::cerr.rdbuf();
      std::cout.rdbuf(stdout.rdbuf());
      std::cerr.rdbuf(stderr.rdbuf());

      /* Run deformetrica */
      deformetrica(argc, argv);

      /* Restore stdout/stderr to previous streaming */
      std::cout.rdbuf(coutbuf);
      std::cout.rdbuf(cerrbuf);
      exit(0);
    };

    auto father = [&]() {
      int returnStatus;
      waitpid(pid, &returnStatus, 0);
    };

    pid = fork();
    pid ? father() : child();
  };

  //check if '--replace-state-file' is present
  auto replace_state_file = std::find(args.begin(), args.end(), "--replace-state-file") != args.end();

  if (replace_state_file)
    TEST_COUT << "Replace state file ON" << std::endl;
  else
    TEST_COUT << "Replace state file OFF" << std::endl;

  bpt::ptree pt;
  bpt::ini_parser::read_ini(FUNCTIONAL_TESTS_DIR"/configure.ini", pt);

  std::string tag;
  std::size_t test_counter = 0;
  bool good_run = true;

  while(++test_counter < 1000) {
    /* EXTRACTION INI DATA */
    std::stringstream ss;
    ss << "Test" << test_counter;
    if (pt.find(ss.str()) == pt.not_found()) continue;

    tag = ss.str() + ".use_cuda";
    std::string ini_use_cuda = def::support::utilities::strtolower(pt.get<std::string>(tag));
    ASSERT_TRUE(ini_use_cuda == "no" || ini_use_cuda == "yes") << "Wrong value for " << tag;

    TEST_COUT << "_____________________________________________________________" << std::endl;

    if (ini_use_cuda == "yes" && use_cuda == false) {
      TEST_COUT << "Skip the [Test" << test_counter << "] : CUDA IS NOT AVAILABLE" << std::endl;
      continue;
    }

    tag = ss.str() + ".use_double_precision";
    std::string ini_use_double_precision = def::support::utilities::strtolower(pt.get<std::string>(tag));
    ASSERT_TRUE(ini_use_double_precision == "no" || ini_use_double_precision == "yes") << "Wrong value for " << tag;

    if (ini_use_double_precision == "yes" && use_double_precision == false) {
      TEST_COUT << "Skip the [Test" << test_counter << "] : DOUBLE PRECISION IS NOT AVAILABLE" << std::endl;
      continue;
    }

    tag = ss.str() + ".tolerance";
    float ini_tolerance = pt.get<float>(tag);
    ASSERT_GE(ini_tolerance, 0.0) << "Wrong value for " << tag;
    EXPECT_LT(ini_tolerance, 10.0) << "Tolerance is too big " << tag;

    tag = ss.str() + ".path";
    std::string path = pt.get<std::string>(tag);

    tag = ss.str() + ".exec";
    std::string exec = pt.get<std::string>(tag);

    tag = ss.str() + ".state-compare";
    std::string state_compare = pt.get<std::string>(tag);

    auto compare_binary = FUNCTIONAL_TESTS_DIR"/" + path + "/" + state_compare;

    if (!file_exist(compare_binary)) {
      TEST_COUT << "Skip the [Test" << test_counter << "] : Missing serialization binary file ["  << compare_binary + "]" << std::endl;
      continue;
    }

    /* CREATING TMP DIRECTORY */
    boost::system::error_code ec;
    bfs::path test_root(bfs::unique_path(bfs::temp_directory_path() / "%%%%-%%%%-%%%%"));
    ASSERT_TRUE(bfs::create_directory(test_root, ec)) << "Failed creating " << test_root << ": " << ec.message() << std::endl;
    std::string output_state_file = test_root.string() + "/output-state-file.bin";
    auto working_directory = FUNCTIONAL_TESTS_DIR"/" + path;
    auto log_file = test_root.string() + "/log.txt";
    auto err_file = test_root.string() + "/error.txt";

    std::stringstream cmdline;
    cmdline << exec
            << " --output-state-file=" << output_state_file
            << " --output-dir=" << test_root.string();

    TEST_COUT << "Running functional tests ["  << ss.str() + "]" << std::endl;
    TEST_COUT << "Tolerance ["  << ini_tolerance << "]" << std::endl;
    TEST_COUT << "Log file: [" << log_file << "]" << std::endl;
    TEST_COUT << "Error file: [" << err_file << "]" << std::endl;
    TEST_COUT << "Temporary Output Directory: [" << test_root.string() << "]" << std::endl;
    TEST_COUT << "Exec: [" << cmdline.str() << "]" << std::endl;

    try {

      try {
        /* Running deformetrica in a fork */
        deformetrica_run(cmdline.str(), working_directory, log_file, err_file);
      } catch(...) {
        TEST_COUT << "Current test has crashed" << std::endl;
        good_run = false;
        continue;
      }

      if (!file_exist(output_state_file)) {
        TEST_COUT << "Current test has not produced the state-file ["  << output_state_file + "]" << std::endl;
        good_run = false;
        continue;
      }

      def::utils::DeformationState df1, df2;

      try {
        df1.load(compare_binary);
      } catch(std::exception& err) {
        TEST_COUT << err.what() << " - Error while reading the [" << compare_binary << "] file" << std::endl;
        throw err;
      }

      try {
        df2.load(output_state_file);
      } catch(std::exception& err) {
        TEST_COUT << err.what() << " - Error while reading the [" << output_state_file << "] file" << std::endl;
        throw err;
      }

      if (df1.compare(df2, ini_tolerance))
        TEST_COUT << "Functional tests PASSED" << std::endl;
      else
        TEST_COUT << "Functional tests NOT PASSED" << std::endl;

    } catch (...){
      good_run = false;
    }

    if (replace_state_file) {
      try {
        bfs::copy_file(output_state_file, compare_binary, bfs::copy_option::overwrite_if_exists);
        TEST_COUT << "binary state file replaced" << std::endl;
      } catch(...) {
        TEST_COUT << "binary state file can't be replaced" << std::endl;
        good_run = false;
      }
    }

  }

  TEST_COUT << "_____________________________________________________________" << std::endl;
  ASSERT_TRUE(good_run);
}

}
}
