from conans import ConanFile, CMake, tools

class L2APConan(ConanFile):
  name = "l2ap"
  version = "0.2.1"
  license = "MIT"
  author = "Jeremy Iverson (jiverson002@csbsju.edu)"
  url = "https://github.com/jiverson002/l2ap"
  homepage = "www.davidanastasiu.net/software/l2ap"
  description = "L2-norm All-Pairs similarity search."
  topics = ("hpc", "sfr", "apss")
  settings = "os", "compiler", "build_type", "arch"
  options = {
    "shared": [True, False],
    "fPIC": [True, False],
    "visibility": ["hidden", "default"]
  }
  default_options = {
    "shared": False,
    "fPIC": True,
    "visibility": "hidden"
  }
  exports = ["LICENSE"]
  exports_sources = "GKlib/*", "*.c", "*.h", "CMakeLists.txt", "L2APConfig.cmake.in"

  def build(self):
    cmake = CMake(self)
    #cmake.verbose = True
    cmake.configure()
    cmake.build()
    #cmake.test()
    cmake.install()

  def package(self):
    self.copy("*.h", dst="include", src="include")
    self.copy("*l2ap.lib", dst="lib", keep_path=False)
    self.copy("*.dll", dst="bin", keep_path=False)
    self.copy("*.so", dst="lib", keep_path=False)
    self.copy("*.dylib", dst="lib", keep_path=False)
    self.copy("*.a", dst="lib", keep_path=False)
