from conan import ConanFile
from conan.tools.cmake import CMake


class KidsCppRecipe(ConanFile):
    """Conan 2 recipe for the trimmed kidscpp timestream package."""

    name = "kidscpp"
    version = "3.1.0"
    description = "TolTEC KIDs timestream processing library"
    license = "BSD-3-Clause"
    url = "https://github.com/toltec-astro/kidscpp"
    package_type = "static-library"
    required_conan_version = ">=2.31"
    python_requires = "tula-cmake/3.1.0"
    python_requires_extend = "tula-cmake.TulaConan"
    settings = ()
    options = {}
    default_options = {
        "tula/*:logging": "conan",
        "tula/*:yaml_cpp": "conan",
        "tula/*:csv_parser": "cpm",
        "tula/*:netcdf_c": "system",
        "tula/*:netcdf_cxx4": "system",
        "tula/*:bitmask": "cpm",
        "tula/*:meta_enum": "cpm",
        "tula/*:perflibs": "system",
        "tula/*:eigen": "conan",
        "tula/*:grppi": "cpm",
    }
    tula_default_options = {
        "logging": "conan",
        "yaml_cpp": "conan",
        "csv_parser": "cpm",
        "netcdf_c": "system",
        "netcdf_cxx4": "system",
        "bitmask": "cpm",
        "meta_enum": "cpm",
        "perflibs": "system",
        "eigen": "conan",
        "grppi": "cpm",
    }
    tula_public_features = tuple(tula_default_options)
    exports_sources = "CMakeLists.txt", "include/*", "src/*", "tests/*"

    def requirements(self) -> None:
        super().requirements()
        self.requires(
            "tula/3.1.0",
            transitive_headers=True,
            transitive_libs=True,
        )

    def build_requirements(self) -> None:
        if not self.conf.get("tools.build:skip_test", default=False, check_type=bool):
            self.test_requires("gtest/1.17.0")

    def build(self) -> None:
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        if not self.conf.get("tools.build:skip_test", default=False, check_type=bool):
            cmake.ctest()

    def package(self) -> None:
        CMake(self).install()

    def package_info(self) -> None:
        self.cpp_info.set_property("cmake_file_name", "kidscpp")
        self.cpp_info.set_property("cmake_target_name", "kids::kids")
        self.cpp_info.libs = ["kids"]
