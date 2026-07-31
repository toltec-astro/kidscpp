"""Spack package for TolTEC KIDs timestream processing."""

from spack.package import (
    depends_on,
    on_package_attributes,
    run_after,
    version,
    working_dir,
)
from spack.util.executable import Executable
from spack_repo.builtin.build_systems.cmake import CMakePackage


class Kidscpp(CMakePackage):
    """Build the Kidscpp library and its real timestream reader/solver tests."""

    homepage = "https://github.com/toltec-astro/kidscpp"

    version("3.1.0")

    depends_on("cmake@3.25:", type="build")
    depends_on("cxx", type="build")
    depends_on("tula-cmake@3.2.0", type="build")
    depends_on(
        "tula@3.1.0+ecsv+netcdf+enum+grppi+openmp",
        type=("build", "link"),
    )
    depends_on("googletest@1.14:~shared", type=("build", "test"))

    def cmake_args(self) -> list[str]:
        """Enable native tests only when Spack requests package testing."""
        return [self.define("KIDS_BUILD_TESTS", self.run_tests)]

    @run_after("build")
    @on_package_attributes(run_tests=True)
    def check(self) -> None:
        """Run the Kidscpp behavior tests before installation."""
        with working_dir(self.build_directory):
            ctest = Executable("ctest")
            listing = ctest("-N", output=str)
            if "Total Tests: 0" in listing:
                raise RuntimeError("Kidscpp configured without tests")
            ctest("--output-on-failure")
