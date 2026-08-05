"""Spack package for TolTEC KIDs timestream processing."""

from spack.package import (
    depends_on,
    on_package_attributes,
    run_after,
    variant,
    version,
    working_dir,
)
from spack.util.executable import Executable
from spack_repo.builtin.build_systems.cmake import CMakePackage


class Kidscpp(CMakePackage):
    """Build the Kidscpp library and its real timestream reader/solver tests."""

    homepage = "https://github.com/toltec-astro/kidscpp"
    git = "https://github.com/toltec-astro/kidscpp.git"

    version("3.1.0", commit="8d80a80ae6683db8391973de02fd3457b90a9a3c")

    variant(
        "openmp",
        default=True,
        description="Enable OpenMP in the transitive Tula performance layer",
    )

    depends_on("cmake@3.25:", type="build")
    depends_on("cxx", type="build")
    depends_on("tula-cmake@3.2.0", type="build")
    depends_on(
        "tula@3.1.0+ecsv+netcdf+enum+grppi+openmp",
        when="+openmp",
        type=("build", "link"),
    )
    depends_on(
        "tula@3.1.0+ecsv+netcdf+enum+grppi~openmp",
        when="~openmp",
        type=("build", "link"),
    )
    depends_on("googletest@1.14:~shared", type=("build", "test"))

    def cmake_args(self) -> list[str]:
        """Enable native tests only when Spack requests package testing."""
        return [
            self.define("KIDS_BUILD_TESTS", self.run_tests),
            self.define("KIDSCPP_PACKAGE_SPEC", str(self.spec)),
            self.define("KIDSCPP_DAG_HASH", self.spec.dag_hash()),
        ]

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
