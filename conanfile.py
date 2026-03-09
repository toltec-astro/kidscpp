#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "conan>=2.10",
#   "tula-cmake",
# ]
#
# [tool.uv.sources]
# # Local editable source (sibling directory, monorepo layout):
# tula-cmake = { path = "../tula/tula_cmake" }
# # Once pyproject.toml is merged to the remote, switch to:
# # tula-cmake = { git = "https://github.com/toltec-astro/tula_cmake.git" }
# ///
#
# Self-bootstrapping Conan recipe.
#
# When loaded by Conan (normal use):
#   conan install . --profile=$(tula-cmake profiles-dir)/linux-gcc14-debug ...
#
# When run directly with uv (zero-install bootstrap):
#   uv run conanfile.py install . --profile=... -o "&:Eigen3=conan" ...
#   ./conanfile.py install . ...          # chmod +x first
#
# uv installs conan + tula-cmake into an isolated venv, then forwards all
# arguments to `conan`.  The `# /// script` block is plain comments to Python
# so Conan never sees it.

try:
    from tula_cmake import TulaConan          # pip install -e /path/to/tula_cmake
except ImportError:
    import os, sys
    from pathlib import Path
    _root = Path(__file__).parent
    for _d in [Path(os.environ.get("TULA_CMAKE_DIR", "")),
               _root.parent / "tula" / "tula_cmake"]:
        if _d and (_d / "tula_conan.py").exists():
            sys.path.insert(0, str(_d)); break
    else:
        raise ImportError(
            "tula_cmake not found.\n"
            "  Install: pip install -e /path/to/tula_cmake\n"
            "  Or set:  TULA_CMAKE_DIR=/path/to/tula_cmake\n"
            "  Or bootstrap: uv run conanfile.py install ."
        )
    from tula_conan import TulaConan


class KidsCppRecipe(TulaConan):
    """
    kidscpp v3: KIDs data processing for TolTEC — timestream solving.

    Heavy port of refs/kidscpp, trimmed to timestream solving only.
    Sweep finding/fitting is now handled by Python.

    Dependencies (via tula v3 conan-centric build):
      - Eigen3   (conan)  - linear algebra
      - logging  (conan)  - spdlog + fmt
      - Yaml     (conan)  - yaml-cpp config
      - Enum     (cpm)    - meta_enum + bitmask
      - Grppi    (cpm)    - parallel patterns
      - Ceres    (conan)  - non-linear optimization (sweep model calibration)
      - NetCDF   (system) - NetCDF file I/O (system-installed 4.9.2)
      - NetCDFCXX4(system)- NetCDF C++ bindings (system-installed 4.3.1)
      - Clipp    (conan)  - CLI parsing

    Usage (from /workspaces/cpp/kidscpp):
        uv run conanfile.py install . \\
          --profile:build=$(tula-cmake profiles-dir)/linux-gcc14-debug \\
          --profile:host=$(tula-cmake profiles-dir)/linux-gcc14-debug \\
          -o "&:Eigen3=conan" -o "&:logging=conan" -o "&:Yaml=conan" \\
          -o "&:Enum=cpm" -o "&:Grppi=cpm" -o "&:Ceres=conan" \\
          -o "&:Clipp=conan" -o "&:NetCDF=system" -o "&:NetCDFCXX4=system" \\
          --build=missing --output-folder=build/gcc14-debug
    """
    pass


if __name__ == "__main__":
    import subprocess, sys
    raise SystemExit(subprocess.run(["conan"] + sys.argv[1:]).returncode)
