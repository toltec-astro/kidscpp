try:
    from tula_cmake import TulaConan          # pip install -e /path/to/tula/tula_cmake
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
            "  Install: pip install -e /path/to/tula/tula_cmake\n"
            "  Or set:  TULA_CMAKE_DIR=/path/to/tula_cmake\n"
            "  Or run:  tula-cmake fetch --project-root ."
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
        conan install . \\
          --profile:build=$(tula-cmake profiles-dir)/linux-gcc14-debug \\
          --profile:host=$(tula-cmake profiles-dir)/linux-gcc14-debug \\
          -o "&:Eigen3=conan" -o "&:logging=conan" -o "&:Yaml=conan" \\
          -o "&:Enum=cpm" -o "&:Grppi=cpm" -o "&:Ceres=conan" \\
          -o "&:Clipp=conan" -o "&:NetCDF=system" -o "&:NetCDFCXX4=system" \\
          --build=missing --output-folder=build/gcc14-debug
    """
    pass
