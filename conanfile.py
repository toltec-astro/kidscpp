import os
import subprocess
import sys
from pathlib import Path


def _find_tula_cmake(project_root: Path = Path(__file__).parent) -> Path:
    """Locate tula_cmake locally or fetch from GitHub (sparse clone).

    Search order:
      1. Sibling directory: <project_root>/../tula/tula_cmake/
      2. TULA_CMAKE_DIR environment variable
      3. Sparse-checkout cache: <project_root>/.tula_bootstrap/tula/tula_cmake/
         (auto-populated via git clone --sparse on first use)

    Override the remote with:
      TULA_GIT_REPO  (default: https://github.com/toltec-astro/tula.git)
      TULA_GIT_TAG   (default: main)
    """
    # 1. Sibling directory (monorepo / side-by-side clone)
    sibling = project_root.parent / "tula" / "tula_cmake"
    if (sibling / "tula_conan.py").exists():
        return sibling
    # 2. Environment variable
    if (env := os.environ.get("TULA_CMAKE_DIR")):
        p = Path(env)
        if (p / "tula_conan.py").exists():
            return p
    # 3. Sparse-checkout cache
    cache = project_root / ".tula_bootstrap" / "tula"
    tula_cmake = cache / "tula_cmake"
    if not (tula_cmake / "tula_conan.py").exists():
        repo = os.environ.get("TULA_GIT_REPO", "https://github.com/toltec-astro/tula.git")
        tag  = os.environ.get("TULA_GIT_TAG",  "main")
        print(f"[tula] fetching tula_cmake ({tag}) from {repo}")
        cache.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "clone", "--depth=1", "--filter=blob:none",
                        "--sparse", "--branch", tag, repo, str(cache)], check=True)
        subprocess.run(["git", "-C", str(cache), "sparse-checkout", "set",
                        "tula_cmake"], check=True)
    return tula_cmake


sys.path.insert(0, str(_find_tula_cmake()))
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
          --profile:build=../tula/tula_cmake/profiles/linux-gcc14-debug \\
          --profile:host=../tula/tula_cmake/profiles/linux-gcc14-debug \\
          -o "&:Eigen3=conan" -o "&:logging=conan" -o "&:Yaml=conan" \\
          -o "&:Enum=cpm" -o "&:Grppi=cpm" -o "&:Ceres=conan" \\
          -o "&:Clipp=conan" -o "&:NetCDF=system" -o "&:NetCDFCXX4=system" \\
          --build=missing --output-folder=build/gcc14-debug
    """
    pass
