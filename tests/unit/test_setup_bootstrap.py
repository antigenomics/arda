"""Exercise real uv environment creation; stop before network/native installation."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

pytestmark = pytest.mark.skipif(os.name == "nt", reason="setup.sh is the POSIX bootstrap")


@pytest.mark.parametrize("existing", [False, True])
def test_setup_supplies_isolated_cmake_and_preserves_existing_environment(tmp_path, existing):
    uv = shutil.which("uv")
    if uv is None:
        pytest.skip("developer setup requires uv")
    project = tmp_path / "project"
    project.mkdir()
    shutil.copyfile(Path(__file__).parents[2] / "setup.sh", project / "setup.sh")
    venv = project / ".venv"
    env = dict(os.environ, UV_PYTHON=sys.executable, UV_OFFLINE="1",
               UV_CACHE_DIR=str(tmp_path / "uv-cache"), UV_NO_CONFIG="1")
    if existing:
        subprocess.run([uv, "venv", str(venv)], env=env, check=True, capture_output=True)
        (venv / "user-marker").write_text("preserve me")
        config = (venv / "pyvenv.cfg").read_bytes()
    tools = tmp_path / "tools"
    tools.mkdir()
    wrapper = tools / "uv"
    # Model the dependency install only. uv itself handles the real venv semantics;
    # editable installation stops deliberately before any downloads or compilation.
    wrapper.write_text(f'''#!{sys.executable}
import os
from pathlib import Path
import shutil
import subprocess
import sys
args = sys.argv[1:]
if args[0] == "venv":
    raise SystemExit(subprocess.call([{uv!r}, *args]))
if "-e" not in args:
    if any(arg == "cmake" or arg.startswith("cmake>") for arg in args):
        cmake = Path(os.environ["VIRTUAL_ENV"]) / "bin" / "cmake"
        cmake.write_text("#!/bin/sh\\necho 'cmake version 3.18.0'\\n")
        cmake.chmod(0o755)
    if any(arg == "airr" or arg.startswith("airr>") for arg in args):
        assert "--no-build-isolation" not in args
        (Path(os.environ["VIRTUAL_ENV"]) / "airr-isolated").touch()
    raise SystemExit(0)
cmake = shutil.which("cmake")
if cmake is None or Path(cmake).parent != Path(os.environ["VIRTUAL_ENV"]) / "bin":
    sys.exit("CMake missing from isolated build environment")
subprocess.run([cmake, "--version"], check=True)
if not (Path(os.environ["VIRTUAL_ENV"]) / "airr-isolated").exists():
    sys.exit("AIRR source build was not prepared with isolation")
raise SystemExit(64)  # reached editable build; do not install or fetch resources
''')
    wrapper.chmod(0o755)
    env["PATH"] = str(tools) + os.pathsep + "/usr/bin:/bin"
    result = subprocess.run(["/bin/bash", "setup.sh"], cwd=project, env=env,
                            capture_output=True, text=True)
    assert result.returncode == 64, result.stdout + result.stderr
    assert (venv / "bin" / "python").is_file()
    if existing:
        assert (venv / "user-marker").read_text() == "preserve me"
        assert (venv / "pyvenv.cfg").read_bytes() == config
