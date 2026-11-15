"""Exercise real uv environment creation; stop before network/native installation."""
import os
from pathlib import Path
import shutil
import shlex
import subprocess
import sys
import zipfile

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
    outer = tmp_path / "other-environment"
    subprocess.run([uv, "venv", str(outer)], env=env, check=True, capture_output=True)
    env["UV_PYTHON"] = str(outer / "bin" / "python")
    wheel = tmp_path / "bootstrap_probe-0.0.0-py3-none-any.whl"
    files = {
        "bootstrap_probe.py": "VALUE = 1\n",
        "bootstrap_probe-0.0.0.dist-info/METADATA": "Metadata-Version: 2.1\nName: bootstrap-probe\nVersion: 0.0.0\n",
        "bootstrap_probe-0.0.0.dist-info/WHEEL": "Wheel-Version: 1.0\nRoot-Is-Purelib: true\nTag: py3-none-any\n",
    }
    files["bootstrap_probe-0.0.0.dist-info/RECORD"] = "".join(name + ",,\n" for name in files)
    with zipfile.ZipFile(wheel, "w") as archive:
        for name, text in files.items():
            archive.writestr(name, text)
    if existing:
        subprocess.run([uv, "venv", str(venv)], env=env, check=True, capture_output=True)
        (venv / "user-marker").write_text("preserve me")
        config = (venv / "pyvenv.cfg").read_bytes()
    tools = tmp_path / "tools"
    tools.mkdir()
    wrapper, driver = tools / "uv", tools / "uv-driver.py"
    # Model the dependency install only. uv itself handles the real venv semantics;
    # editable installation stops deliberately before any downloads or compilation.
    driver.write_text(f'''
import os
from pathlib import Path
import shutil
import subprocess
import sys
args = sys.argv[1:]
if args[0] == "venv":
    raise SystemExit(subprocess.call([{uv!r}, *args]))
target = ["--python", args[args.index("--python")+1]] if "--python" in args else []
subprocess.run([{uv!r}, "pip", "install", "--no-index", *target, {str(wheel)!r}],
               check=True, capture_output=True)
subprocess.run([str(Path(os.environ["VIRTUAL_ENV"]) / "bin" / "python"),
                "-c", "import bootstrap_probe"], check=True)
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
    wrapper.write_text("#!/bin/sh\nexec " + shlex.quote(sys.executable) + " " +
                       shlex.quote(str(driver)) + ' "$@"\n')
    wrapper.chmod(0o755)
    env["PATH"] = str(tools) + os.pathsep + "/usr/bin:/bin"
    result = subprocess.run(["/bin/bash", "setup.sh"], cwd=project, env=env,
                            capture_output=True, text=True)
    assert result.returncode == 64, result.stdout + result.stderr
    assert (venv / "bin" / "python").is_file()
    assert not list(outer.glob("lib/python*/site-packages/bootstrap_probe*"))
    if existing:
        assert (venv / "user-marker").read_text() == "preserve me"
        assert (venv / "pyvenv.cfg").read_bytes() == config
