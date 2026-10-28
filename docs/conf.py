# Configuration for the Sphinx documentation of GWFrames.  Build it from the
# top of the repository with
#
#     uv sync --no-editable --group docs
#     uv run --no-sync sphinx-build -W --keep-going -b html docs docs/_build/html
#
# The API pages are generated from the docstrings of the installed package,
# so the package must be built and installed first.

import json
import os
import sys
from importlib.metadata import version as _version
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

# Run the notebooks in the interpreter running Sphinx, where GWFrames is
# installed, rather than in whatever kernel named "python3" the user may have
# registered.  Directories in JUPYTER_PATH are searched before the user's.
_jupyter_path = Path(__file__).parent / "_build" / "jupyter"
_kernel_dir = _jupyter_path / "kernels" / "python3"
_kernel_dir.mkdir(parents=True, exist_ok=True)
(_kernel_dir / "kernel.json").write_text(
    json.dumps(
        {
            "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
            "display_name": "Python 3",
            "language": "python",
        }
    )
)
os.environ["JUPYTER_PATH"] = os.pathsep.join(filter(None, [str(_jupyter_path), os.environ.get("JUPYTER_PATH")]))

project = "GWFrames"
author = "Michael Boyle"
copyright = "2013–2026, Michael Boyle"
release = _version("GWFrames")
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "myst_nb",
    "pybind11_signatures",
]

exclude_patterns = ["_build", "jupyter_execute", "**.ipynb_checkpoints"]
templates_path = []
default_role = "py:obj"
nitpicky = False

# Markdown and notebooks
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist"]
myst_heading_anchors = 3
nb_execution_mode = "cache"
nb_execution_timeout = 300
nb_execution_raise_on_error = True
nb_execution_show_tb = True
# The C++ code reports its progress on stderr, which would clutter the pages
nb_output_stderr = "remove"

# API documentation
autodoc_class_signature = "separated"
autodoc_default_options = {"exclude-members": "__init__"}
autodoc_typehints = "signature"
python_use_unqualified_type_names = True
python_maximum_signature_line_length = 100
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False
napoleon_preprocess_types = True
napoleon_type_aliases = {
    "Quaternion": "~GWFrames.Quaternions.Quaternion",
    "Waveform": "~GWFrames.Waveform",
    "PNWaveform": "~GWFrames.PNWaveform",
    "WaveformAtAPointFT": "~GWFrames.WaveformAtAPointFT",
    "Matrix": "~GWFrames.Matrix",
    "DataGrid": "~GWFrames.DataGrid",
    "Modes": "~GWFrames.Modes",
    "ScriFunctor": "~GWFrames.ScriFunctor",
    "SliceOfScriGrid": "~GWFrames.SliceOfScriGrid",
    "SliceOfScriModes": "~GWFrames.SliceOfScriModes",
    "SliceModes": "~GWFrames.SliceModes",
    "Scri": "~GWFrames.Scri",
    "SuperMomenta": "~GWFrames.SuperMomenta",
    "array_like": ":term:`array_like <numpy:array_like>`",
    "numpy.ndarray": "~numpy.ndarray",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "h5py": ("https://docs.h5py.org/en/stable", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

# HTML output
html_theme = "furo"
html_title = f"GWFrames {release}"
html_static_path = []
html_theme_options = {
    "source_repository": "https://github.com/moble/GWFrames",
    "source_branch": "main",
    "source_directory": "docs/",
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/moble/GWFrames",
            "html": (
                '<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">'
                '<path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17'
                ".55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68"
                "-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64"
                "-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27"
                " 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15"
                " 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013"
                ' 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>'
            ),
            "class": "",
        },
    ],
}
