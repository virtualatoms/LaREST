import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

project = "LaREST"
author = "Ryan Reese, Alex Ganose, Charles Romain"
copyright = "2025, Ryan Reese, Alex Ganose, Charles Romain"  # noqa: A001
release = "1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
]

myst_enable_extensions = ["dollarmath"]

autodoc_mock_imports = ["rdkit", "numpy", "pandas", "tqdm", "censo", "tomllib"]

html_theme = "furo"
html_title = "LaREST"

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

exclude_patterns = ["_build"]
