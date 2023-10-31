import sys

from setuptools import setup, find_packages

# Some Python installations don't add the current directory to path.
if '' not in sys.path:
    sys.path.insert(0, '')

setup(
      name="radqm9pipeline",
      version="0.1.0",
      description="...",
      author="Matthew Avaylon",
      author_email='mavaylon@lbl.gov',
      platforms=["Mac OS and Linux"],
      license="",
      url="",
      packages=find_packages('src'),
      package_dir={'': 'src'},
      install_requires=[
            "tqdm",
            "ase",
            "numpy",
            "monty",
            "PyAstronomy"
            ],
      python_requires=">=3.9"
)