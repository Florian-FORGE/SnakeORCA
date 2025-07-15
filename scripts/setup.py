from setuptools import setup, find_packages

setup(
    name="Stage-IA3D",
    version="0.1",
    packages=find_packages(where="scripts"),
    package_dir={"": "scripts"},
)
