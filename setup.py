import os
import subprocess
import shutil
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

class BuildWithCMake(build_py):
    def run(self):
        
        # 1. Build C++ executables using CMake
        build_dir = os.path.abspath("build")
        os.makedirs(build_dir, exist_ok=True)
        

        # Allow passing custom CMake flags via CMAKE_ARGS env var
        cmake_args = os.environ.get("CMAKE_ARGS", "").split()

        # builds the executables using cmake
        subprocess.check_call(["cmake", ".."] + cmake_args, cwd=build_dir)
        subprocess.check_call(["cmake", "--build", "."], cwd=build_dir)

        # 2. Copy executables into the Python package
        output_bin_dir = os.path.join("GFoil", "bin")
        os.makedirs(output_bin_dir, exist_ok=True)
        shutil.copy(os.path.join(build_dir, "GFoil_fwd_codi"), output_bin_dir)
        shutil.copy(os.path.join(build_dir, "GFoil_AD"), output_bin_dir)
        
        # 3. Continue with the normal Python build
        super().run()


setup(
    name="GFoil",
    version="0.1.0",
    packages=find_packages(),
    cmdclass={
        'build_py': BuildWithCMake,
    },
    include_package_data=True,
    package_data={
        'GFoil': ['bin/*'],
    },
    install_requires=[],
    python_requires=">=3.7",
    author="Kieran Barry",
    description="Python interface to CFoil tools",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
)