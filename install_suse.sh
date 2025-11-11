#!/bin/bash
#
# Script to install all dependencies for MSCD_AUTOMATOR
# on an openSUSE system.
#
# Run this script using: bash install_suse.sh
#

echo "--- Installing MSCD_AUTOMATOR dependencies (openSUSE) ---"

# 1. Refresh repositories
echo "Refreshing repositories..."
sudo zypper refresh

# 2. Install the C/C++ development "patterns"
echo "Installing C/C++ Development Patterns (devel_basis, devel_C_C++)..."
sudo zypper install -t pattern -y devel_basis devel_C_C++

# 3. Install Fortran, MPI, Python/Numpy, and 32-bit libraries
echo "Installing Fortran (gcc-fortran), MPI, Python/Numpy, and 32-bit libraries..."
# 'glibc-32bit' is the 32-bit C library (our libc6:i386 equivalent)
# 'gcc-fortran' provides the gfortran compiler
sudo zypper install -y gcc-fortran openmpi-devel python3-numpy glibc-32bit

echo "--------------------------------------------------"
echo "✅ Dependency installation complete!"
echo "You can now run the main script: ./run_all.sh"
echo "--------------------------------------------------"
