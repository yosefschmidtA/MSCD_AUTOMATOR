#!/bin/bash
#
# Script to install all dependencies for MSCD_AUTOMATOR
# on a Red Hat-based system (Fedora, RHEL, CentOS).
#
# Run this script using: bash install_fedora.sh
#

echo "--- Installing MSCD_AUTOMATOR dependencies (Fedora/RHEL) ---"

# 1. Install the "Development Tools" group
echo "Installing Development Tools (g++, make, gfortran...)"
sudo dnf groupinstall -y "Development Tools"

# 2. Install MPI, Python/Numpy, and 32-bit libraries
echo "Installing MPI, Python data libraries, and 32-bit libraries..."
# 'glibc.i686' is the 32-bit C library (our libc6:i386 equivalent)
sudo dnf install -y openmpi openmpi-devel python3 python3-numpy python3-pandas python3-matplotlib python3-scipy glibc.i686

echo "--------------------------------------------------"
echo "✅ Dependency installation complete!"
echo "You can now run the main script: ./run_all.sh"
echo "--------------------------------------------------"
