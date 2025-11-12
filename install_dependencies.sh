#!/bin/bash
#
# Script to install all dependencies for MSCD_AUTOMATOR
# on a Debian/Ubuntu-based system.
#
# Run this script using: bash install_dependencies.sh
#

echo "--- Starting MSCD_AUTOMATOR dependency installation ---"

# 1. Update the package list
echo "Updating package lists (apt update)..."
sudo apt-get update

# 2. Add 32-bit architecture support (required for psrm.x)
echo "Configuring 32-bit architecture (i386)..."
sudo dpkg --add-architecture i386
sudo apt-get update -q

# 3. Install essential compilers (C, C++, Fortran, MPI)
echo "Installing compilers (build-essential, gfortran, openmpi)..."
sudo apt-get install -y build-essential openmpi-bin libopenmpi-dev gfortran

# 4. Install 32-bit compatibility libraries (for psrm.x)
echo "Installing 32-bit compatibility libraries (libc6:i386)..."
sudo apt-get install -y libc6:i386

# 5. Install Python 3 and required libraries
echo "Installing Python 3 and data science libraries..."
sudo apt-get install -y python3 python3-numpy python3-pandas python3-matplotlib python3-scipy

echo "--------------------------------------------------"
echo "✅ Dependency installation complete!"
echo "You can now run the main script: ./run_all.sh"
echo "--------------------------------------------------"
