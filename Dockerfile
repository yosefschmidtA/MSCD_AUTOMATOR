# Start with a clean Ubuntu 22.04
FROM ubuntu:22.04

# Prevent interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# --- MPI Configuration (Allow running as Root) ---
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# 1. Install Dependencies
# Enable 32-bit architecture (Crucial for psrm.x binary)
RUN dpkg --add-architecture i386
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \
    libc6:i386 \
    libstdc++6:i386 \
    python3 \
    python3-numpy \
    python3-pandas \
    python3-matplotlib \
    python3-scipy \
    git \
    && rm -rf /var/lib/apt/lists/*

# 2. Setup Application Directory
WORKDIR /app
COPY . /app

# 3. Set Executable Permissions
# Ensure scripts and binaries are executable
RUN chmod +x run_all.sh docker_entrypoint.sh arquivos/psrm.x arquivos/phsh1

# 4. Compile Required Binaries
WORKDIR /app/arquivos

# Compile poconv (C++)
RUN rm -f poconv *.o && \
    cd poconv_src && \
    mpic++ -Wno-write-strings -o ../poconv poconv.cpp polation.cpp potentia.cpp userinfo.cpp userutil.cpp

# Compile phsh0 (Fortran)
RUN cd phsh0_src && \
    rm -f ../phsh0 *.o && \
    make FC=gfortran && \
    cp phsh0 ..

# 5. Create Mount Point for User Data
VOLUME /io

# 6. Set Entrypoint
WORKDIR /app
CMD ["./docker_entrypoint.sh"]
