FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN dpkg --add-architecture i386 && apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \
    libc6:i386 \
    libstdc++6:i386 \
    python3 \
    python3-pip \
    python3-numpy \
    python3-pandas \
    python3-matplotlib \
    python3-scipy \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app

RUN pip3 install --no-cache-dir Flask gunicorn

RUN cd arquivos && \
    rm -f poconv *.o && \
    cd poconv_src && \
    mpic++ -Wno-write-strings -o ../poconv poconv.cpp polation.cpp potentia.cpp userinfo.cpp userutil.cpp

RUN cd arquivos/phsh0_src && \
    rm -f ../phsh0 *.o && \
    make FC=gfortran && \
    cp phsh0 ..

RUN chmod +x run_half.sh arquivos/psrm.x arquivos/phsh1 arquivos/poconv arquivos/phsh0

ENV LD_LIBRARY_PATH=/app/arquivos:$LD_LIBRARY_PATH

EXPOSE 5000

CMD ["gunicorn", "--workers", "1", "--timeout", "120", "--bind", "0.0.0.0:5000", "app:app"]
