FROM nvidia/cuda:11.4.3-devel-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG OPENCV_VERSION=4.5.3

# Install system dependencies and Python packages
RUN apt-get update && apt-get upgrade -y &&\
    apt-get install -y \
        python3-pip \
        build-essential \
        cmake \
        git \
        wget \
        unzip \
        yasm \
        pkg-config \
        libswscale-dev \
        libtbb2 \
        libtbb-dev \
        libjpeg-dev \
        libpng-dev \
        libtiff-dev \
        libavformat-dev \
        libpq-dev \
        libxine2-dev \
        libglew-dev \
        libtiff5-dev \
        zlib1g-dev \
        libjpeg-dev \
        libavcodec-dev \
        libavformat-dev \
        libavutil-dev \
        libpostproc-dev \
        libswscale-dev \
        libeigen3-dev \
        libtbb-dev \
        libgtk2.0-dev \
        pkg-config \
        python3-dev \
        python3-numpy \
        gosu \
    && rm -rf /var/lib/apt/lists/*

# Install Jupyter Lab
RUN pip3 install jupyterlab

# Install Python packages: cupy, numpy, matplotlib
RUN pip3 install cupy numpy matplotlib

# Set up OpenCV
RUN cd /opt/ &&\
    wget https://github.com/opencv/opencv/archive/$OPENCV_VERSION.zip &&\
    unzip $OPENCV_VERSION.zip &&\
    rm $OPENCV_VERSION.zip &&\
    wget https://github.com/opencv/opencv_contrib/archive/$OPENCV_VERSION.zip &&\
    unzip ${OPENCV_VERSION}.zip &&\
    rm ${OPENCV_VERSION}.zip &&\
    mkdir /opt/opencv-${OPENCV_VERSION}/build && cd /opt/opencv-${OPENCV_VERSION}/build &&\
    cmake \
        -DOPENCV_EXTRA_MODULES_PATH=/opt/opencv_contrib-${OPENCV_VERSION}/modules \
        -DWITH_CUDA=ON \
        -DCUDA_ARCH_BIN=7.5,8.0,8.6 \
        -DCMAKE_BUILD_TYPE=RELEASE \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        .. &&\
    make -j"$(nproc)" && \
    make install && \
    ldconfig &&\
    rm -rf /opt/opencv-${OPENCV_VERSION} && rm -rf /opt/opencv_contrib-${OPENCV_VERSION}

# Copy entrypoint script
COPY createuser.sh /usr/bin/createuser.sh
RUN chmod +x /usr/bin/createuser.sh

# Expose port for Jupyter Lab
EXPOSE XXXX

# Use entrypoint script
ENTRYPOINT ["/usr/bin/createuser.sh"]
CMD ["jupyter-lab", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
