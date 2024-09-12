FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get --fix-missing update -y && apt-get install -y \
    liberfa-dev \
    linux-tools-generic \
    pkg-config \
    python3-pip \
&& python3 -m pip install meson ninja pyuvdata

WORKDIR /work

COPY . /work/radiointerferometryc99

RUN cd /work/radiointerferometryc99 \
&& meson setup /work/radiointerferometryc99_build \
&& cd /work/radiointerferometryc99_build \
&& ninja install