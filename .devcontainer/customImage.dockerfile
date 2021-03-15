FROM mcr.microsoft.com/vscode/devcontainers/universal:linux

ENV VCPKG_ROOT /vcpkg
ENV CMAKE_TOOLCHAIN_FILE ${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake

RUN pwd
RUN git clone https://github.com/microsoft/vcpkg.git
