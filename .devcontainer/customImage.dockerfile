FROM mcr.microsoft.com/vscode/devcontainers/universal:linux

ENV VCPKG_ROOT /vcpkg
ENV CMAKE_TOOLCHAIN_FILE ${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake

RUN sudo git clone https://github.com/microsoft/vcpkg.git && \
	.${VCPKG_ROOT}/bootstrap-vcpkg.sh -disableMetrics && \
	sudo .${VCPKG_ROOT}/vcpkg integrate install && \
	rm -rf ${VCPKG_ROOT}/buildtrees/* && \
	rm -rf ${VCPKG_ROOT}/downloads/*
