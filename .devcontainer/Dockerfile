#-------------------------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See https://go.microsoft.com/fwlink/?linkid=2090316 for license information.
#-------------------------------------------------------------------------------------------------------------

FROM mcr.microsoft.com/vscode/devcontainers/base:0-buster

# This Dockerfile's base image has a non-root user with sudo access. Use the "remoteUser"
# property in devcontainer.json to use it. On Linux, the container user's GID/UIDs
# will be updated to match your local UID/GID (when using the dockerFile property).
# See https://aka.ms/vscode-remote/containers/non-root-user for details.
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

ENV VCPKG_ROOT /usr/vcpkg

# Configure apt and install packages
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    #
    # Install C++ tools and Git
    && apt-get -y install build-essential gdb cmake cppcheck git zip unzip \
    #
    # [Optional] Update UID/GID if needed
    && if [ "$USER_GID" != "1000" ] || [ "$USER_UID" != "1000" ]; then \
    groupmod --gid $USER_GID $USERNAME \
    && usermod --uid $USER_UID --gid $USER_GID $USERNAME \
    && chown -R $USER_UID:$USER_GID /home/$USERNAME; \
    fi \
    #
    # Clean up
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Install vcpkg dependency manager
RUN git clone https://github.com/Microsoft/vcpkg ${VCPKG_ROOT} \
    && cd ${VCPKG_ROOT} \
    # This is (for now) the way vcpkg handles versioning: freezing the cloned repo at a specific commit!
    # See https://github.com/microsoft/vcpkg/issues/1681.
    && git reset --hard ec6fe06e8da05a8157dc8581fa96b36b571c1bd5 \
    && ./bootstrap-vcpkg.sh

ENV PATH "${VCPKG_ROOT}:${PATH}"

# Install C++ dependencies using vcpkg
RUN vcpkg integrate install \
    && vcpkg install gtest