FROM julia:latest

RUN apt-get update && apt-get install make gfortran -y

ENV USER pluto
ENV USER_HOME_DIR /home/${USER}
ENV JULIA_DEPOT_PATH ${USER_HOME_DIR}/.julia
ENV NOTEBOOK_DIR ${USER_HOME_DIR}/notebooks
ENV JULIA_NUM_THREADS 100

RUN useradd -m -d ${USER_HOME_DIR} ${USER} \
    && mkdir ${NOTEBOOK_DIR}

COPY AIOMFAC ${USER_HOME_DIR}/AIOMFAC/  
COPY install.sh ${USER_HOME_DIR}/
COPY prestartup.jl ${USER_HOME_DIR}/
COPY startup.jl ${USER_HOME_DIR}/
COPY notebooks/webapp.jl ${NOTEBOOK_DIR}/
COPY notebooks/Project.toml ${NOTEBOOK_DIR}/
COPY notebooks/Manifest.toml ${NOTEBOOK_DIR}/

RUN mkdir -p ${NOTEBOOK_DIR}/Inputfiles \
	&& mkdir -p ${NOTEBOOK_DIR}/Outputfiles

RUN julia ${USER_HOME_DIR}/prestartup.jl \
    && chown -R ${USER} ${USER_HOME_DIR}

USER ${USER}

RUN mkdir -p ${NOTEBOOK_DIR}/Inputfiles \
	&& mkdir -p ${NOTEBOOK_DIR}/Outputfiles

WORKDIR ${USER_HOME_DIR}
RUN sh ${USER_HOME_DIR}/install.sh

EXPOSE 1234
VOLUME ${NOTEBOOK_DIR}
WORKDIR ${NOTEBOOK_DIR}

CMD [ "julia", "/home/pluto/startup.jl" ]
