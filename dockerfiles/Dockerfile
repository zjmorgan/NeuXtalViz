ARG SOURCE_IMAGE=source

FROM --platform=amd64 mambaorg/micromamba:2.0-ubuntu20.04 AS source

USER root

COPY . /src
WORKDIR /src

RUN micromamba env create --yes -f environment.yml

ENTRYPOINT [ "micromamba", "run", "-n", "nxv" ]


# This is a workaround that allows the COPY --from location to be defined as a build argument.
# With this, we can reference the built source image from Harbor in our pipelines, while still
# defaulting to using the previous stage for local builds.
FROM $SOURCE_IMAGE AS source_image


FROM --platform=amd64 kitware/trame:py3.10-conda-glvnd-2025-04 AS run

COPY --from=source /src /src
WORKDIR /src

RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba
RUN conda env create --yes -f environment.yml
SHELL ["conda", "run", "-n", "nxv", "/bin/bash", "-c"]
ENTRYPOINT [ "conda", "run", "--no-capture-output", "-n", "nxv" ]

ENV VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow

RUN apt update && apt install -y libxrender1 nginx xvfb

RUN chmod og+rwX -R /var/lib/nginx
RUN chmod og+rwX -R /var/log/nginx
RUN chmod og+rwX -R /etc/nginx

COPY dockerfiles/nginx.conf.template /etc/nginx/nginx.conf.template
COPY dockerfiles/supervisord.conf /etc/supervisord.conf
COPY scripts/run_nginx.sh /
COPY scripts/run_trame.sh /

RUN pip install supervisor

RUN python -m trame.tools.www --client-type vue3 --output /app/www-content

# Get the static JS from nova-trame and add it to our client build
RUN mkdir -p /app/www-content/assets/js
RUN export NOVA_TRAME_LOCATION=$(pip show nova-trame | grep "Location" | awk '{ print $2 }')/nova/trame; cp $NOVA_TRAME_LOCATION/view/theme/assets/js/*.js /app/www-content/assets/js/

USER trame-user

ENV EP_PATH=/app
ENV NOVA_TRAME_DEFAULT_DEBOUNCE=200

CMD ["supervisord"]
