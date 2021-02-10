FROM jupyter/datascience-notebook:r-4.0.3

USER root
RUN mkdir /GEOImporter
RUN chown $NB_USER /GEOImporter
#USER $NB_USER

COPY src/*.R /GEOImporter/
COPY lib/*.tar.gz /GEOImporter/