FROM jupyter/datascience-notebook:r-4.0.3

USER root
RUN mkdir /GEOImporter
COPY src/*.R /GEOImporter/
COPY lib/*.tar.gz /GEOImporter/
USER $NB_USER