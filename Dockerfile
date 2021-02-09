FROM jupyter/datascience-notebook:r-4.0.3

USER root
RUN mkdir /GEOImporter
COPY src/*.R /GEOImporter/
USER $NB_USER