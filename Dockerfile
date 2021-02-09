FROM jupyter/datascience-notebook:r-4.0.3

RUN mkdir /GEOImporter
COPY src/*.R /GEOImporter/