FROM satijalab/seurat:3.1.4 

COPY DESCRIPTION .

# RUN Rscript -e "install.packages('remotes')" \
#             -e "remotes::install_cran('tinytex'); tinytex::install_tinytex()" \
#             -e "remotes::install_cran('lintr')" \
#             -e "remotes::install_deps(dependencies = TRUE)"

RUN Rscript -e "install.packages('lintr')"