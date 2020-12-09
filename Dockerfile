FROM satijalab/seurat:3.1.4 

RUN Rscript -e "install.packages(c('lintr', 'ggsignif'))"