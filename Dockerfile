# This runs an R script belonging to a local library and saves its plotly graphs as image files

FROM rocker/r-ver:4.4

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  apt-utils \
  ed \
  libnlopt-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/

RUN apt-get update && apt-get install -y  --no-install-recommends \
  software-properties-common \
  dirmngr \
  git-core \
  libcurl4-openssl-dev \
  libgit2-dev \
  libssh2-1-dev \
  libicu-dev \
  libpng-dev \
  libudunits2-dev \
  zlib1g-dev \
  libgdal-dev \
  libproj-dev \
  xml2 \
  openssl \
  python3  # needed for `save_image()`

# needed to install devtools
RUN apt-get update && apt-get install -y  --no-install-recommends \
  libssl-dev \
  libfontconfig1-dev \
  libxml2-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev

RUN mkdir /epiline/

COPY R/ /epiline/R/
COPY src/ /epiline/src/
COPY examples/ /epiline/examples/
COPY DESCRIPTION /epiline/
COPY NAMESPACE /epiline/
COPY tools/ /epiline/tools/
COPY examples/ /epiline/examples/

WORKDIR /epiline

# required by local library
# not using urls to .tar.gz files due to `{packages} are not available for package ‘EpiLine’`
RUN Rscript -e "install.packages('http://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.8.2.tar.gz', repos=NULL, type='source')"
RUN Rscript -e "install.packages('rstan', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('rstantools', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('data.table', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('plotly', repos='http://cran.rstudio.com/')"
RUN Rscript -e "install.packages('moments', repos='http://cran.rstudio.com/')"

# needed for `save_image()`
RUN Rscript -e "install.packages('reticulate', repos='http://cran.rstudio.com/')"

# # needed to use local library
RUN Rscript -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
RUN Rscript -e "devtools::load_all()"
RUN Rscript -e "devtools::install()"

# needed for `save_image()`
RUN Rscript -e "reticulate::install_miniconda()"
RUN Rscript -e "reticulate::conda_install('r-reticulate', 'python-kaleido')"
RUN Rscript -e "reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')"
RUN Rscript -e "reticulate::use_miniconda('r-reticulate')"

CMD Rscript /epiline/examples/linear_r_dist.R
