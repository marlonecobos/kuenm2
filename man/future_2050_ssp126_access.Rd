% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{future_2050_ssp126_access}
\alias{future_2050_ssp126_access}
\title{SpatRaster Representing Future Conditions (2041-2060, SSP126, GCM: ACCESS-CM2)}
\format{
A \code{SpatRaster} object.
}
\value{
No return value. Used with function \code{\link[terra]{rast}} to
bring raster variables to analysis.
}
\description{
A raster layer containing bioclimatic variables representing future climatic
conditions (2041-2060) based on the ACCESS-CM2 General Circulation Model
under the SSP126 scenario. The variables were obtained at a 10 arc-minute
resolution and masked using the \code{m} region provided in the package. Data
sourced from WorldClim: \url{https://worldclim.org/data/cmip6/cmip6climate.html}
}
\examples{
future_2050_ssp126_access <- terra::rast(system.file("extdata",
                                    "wc2.1_10m_bioc_ACCESS-CM2_ssp126_2041-2060.tif",
                                     package = "kuenm2"))
terra::plot(future_2050_ssp126_access)
}
