###################################################
#    This file is part of RPAWL.
#
#    RPAWL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RPAWL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RPAWL.  If not, see <http://www.gnu.org/licenses/>.
###################################################

#'PARALLEL ADAPTIVE WANG-LANDAU
#'
#'The package implements the Parallel Adaptive Wang-Landau algorithm on various
#'examples.  The provided demos allow to reproduce the figures of the article.
#'
#'\tabular{ll}{ Package: \tab PAWL\cr Type: \tab Package\cr Version: \tab
#'1.0\cr Date: \tab 2011-08-11\cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
#'yes\cr Depends: \tab mvtnorm\cr Suggests: \tab ggplot2\cr } The main function
#'is \code{pawl}. It takes algorithmic parameters in arguments (see the help of
#'the \code{pawl} function), as well a target distribution. Look at the demos
#'to learn how to specify a target distribution.
#'
#'@name PAWL-package
#'@aliases PAWL-package PAWL
#'@docType package
#'@author Luke Bornn <bornn@@stat.harvard.edu>, Pierre E. Jacob
#'<pierre.jacob.work@@gmail.com>
#'@keywords package
#'@examples
#'
#'  demo(discreteexample)
#'  demo(gaussianexample)
#'  demo(mixture2kexample)
#'
NULL

#'Image of ice floes
#'
#'This data represents a binary matrix, representing an image of ice floes.
#'
#'
#'@name IceFloe
#'@aliases IceFloe icefloe
#'@docType data
#'@format A matrix containing 40 rows and 40 columns
#'@source Banfield, J. and Raftery, A. (1992). Ice floe identification in
#'satellite images using mathematical morphology and clustering about principal
#'curves. Journal of the American Statistical Association, 87(417):7-16.
#'@keywords datasets
NULL

#'Pollution Data
#'
#'This data contains 1 response (mortality, normalized to have mean zero) along
#'with 15 pollution-related explanatory variables.
#'
#'
#'@name Pollution
#'@aliases Pollution pollution
#'@docType data
#'@format A matrix containing 60 rows and 16 columns
#'@source McDonald, G.C. and Schwing, R.C. (1973) 'Instabilities of regression
#'estimates relating air pollution to mortality', Technometrics, vol.15,
#'463-482.
#'@keywords datasets
NULL






