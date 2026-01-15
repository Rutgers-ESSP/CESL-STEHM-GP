# CESL-STEHM-GP: Spatio-Temporal Empirical Hierarchical Modeling of sea-level data with Gaussian Processes 

README file last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2026-01-15

## Citation

Version 1.0 of this code was released on 30 January 2016 to accompany

	Kopp, R. E., A. C. Kemp, K. Bittermann, B. P. Horton, 	J. P. Donnelly,
	W. R. Gehrels, C. C. Hay, J. X. Mitrovica, 	E. D. Morrow, and S. Rahmstorf
	(2016). Temperature-driven 	global sea-level variability in the Common
	Era. Proceedings 	of the National Academy of Sciences 113: E1434-E1441.
	doi: 10.1073/pnas.1517056113.
	
Please cite this source when using this code.

Version 2.0 of this code was updated to accompany

    Kemp, A. C., Wright, A. J., Edwards, R. J., Barnett, R. L., Brain, M. J.,
    Kopp, R. E., Cahill, N., Horton, B. P., Charman, D. J., Hawkes, A. D.,
    Hill, T. D, & van de Plassche, O. (2018). . Quaternary
    Science Reviews, 201, 89-110.
	
Version 2.1 of this code was updated to accompany

    Walker, J. S., R. E. Kopp, C. M. Little, and B. P. Horton (2022).
	Timing of emergence of modern rates of sea-level rise by 1863.
	Nature Communications. https://doi.org/10.1038/s41467-022-28564-6
	
Note that this code base is not actively maintained. 
New users are encouraged to instead employ [PaleoSTeHM](https://github.com/radical-collaboration/PaleoSTeHM).

## Overview

This code is intended to be used to generate Gaussian process-based, spatio-temporal empirical hierarchical models of sea-level data with both vertical and temporal uncertainty. It assumes uncertainties are normally distributed, and uses the noisy-input Gaussian process method of McHutchon and Ramussen (2011) to translate temporal uncertainties into comparable vertical uncertainties.    

This code requires MATLAB with the Statistics, Optimization, and Global Optimization Toolboxes. It is known to run with MATLAB 2015b.

The MFILES directory contains the core functions of this code. The scripts-KoppEtAl2016 directory contains the scripts using these core functions to generate the analyses described in Kopp et al. (2016). The scripts-KempEtAl2018 directory contains the scripts associated with Kemp et al. (2018). The scripts-WalkerEtAl2022 directory contains the scripts associated with Walker et al. (2022). The scripts-skeleton directory contains a template that can be used for regional applications.

The master script for these applications is runAnalyzeCESL, which calls various subsidiary scripts. Particularly attention should be paid to runSetupHoloceneCovariance, which defines the covariance functions with characterize the different priors, including the hyperparameter bounds used by the optimization code. This is the most intellectually challenging part of the exercise, as it is effectively the model specification step. Attention should also be paid to runImportCESLDataSets, which imports the data files. 

## Current status

This code is no longer being actively developed. Please see the [PaleoSTeHM project](https://github.com/radical-collaboration/PaleoSTeHM/) for new developments.

## Acknowledgements

The development of version 1.0 of this code was supported in part by the US National Science Foundation (grants ARC-1203415), the National Oceanic and Atmospheric Administration (grant NA14OAR4170085) and the New Jersey Sea Grant Consortium. Subsequent developmentGP was supported by the National Science Foundation (awards 1458904, 1702587, 1804999, and 2002437).

----

    Copyright (C) 2022 by Robert E. Kopp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
