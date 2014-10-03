Code files 
==========


## Code in this repository.

The final code versions and their relations are as follows:

1. `code/pdef_dynamics_2.0.c` : was started in multiple instances to run the simulations over the entire gradient of 1 to 100 predator species. outputs single files for each instance. 
2. `code/pdef_analysis.r` : reads inthe output files of 1., calculates different metrics (e.g. the final functional diversity) and merges everything into the object `webstats`. Generates result figure and fits the statistical models. The linear model results are displayed in the figure. 

## License

```
    Project: Predator Diversity and Ecosystem Functioning

    Copyright (C) 2014 Christian Guill & Florian D. Schneider

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
      
    Contact the authors: [Christian Guill](C.P.Guill@uva.nl), [Florian D. Schneider](florian.schneider@univ-montp2.fr)
```