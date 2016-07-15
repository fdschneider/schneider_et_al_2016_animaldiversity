Code release for the paper 'Animal diversity and ecosystem functioning in dynamic food webs'
===========================================================================================

by Florian D. Schneider, Ulrich Brose, Björn C. Rall & Christian Guill

This repository contains simulation code and code for statistical analysis of a research project. It applies allometric models of animal foraging and network structure to simulate the ecosystem-level stocks and rates of animal and plant biomass in relation to the number of animal species present in the ecosystem. 

## Code in this repository.

1. `code/pdef_dynamics_2.4.c` & `code/pdef_dynamics_1.1.h`: Original simulation code. This was compiled and run in multiple instances for simulations over the entire gradient of 10 to 100 animal species with a single file output file for each instance found in `data/pdef_2.4/`. Additional instances were run for sensitivity analyses and are found in the same location. 
2. `code/data.r` : Data collect and compilation. The file reads in all output files, corrects some columns, calculates additional columns and saves everything into a single output file `data/webstats.txt`. Seperate files are created for the sensititvity analyses. Column names correspond do parameters as used in simulation code and slightly differ from nomenclature in the article.
3. `code/analyse.r` : contains function `analyse()` which fits the statistical models describing the relationships between animal diversity and ecosystem function. Call using `analyse(webstats$S_c_fin, webstats$meanB_c, equation = TRUE, ylab = "")`
4. `code/sensitivity.r` : contains function `sensitivity()` for calling the sensitivity analysis for any parameter x, as well as function `sensitivity_plus()` for visualising the simulations of the alternative model runs. 
4. `code/analysis.r` : applies function `analyse()` to return the result figure 4 and 5 of the paper and an overview output table. applies functions in sensitivity.r to return sensitivity analysis of Supplementary Material. 

## Comments welcome

We encourage testing or reviewing this simulation code. Any comments, bug reports or mistakes should be reported in the [Github issues](https://github.com/fdschneider/schneider_et_al_2016_animaldiversity/issues).

## License

The MIT License (MIT)

Copyright (C) 2016 Christian Guill & Florian D. Schneider

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


