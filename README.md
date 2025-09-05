# RIsearch
RIsearch consists of two programs, RIsearch1 and RIsearch2

* RIsearch: a tool for large-scale RNA–RNA, RNA-DNA, and DNA-DNA interaction prediction
* RIsearch2 is a program for fast RNA-RNA interaction prediction using a seed-and-extend approach

## Installation

First, make sure that you have the following programs installed:

* gcc
* make
* CMake
* libpcre3-devel

To build and install

	./readybuild.sh
	./configure
	make
	make install

## Executing RIsearch
Here we provide examples on how to run RIsearch2 and RIsearch1 with default parameters.

See the manuals for RIsearch2 and RIsearch1 for the full list of options and examples.

### RIsearch1
To run RIsearch1 in default settings use the following command:

	RIsearch -q query.fa -t target.fa

Alternatively, single sequences can be given directly on commandline with
-Q acgu -T acgu

Note: the files query.fa and target.fa may contain several sequences, 
RIsearch will scan all vs. all


### RIsearch2
First, generate the index structure for the target sequence(s) in the file target.fa
and store them in the file target.suf. 

	RIsearch2 -c target.fa -o target.suf

To run RIsearch2 in default settings use the following command:

	RIsearch2 -q query.fa -i target.suf

Note: the files query.fa and target.fa may contain several sequences, 
RIsearch will scan all vs. all

## RIsearch2 siRNA off-target prediction

siRNA off-target Discovery Pipeline.

	./scripts/pipeline.py

Script that divides the input target sequence in 100kb pieces, runs RNAplfold
with them and packs the results into a binary format to be able to pass them
into the pipeline.

	 ./src/run_RNAplfold_and_pack_results.py

## RIsearch manuals

Please read the full manuals for RIsearch1 and RIsearch2, included in the
documentation folder.

## Copyright

Copyright 2015-2025 by the contributors (see AUTHORS and RIsearch1/README files)

RIsearchis released under the GNU General Public License version 3. Note that
libdivsufsort is packaged along with RIsearch. Libdivsufsort comes with its own
authors and copyright under libdivsufsort-2.0.1/{AUTHORS,COPYING}

GNU GENERAL PUBLIC LICENSE

This is a free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License, either version 3 of the License, or
(at your option) any later version. You should have received a copy of the GNU General Public License
along with RIsearch, see file COPYING. If not, see <http://www.gnu.org/licenses/>.

This software is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citations

If you use RIsearch in your publication please cite:

**RIsearch2: suffix array-based large-scale prediction of RNA-RNA interactions and siRNA off-targets**
Alkan F, Wenzel A, Palasca O, Kerpedjiev P, Rudebeck AF, Stadler PF, Hofacker IL, Gorodkin J *Nucleic Acids Res*. 2017 May 5;45(8):e60

**RIsearch: fast RNA-RNA interaction search using a simplified nearest-neighbor energy model**
Wenzel A, Akbasli E, Gorodkin J. *Bioinformatics*. 2012 Nov 1;28(21):2738-46. Epub 2012 Aug 24

## Contact

In case of problems or bug reports, please contact: <software@rth.dk>

