#!/usr/bin/env bash

git clone https://github.com/rgcca-factory/RGCCA.git ~/bin/RGCCA
cd ~/bin/RGCCA
git checkout nucleipark
git clone https://github.com/BrainAndSpineInstitute/rgcca_Rpackage.git ~/bin/rgccaLauncher/
cd ~/bin/rgccaLauncher
git checkout nucleipark3