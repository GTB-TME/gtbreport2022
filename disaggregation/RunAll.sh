#!/bin/bash

# NOTE currently being changed! 
# a script to run all of disaggregation analyses from the top
# NB this will be slow!
# NB this should be run from the directory where it lives


# age --------------
# prevalence survey: 01prev_survey.R
echo 'runnning prevalence survey analysis...'
cd R                           # change directory
R -q -e "rmarkdown::render(\"01prev_survey.R\",output_dir=\"../html\")"  # run & report
cd ..
echo '...finished!'


echo 'runnning child incidence analysis...'
cd R/02children/                           # change directory
R --slave <021paedmechrun.R
cd ../..
echo '...finished!'

# SLOW!
# adult incidence splits: 03inc_split.R
echo 'runnning incidence splits analysis...'
cd R                           # change directory
R -q -e "rmarkdown::render(\"03inc_split.R\",output_dir=\"../html\")" # run & report
cd ..
echo '...finished!'

#renaming to avoid windows issues
mv output/incsplits/plots/pass1/COM1.pdf  output/incsplits/plots/pass1/COMisISO1.pdf

# incidence splits finishing: 04inc_split_finish.R
echo 'runnning final incidence splits analysis...'
cd R                           # change directory
R -q -e "rmarkdown::render(\"04inc_split_finish.R\",output_dir=\"../html\")"
cd ..
echo '...finished!'

echo 'runnning child mortality data preparation (if not done)...'
cd R/02children/mortinput                           # change directory
R --slave <0220dataprep.R
cd ../../..
echo '...finished!'

# child mortality: 023cfrmortality.R
echo 'runnning child mortality analysis...'
cd R/02children/                           # change directory
R --slave <022cfrmortality.R
cd ../..
echo '...finished!'

# # mortality --------------
# final mortality splits
echo 'runnning mortality splits...'
cd R                           # change directory
R -q -e "rmarkdown::render(\"05mortsplit.R\",output_dir=\"../html\")"
cd ../..
echo '...finished!'


# finished!!!!!!
echo 'ALL DONE!'

