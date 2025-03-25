#!/bin/bash
R --slave < ../R/5update_fr.R
R --slave < ../R/6incidencecomparison.R
R --slave < ../R/7incidenceoutput.R
