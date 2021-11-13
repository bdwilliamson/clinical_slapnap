#!/usr/local/bin/Rscript
# clinical slapnap simulation 2: improving power using slapnap -----------------
# load required functions and packages -----------------------------------------
library("here")
library("tidyverse")
library("data.table")
library("argparse")
library("sievePH") # for Lunn and McNeil (1995) method

source(here("R", "00_sim_utils.R"))

# get command-line args --------------------------------------------------------
# these specify: the number of total replicates, the number of replicates per job
parser <- ArgumentParser()
parser$add_argument("--nreps-total", default = 1000, type = "double", help = "the total number of replicates")
parser$add_argument("--nreps-per-job", default = 1000, type = "double", help = "the total number of replicates per job")
parser$add_argument("--output-dir", default = here::here("R_output", "simulation_1"), help = "the output directory")
args <- parser$parse_args()

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  job_id <- 1
} else {
  job_id <-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

# define "important" AA positions in gp120 (from AMP sieve analysis plan)

# read in the gamma values
