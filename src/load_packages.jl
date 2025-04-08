# Load packages
using XLSX
using DataFrames
using StatsBase:mean,std
using Statistics
using LinearAlgebra
using Accessors # to mutate NamedTuple
using Distributions
using Combinatorics
using FreqTables
using DataFramesMeta
using CategoricalArrays
using GLM
using CairoMakie
using Random
import TidierPlots as gg
# set plot_log to false
gg.TidierPlots_set("plot_log",false)