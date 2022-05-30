#

# Load Packages
using CSV
using DataFrames
using LinearAlgebra
using Plots
using Random

# Initialization

data = eachrow(CSV.read("example.txt", DataFrame, header=3))

# Find acceleration

# update loop
