# Generates a data file of N particles of Species S and box size L
# From the terminal, enter N as the first argument and S as the second
# The fourth argument is the name of the file to be saved to. Careful!!

using CSV
using DataFrames
using Random

if isinteractive() 
	print("Enter desired number of particles: ")
	const N = parse(Int64, readline())
	print("Enter desired particle species: ")
  const S = readline()
	print("Enter desired particle box size: ")
	const L = parse(Float64,readline())
	print("Enter save-file: ")
  const savefile = readline()
else
	const N = parse(Int64,ARGS[1]) 
	const S = ARGS[2]
	const L = parse(Float64, ARGS[3])
	const savefile = ARGS[4]
end	
println("N: ", N, " S: ", S, " L: ", L)
const gridres = Float64(1/10000)
const Vmax    = 10		                          # Hardcoded :(	
println("gridres: ", gridres)

spc = [S for i in 1:N]
pos = rand(0:gridres:L,(N,3)) #rand(0:gridres:Vmax,(N,3)))
vel = rand(0:gridres:Vmax,(N,3))
Data = hcat(spc,pos,vel)
df = DataFrame(Data, :auto)
println(df)


touch(savefile)
open(savefile, "w") do file 
	CSV.write(file, df, delim='	')
end




