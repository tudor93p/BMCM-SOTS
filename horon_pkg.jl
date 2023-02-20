#using BMCMSOTS: PATH_SNAKE, MBPATH 

PATH_SNAKE = "/net/horon/scratch/pahomit/SnakeStates" 

ps = ["Constants","myPlots"]

for p in ps 

	Pkg.rem(p)

end 

for p in ps 

	Pkg.add(url=joinpath(PATH_SNAKE,p))

end 
