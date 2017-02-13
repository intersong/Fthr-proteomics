from sys import argv
file = open(argv[1], "r")
for line in file:
   if line[0] == ">":
      print(line.split("embl-cds:")[1].split(" ")[0])
      print(line.split("GN=")[1].split(" ")[0])
