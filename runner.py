from runner_script.helper import *
import argparse

# command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-m", "--mode", type=str, required=True, help="mode")
parser.add_argument(
    "-t", "--timelimit", type=int, required=True, help="time limit in seconds"
)
parser.add_argument("-o", "--output", type=str, required=True, help="output file")
args = parser.parse_args()
timelimit = args.timelimit
tableName = args.output
mode = args.mode

models = [
    " --node --partial_route --sri --lagrangian --mip --flow ",
    " --node --partial_route --sri --lagrangian --mip --flow --dualset ",
]

# set parameters according to mode
if mode == "Jabali":
    folders = ["./data/Jabali/k3"]
elif mode == "parada":
    folders = ["./data/parada/scenarios"]
elif mode == "ccvrp":
    folders = ["./data/ccvrp"]
else:
    raise Exception("Unknown mode")

ks = [" --readk "]

# get total number of instances
total = 0
for folder in folders:
    all_instances = getAllInstances(folder)
    total += len(all_instances) * len(models) * len(ks)

# create empty csv table
createEmptyTable(tableName)

# run algorithms for every instance below folder
for folder in folders:
    all_instances = sorted(getAllInstances(folder))

    for instance in all_instances:
        runModels(models, ks, instance, timelimit, total, tableName)
