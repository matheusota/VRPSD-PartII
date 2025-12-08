from runner_script.helper import *
import argparse

# command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-m", "--mode", type=str, required=True, help="mode")
parser.add_argument(
    "-t", "--timelimit", type=int, required=True, help="time limit in seconds"
)
parser.add_argument("-o", "--output", type=str, required=True, help="output file")
parser.add_argument("-b", "--basic", action="store_true", help="Run basic variant.")
args = parser.parse_args()
timelimit = args.timelimit
tableName = args.output
mode = args.mode

if args.basic:
    models = [
        " --pr --ps --basic ",
        " --pr --sri --flow --basic ",
        " --pr --sri --flow --basic --prScenarioOptimal ",
    ]
else:
    models = [
        " --pr --ps ",
        " --pr --sri --flow ",
        " --pr --sri --flow  --prScenarioOptimal ",
    ]

# set parameters according to mode
if mode == "jabali2":
    folders = ["./data/Jabali/k2"]
elif mode == "jabali3":
    folders = ["./data/Jabali/k3"]
elif mode == "jabali4":
    folders = ["./data/Jabali/k4"]
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
