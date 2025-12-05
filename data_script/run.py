import argparse
from dataProcessor.dataProcessor import *
from dataProcessor.graphProcessor import *

# Parse arguments.
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input", type=str, help="input file name", required=True)
parser.add_argument("-o", "--output", type=str, help="output file name", required=False)
parser.add_argument("-t", "--timelimit", type=int, help="timelimit", required=True)
parser.add_argument("-m", "--mode", type=str, help="mode: time/gap/node", required=True)
args = parser.parse_args()

# Process input table.
table = getProcessedTable(args.input)
instances = list(table.keys())
algorithms = list(table[instances[0]].keys())

# Delete problematic instances.
to_remove = set()
for instance in instances:
    for alg in algorithms:
        if alg not in table[instance]:
            print("Deleting instance: " + instance)
            to_remove.add(instance)

for instance in to_remove:
    instances.remove(instance)

algorithms = ["ILS", "ILS+SRI"]
setGapInTables(table, instances, algorithms, args.timelimit)

# algorithms = [
#     "$\\textsc{pr}^*$+$D1$",
#     "$\\textsc{pr}^*$+$D2$",
#     "$\\textsc{pr}^*$+$D2$+\\textsc{set}",
# ]
# algorithms = [
#     "\\textsc{pr}+$\\mathcal{W}_{HS}$",
#     "\\textsc{pr}+$\\mathcal{W}_{OF}$",
#     "\\textsc{pr}+$\\mathcal{W}_{HS}$+\\textsc{sep}",
#     "\\textsc{pr}+$\\mathcal{W}_{OF}$+\\textsc{sep}",
# ]
mode = args.mode
if mode == "avg":
    printAverages(table, instances, algorithms, args.timelimit)
elif mode == "time-gap":
    createTimeAndGapGraph(args.output, table, instances, algorithms, args.timelimit)
elif mode == "time" or mode == "gap" or mode == "node":
    createGraph(args.output, table, instances, algorithms, mode, args.timelimit)
else:
    print("Invalid mode!")
