import argparse
from dataProcessor.dataProcessor import *
from dataProcessor.graphProcessor import *
from dataProcessor.verifier import *

# Parse arguments.
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-t", "--timelimit", type=int, help="timelimit", required=True)
parser.add_argument(
    "-m", "--mode", type=str, help="mode: avg/paper/appendix", required=True
)
parser.add_argument("-d", "--date", type=str, help="date string", required=True)
args = parser.parse_args()
mode = args.mode
timelimit = args.timelimit
date = args.date

# Table file names.
ccvrp_input = "./" + date + "/" + date + "_ccvrp.csv"
jabali_input = "./" + date + "/" + date + "_jabali.csv"
ccvrp_input_basic = "./" + date + "/" + date + "_ccvrp_basic.csv"
jabali_input_basic = "./" + date + "/" + date + "_jabali_basic.csv"

# Process tables.
ccvrp_table = getProcessedTable(ccvrp_input)
jabali_table = getProcessedTable(jabali_input)
ccvrp_table_basic = getProcessedTable(ccvrp_input_basic)
jabali_table_basic = getProcessedTable(jabali_input_basic)
ccvrp_instances = list(ccvrp_table.keys())
jabali_instances = list(jabali_table.keys())
algorithms = list(ccvrp_table[ccvrp_instances[0]].keys())
algorithms = ["ILS", "ILS+SRI"]

# Verify if everything is ok.
# print("Deleting Normal instances")
# deleteInstances(ccvrp_table, ccvrp_instances, algorithms)
# deleteInstances(jabali_table, jabali_instances, algorithms)
# print("Deleting Basic instances")
# deleteInstances(ccvrp_table_basic, ccvrp_instances, algorithms)
# deleteInstances(jabali_table_basic, jabali_instances, algorithms)
# print("Verifying optimal costs for CCVRP.")
# checkCosts(ccvrp_table, ccvrp_instances, algorithms)
# checkCosts(ccvrp_table_basic, ccvrp_instances, algorithms)
# print("Verifying optimal costs for Jabali.")
# checkCosts(jabali_table, jabali_instances, algorithms)
# checkCosts(jabali_table_basic, jabali_instances, algorithms)

# Set gap based on best solution found.
setGapInTables(ccvrp_table, ccvrp_instances, algorithms)
setGapInTables(jabali_table, jabali_instances, algorithms)
setGapInTables(ccvrp_table_basic, ccvrp_instances, algorithms)
setGapInTables(jabali_table_basic, jabali_instances, algorithms)

# Get plots.
if mode == "avg":
    print("Averages for CCVRP:")
    printAverages(ccvrp_table, ccvrp_instances, algorithms, timelimit)
    print("Averages for CCVRP-Basic:")
    printAverages(ccvrp_table_basic, ccvrp_instances, algorithms, timelimit)
    print("Averages for Jabali:")
    printAverages(jabali_table, jabali_instances, algorithms, timelimit)
    print("Averages for Jabali-Basic:")
    printAverages(jabali_table_basic, jabali_instances, algorithms, timelimit)
elif mode == "paper1":
    print("Creating paper1 plot normal")
    createPaper1Graph(
        ccvrp_table,
        jabali_table,
        ccvrp_instances,
        jabali_instances,
        algorithms,
        timelimit,
        "paper_classical.png",
    )
    print("Creating paper1 plot basic")
    createPaper1Graph(
        ccvrp_table_basic,
        jabali_table_basic,
        ccvrp_instances,
        jabali_instances,
        algorithms,
        timelimit,
        "paper_classical_basic.png",
    )
elif mode == "paper2":
    print("Creating paper2 plot")
    createPaper2Graph(
        ccvrp_table,
        ccvrp_table_basic,
        jabali_table,
        jabali_table_basic,
        ccvrp_instances,
        jabali_instances,
        algorithms,
        timelimit,
        "paper2_classical.png",
    )
else:
    print("Invalid mode!")
