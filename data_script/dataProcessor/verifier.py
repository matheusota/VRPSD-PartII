# Check if optimal solutions match.
def checkCosts(table, instances, algorithms):
    for instance in instances:
        opt_cost = 0
        for alg in algorithms:
            if table[instance][alg]["Solved"]:
                opt_cost = table[instance][alg]["Model Solution Cost"]

        for alg in algorithms:
            if (
                table[instance][alg]["Solved"]
                and table[instance][alg]["Model Solution Cost"] <= opt_cost - 0.01
            ):
                print("Error in instance: " + instance)
                print(
                    "Cost of "
                    + alg
                    + " = "
                    + str(table[instance][alg]["Model Solution Cost"])
                )
                print("Optimal cost = " + str(opt_cost))


# Delete instances that were not solved.
def deleteInstances(table, instances, algorithms):
    to_remove = set()
    for instance in instances:
        for alg in algorithms:
            if alg not in table[instance]:
                print(instance + " not solved by " + alg)
                print(table[instance])
                print("Deleting instance: " + instance)
                to_remove.add(instance)
                exit(1)

    for instance in to_remove:
        instances.remove(instance)
