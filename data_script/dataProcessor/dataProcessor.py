import pandas as pd


# Given a string of the format filepath/instance_name.svrp we return the instance name only.
def getInstanceName(S):
    S = S.split("/")
    aux = S[-1]
    aux = aux[:-6]

    return aux


# Given an input file, get the corresponding table data structure.
def getProcessedTable(inputFile):
    f = open(inputFile, "r")
    table = {}
    indexMap = {}

    i = 0
    for line in f:
        # First line.
        if i == 0:
            line = line.split(",")
            for j, field in enumerate(line):
                indexMap[field.strip()] = j

        else:
            # Get run info.
            line = line.replace(",  ,", ",")
            line = line.strip()
            line = line.split(",")
            algorithm = line[indexMap["Algorithm"]].strip()
            instance = getInstanceName(line[indexMap["Instance"]]).strip()

            # Create new hashmap/row if necessary.
            if instance not in table:
                table[instance] = {}

            # Add entry to table.
            table[instance][algorithm] = {}

            for key in indexMap:
                aux = line[indexMap[key]]
                aux = aux.strip()

                if key not in ["Algorithm", "Instance"]:
                    table[instance][algorithm][key] = float(aux)

        i += 1

    return table


def setGapInTables(table, instances, algorithms, timelimit):
    # For each instance, we get the best UB and set the gap accordingly.
    for instance in instances:
        bestUB = 1e9

        for alg in algorithms:
            if alg in table[instance]:
                bestUB = min(bestUB, table[instance][alg]["Model Solution Cost"])

        for alg in algorithms:
            if (
                table[instance][alg]["Node Count"] <= 1
                and table[instance][alg]["Time (s)"] <= timelimit
            ):
                table[instance][alg]["Root Lower Bound"] = bestUB

            table[instance][alg]["Gap"] = (
                (bestUB - table[instance][alg]["Root Lower Bound"]) / bestUB
            ) * 100.0

            table[instance][alg]["Final Gap"] = (
                (bestUB - table[instance][alg]["Lower Bound"]) / bestUB
            ) * 100.0


def printAverages(table, instances, algorithms, timeLimit):
    # set pandas dataframe
    df = pd.DataFrame(
        columns=[
            "Instance",
            "Algorithm",
            "Time (s)",
            "Root Gap (\\%)",
            "Node Count",
            "CVRPSEP Cuts",
            "Partial Route Cuts",
            "Partial Route Time (s)",
            "Root Time (s)",
            "Final Gap",
        ]
    )

    i = 0
    for instance in instances:
        for alg in algorithms:
            if alg not in table[instance]:
                continue

            df.loc[i] = [
                instance,
                alg,
                table[instance][alg]["Time (s)"],
                table[instance][alg]["Gap"],
                table[instance][alg]["Node Count"],
                table[instance][alg]["CVRPSEP Cuts"],
                table[instance][alg]["Partial Route Cuts"],
                table[instance][alg]["Partial Route Time (s)"],
                table[instance][alg]["Root Time (s)"],
                table[instance][alg]["Final Gap"],
            ]
            i += 1

    print(df)
    print("\nAverages\n")
    print(
        df.groupby(["Algorithm"])[
            [
                "Time (s)",
                "Root Gap (\\%)",
                "Node Count",
                "CVRPSEP Cuts",
                "Partial Route Cuts",
                "Partial Route Time (s)",
                "Root Time (s)",
            ]
        ].mean()
    )

    print("\nCount\n")
    print("Total: " + str(len(instances)))
    df2 = df.loc[df["Time (s)"] <= timeLimit]
    print(df2.groupby(["Algorithm"]).size())

    print("\nFinal Gap\n")
    df3 = df.loc[df["Time (s)"] >= timeLimit - 1]
    print(df3.groupby(["Algorithm"])["Final Gap"].mean())
